#include "pybind11/pybind11.h"
#include "pybind11/numpy.h"
#include "bonsai/encoder.h"
#include "hll/include/flat_hash_map/flat_hash_map.hpp"
namespace py = pybind11;
using namespace pybind11::literals;

using namespace bns;

template<typename Enc>
py::array_t<uint64_t> get_kmers(Enc &enc, const char *str, size_t n, size_t reservation, bool unique) {
    if(unique) {
        ska::flat_hash_set<uint64_t> kmerset;
        auto func = [&kmerset](auto x) {kmerset.insert(x);};
        if(n)
            enc.for_each(func, str, n);
        else
            enc.for_each(func, str);
        py::array_t<uint64_t> arr(kmerset.size());
        auto ptr = (uint64_t *)arr.request().ptr, p = ptr;
        for(const auto v: kmerset) *p++ = v;
        std::sort(ptr, p);
        return arr;
    } else {
        Py_ssize_t i = 0;
        py::array_t<uint64_t> arr(reservation);
        auto p = (uint64_t *)arr.request().ptr;
        auto func = [&](auto x) {
            if(i == arr.size()) {
                arr.resize({std::max(arr.size() << 1, Py_ssize_t(4))});
                p = (uint64_t *)arr.request().ptr;
            }
            p[i++] = x;
        };
        if(n) enc.for_each(func, str, n);
        else  enc.for_each(func, str);
        arr.resize({i});
        return arr;
    }
}

py::array_t<uint64_t> get_kmers(RollingHasher<uint64_t> &enc, const char *str, size_t n, size_t reservation, bool unique) {
    if(unique) {
        ska::flat_hash_set<uint64_t> kmerset;
        auto func = [&kmerset](auto x) {kmerset.insert(x);};
        if(n) {
            if(enc.canon_)
                enc.for_each_canon(func, str, n);
            else
                enc.for_each_uncanon(func, str, n);
        } else {
            enc.for_each_hash(func, str);
        }
        py::array_t<uint64_t> arr(kmerset.size());
        auto ptr = (uint64_t *)arr.request().ptr, p = ptr;
        for(const auto v: kmerset) *p++ = v;
        std::sort(ptr, p);
        return arr;
    } else {
        Py_ssize_t i = 0;
        py::array_t<uint64_t> arr(reservation);
        auto p = (uint64_t *)arr.request().ptr;
        auto func = [&](auto x) {
            if(i == arr.size()) {
                arr.resize({std::max(arr.size() << 1, Py_ssize_t(4))});
                p = (uint64_t *)arr.request().ptr;
            }
            p[i++] = x;
        };
        if(n) {
            if(enc.canon_) enc.for_each_canon(func, str, n);
            else           enc.for_each_uncanon(func, str, n);
        } else enc.for_each_hash(func, str);
        arr.resize({i});
        return arr;
    }
}


PYBIND11_MODULE(bns, m) {
    m.doc() = "Bonsai Python bindings";
    m.def("from_fasta_r", [](const char *path, int k, bool canon, size_t reserve, bool unique) {
        RollingHasher<uint64_t> enc(k, canon);
        return get_kmers(enc, path, 0, reserve, unique);
    }, py::return_value_policy::take_ownership, "return a numpy array of integers from the fasta.\nk: k-mer length\ncanon: Whether or not to canonicalize [True]\nreserve: initial buffer reserve [1024].\nunique: to make unique first [False]\n",
       "path"_a=nullptr, "k"_a=31, "canon"_a=true, "reserve"_a=1024, "unique"_a=false);
    m.def("seqlist", [](const char *path, int k,  const char *spacing, int w, bool canonicalize, size_t reservation, bool unique, bool rolling) {
        Spacer sp(k, w, spacing && *spacing ? parse_spacing(spacing, k): spvec_t(k - 1, 0));
        Encoder<> enc(sp, canonicalize);
        RollingHasher<uint64_t> rollingenc(k);
        gzFile fp = gzopen(path, "rb");
        if(!fp) throw std::runtime_error(std::string("Failed to open file at ") + path);
        py::list ret;
        kseq_t *ks = kseq_init(fp);
        while(kseq_read(ks) >= 0) {
            if(rolling)
                ret.append(get_kmers(rollingenc,  ks->seq.s, ks->seq.l, reservation, unique));
            else
                ret.append(get_kmers(enc, ks->seq.s, ks->seq.l, reservation, unique));
        }
        kseq_destroy(ks);
        gzclose(fp);
        return ret;
    }, py::return_value_policy::take_ownership, "return a numpy array of integers from the fasta.\nk: k-mer length\ncanon: Whether or not to canonicalize [True]\nreserve: initial buffer reserve [1024].\nunique: to make unique first [False]\n",
        "path"_a, "k"_a=31, "spacing"_a="", "w"_a=0, "canon"_a=true, "reserve"_a=1024, "unique"_a=false, "rolling"_a=false);
    m.def("from_fasta", [](const char *path, int k,  const char *spacing, int w, bool canonicalize, size_t reservation, bool unique) {
        Spacer sp(k, w, spacing && *spacing ? parse_spacing(spacing, k): spvec_t(k - 1, 0));
        Encoder<> enc(sp, canonicalize);
        return get_kmers(enc, path, 0, reservation, unique);
    }, py::return_value_policy::take_ownership, "return a numpy array of integers from the fasta.\nk: k-mer length\ncanon: Whether or not to canonicalize [True]\nreserve: initial buffer reserve [1024].\nunique: to make unique first [False]\n",
        "path"_a, "k"_a=31, "spacing"_a="", "w"_a=0, "canon"_a=true, "reserve"_a=1024, "unique"_a=false);
    m.def("from_str", [](const std::string &str, int k,  const char *spacing, int w, bool canon, size_t reserve) {
        Spacer sp(k, w, spacing && *spacing ? parse_spacing(spacing, k): spvec_t(k - 1, 0));
        Encoder<> enc(sp, canon);
        py::array_t<uint64_t> ret({reserve});
        auto bi = ret.request();
        ssize_t i = 0;
        uint64_t *p = static_cast<uint64_t *>(bi.ptr);
        enc.for_each([&](uint64_t x) {
            if(i == ret.size()) {
                ret.resize({std::max(ret.size() << 1, Py_ssize_t(4))});
                bi = ret.request();
                p = static_cast<uint64_t *>(bi.ptr);
            }
            p[i++] = x;
        }, str.data(), str.size());
        ret.resize({i});
        return ret;
    }, py::return_value_policy::take_ownership, "return a numpy array of integers from the fasta", "str"_a, "k"_a=31, "spacing"_a=nullptr, "w"_a=0, "canon"_a=true, "reserve"_a=1024);
    m.def("repack", [](py::list vallist, size_t ngenomes) -> py::array_t<float> {
        size_t nks = vallist.size();
        py::array_t<float> ret({(ngenomes * (ngenomes - 1)) >> 1, nks});
        float *retp = (float *)ret.request().ptr;
        assert(vallist.size() == nks);
        size_t kind = 0;
        for(py::handle ob: vallist) {
            auto vals = ob.cast<py::array_t<float>>();
            assert(vals.size() == ssize_t(nks * ((ngenomes * (ngenomes - 1)) >> 1)));
            float *valp = (float *)vals.request().ptr;
            for(size_t i = 0; i < ngenomes; ++i) {
                for(size_t j = i + 1; j < ngenomes; ++j) {
                    size_t offset = (i * (ngenomes * 2 - i - 1)) / 2 + j - (i + 1);
                    size_t o_offset = (ngenomes * (ngenomes - 1) / 2) * kind;
                    retp[offset + o_offset] = valp[offset];
                }
            }
            ++kind;
        }
        return ret;
    }, py::return_value_policy::take_ownership);
    m.def("seqdict_r", [](std::string path, int k) -> py::dict {
        py::dict ret;
        RollingHasher<uint64_t> enc(k);
        gzFile fp = gzopen(path.data(), "rb");
        if(!fp) throw std::runtime_error(std::string("Couldn't open file at ") + path);
        kseq_t *ks = kseq_init(fp);
        while(kseq_read(ks) >= 0) {
            py::array_t<uint64_t> arr({1u << 8});
            ssize_t i = 0;
            uint64_t *p = static_cast<uint64_t *>(arr.request().ptr);
            enc.for_each_hash([&](uint64_t x) {
                if(i == arr.size()) {
                    arr.resize({arr.size() << 1});
                    p = static_cast<uint64_t *>(arr.request().ptr);
                }
                p[i++] = x;
            }, path.data());
            arr.resize({i});
            ret[ks->name.s] = arr;
        }
        gzclose(fp);
        kseq_destroy(ks);
        return ret;
    }, py::arg("path"), py::arg("k") = 31, py::return_value_policy::take_ownership);
    m.def("seqdict", [](std::string path, int k, const char *spacing, int w) -> py::dict {
        py::dict ret;
        Spacer sp(k, w, spacing && *spacing ? parse_spacing(spacing, k): spvec_t(k - 1, 0));
        Encoder<> enc(sp);
        gzFile fp = gzopen(path.data(), "rb");
        if(!fp) throw std::runtime_error(std::string("Couldn't open file at ") + path);
        kseq_t *ks = kseq_init(fp);
        while(kseq_read(ks) >= 0) {
            py::array_t<uint64_t> arr({1u << 8});
            ssize_t i = 0;
            uint64_t *p = static_cast<uint64_t *>(arr.request().ptr);
            enc.for_each([&](uint64_t x) {
                if(i == arr.size()) {
                    arr.resize({arr.size() << 1});
                    p = static_cast<uint64_t *>(arr.request().ptr);
                }
                p[i++] = x;
            }, path.data());
            arr.resize({i});
            ret[ks->name.s] = arr;
        }
        gzclose(fp);
        kseq_destroy(ks);
        return ret;
    }, py::arg("path"), py::arg("k") = 31, py::arg("spacing") = "", py::arg("w") = 0, py::return_value_policy::take_ownership);
}
