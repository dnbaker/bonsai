#include "pybind11/pybind11.h"
#include "pybind11/numpy.h"
#include "bonsai/include/encoder.h"
namespace py = pybind11;
using namespace pybind11::literals;

using namespace bns;


PYBIND11_MODULE(bns, m) {
    m.doc() = "Bonsai Python bindings";
    m.def("from_fasta_r", [](const char *path, int k) {
        RollingHasher<uint64_t> enc(k);
        py::array_t<uint64_t> ret({1u << 16});
        ssize_t i = 0;
        uint64_t *p = static_cast<uint64_t *>(ret.request().ptr);
        enc.for_each_hash([&](uint64_t x) {
            if(i == ret.size()) {
                ret.resize({ret.size() << 1});
                p = static_cast<uint64_t *>(ret.request().ptr);
            }
            p[i++] = x;
        }, path);
        ret.resize({i});
        return ret;
    }, "return a numpy array of integers from the fasta", "path"_a=nullptr, "k"_a=31);
    m.def("from_fasta", [](const char *path, int k,  const char *spacing, int w) {
        Spacer sp(k, w, spacing && *spacing ? parse_spacing(spacing, k): spvec_t(k - 1, 0));
        Encoder<> enc(sp);
        py::array_t<uint64_t> ret({1u << 16});
        ssize_t i = 0;
        uint64_t *p = static_cast<uint64_t *>(ret.request().ptr);
        enc.for_each([&](uint64_t x) {
            if(i == ret.size()) {
                ret.resize({ret.size() << 1});
                p = static_cast<uint64_t *>(ret.request().ptr);
            }
            p[i++] = x;
        }, path);
        ret.resize({i});
        return ret;
    }, "return a numpy array of integers from the fasta", "path"_a=nullptr, "k"_a=31, "spacing"_a=nullptr, "w"_a=0);
    m.def("from_str", [](const std::string &str, int k,  const char *spacing, int w) {
        Spacer sp(k, w, spacing && *spacing ? parse_spacing(spacing, k): spvec_t(k - 1, 0));
        Encoder<> enc(sp);
        py::array_t<uint64_t> ret({1u << 16});
        ssize_t i = 0;
        uint64_t *p = static_cast<uint64_t *>(ret.request().ptr);
        enc.for_each([&](uint64_t x) {
            if(i == ret.size()) {
                ret.resize({ret.size() << 1});
                p = static_cast<uint64_t *>(ret.request().ptr);
            }
            p[i++] = x;
        }, str.data(), str.size());
        ret.resize({i});
        return ret;
    }, "return a numpy array of integers from the fasta", "str"_a, "k"_a=31, "spacing"_a=nullptr, "w"_a=0);
    m.def("repack", [](py::list vallist, size_t ngenomes) -> py::array_t<float> {
        size_t nks = vallist.size();
        py::array_t<float> ret({(ngenomes * (ngenomes - 1)) >> 1, nks});
        float *retp = (float *)ret.request().ptr;
        assert(vallist.size() == nks);
        size_t kind = 0;
        for(py::handle ob: vallist) {
            auto vals = ob.cast<py::array_t<float>>();
            assert(vals.size()  == nks * ((ngenomes * (ngenomes - 1)) >> 1));
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
}
