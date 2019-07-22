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
        Spacer sp(k, w, spacing && *spacing ? parse_spacing(spacing, k): spvec_t(1, k));
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
        Spacer sp(k, w, spacing && *spacing ? parse_spacing(spacing, k): spvec_t(1, k));
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
}
