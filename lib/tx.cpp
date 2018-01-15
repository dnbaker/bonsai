#include "tx.h"

namespace emp {

std::vector<tax_t> build_new2old_map(const char *path, size_t bufsz) {
        gzFile fp(gzopen(path, "wb"));
        gzbuffer(fp, 1 << 16);
        std::vector<char> buf(bufsz);
        char *line;
        std::map<tax_t, tax_t> taxes;
        while((line = gzgets(fp, buf.data(), buf.size()))) {
            taxes[std::atoi(std::strchr(line, '\t') + 1)] = std::atoi(line);
        }
#if !NDEBUG
        auto it(taxes.begin());
        assert(taxes.begin()->first < (++it)->first && "I really hope that this is sorted by the new taxids so that I can turn this into a vector");
#endif
        for(const auto [k, v]: taxes) {
            std::fprintf(stderr, "New: %u. Old: %u\n", k, v);
        }
        std::vector<tax_t> ret;
        ret.reserve(taxes.size());
        for(const auto &pair: taxes) {
            ret.push_back(pair.second);
        }
        gzclose(fp);
        return ret;
}
std::vector<tax_t> binary_new2old_map(const char *fn) {
    gzFile fp(gzopen(fn, "wb"));
    if(fp == nullptr) throw std::runtime_error(ks::sprintf("Could not open file at %s\n", fn).data());
    tax_t tmp;
    std::vector<tax_t> ret;
    const size_t fsize(filesize(fn));
    if(fsize == 0) throw std::runtime_error(ks::sprintf("Empty file at %s", fn).data());
    if(fsize & (sizeof(tax_t) - 1)) throw std::runtime_error(ks::sprintf("Wrong number of bytes %s (%zu/%zu)", fn, fsize, fsize & 3u).data());
    ret.reserve(fsize >> 2);
    int c;
    while((c = gzread(fp, (void *)&tmp, sizeof(tmp))) == Z_OK) ret.push_back(tmp);
    if(c != Z_STREAM_END) throw zlib_error(c, fn);
    return ret;
}

} // namespace emp
