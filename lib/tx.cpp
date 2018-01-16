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


void bitmap_filler_helper(void *data_, long index, int tid) {
    bf_helper_t &data(*(bf_helper_t *)data_);
    ba::MMapTaxonomyBitmap &map(data.bm_);
    const tax_t taxid(data.taxes_[index]);
    u64 val;
    Encoder enc(data.sp_);
    gzFile fp(gzopen(data.paths_[index].data(), "rb"));
    kseq_t *ks(kseq_init(fp));
    while(kseq_read(ks) >= 0) {
        enc.assign(ks);
        while(enc.has_next_kmer())
            if((val = enc.next_minimizer()) != BF)
                map.set_kmer_ts(data.h_, val, taxid);
    }
    kseq_destroy(ks);
    gzclose(fp);
}
#if 0
struct bf_helper_t {
    const Spacer &sp_;
    const std::vector<std::string> &paths_;
    const khash_t(64) *h_;
};
#endif

} // namespace emp
