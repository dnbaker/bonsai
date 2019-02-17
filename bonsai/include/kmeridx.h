#ifndef BONSAI_KMER_IDX_H__
#define BONSAI_KMER_IDX_H__
#include "encoder.h"
#include "lazy/vector.h"
#include "flat_hash_map/flat_hash_map.hpp"

namespace bns {

template<typename KmerType=u64, typename IT1=u64, typename IT2=u32>
struct KmerIdx {
    static_assert(std::is_integral<KmerType>::value && std::is_integral<IT1>::value && std::is_integral<IT2>::value, "IT1 and IT2 and KmerType must both be integral");
    static_assert(std::is_unsigned<KmerType>::value && std::is_unsigned<IT1>::value && std::is_unsigned<IT2>::value, "IT1 and IT2 and KmerType must both be unsigned");
    static constexpr KmerType BAD_KMER = std::numeric_limits<KmerType>::max();

    const unsigned k_;
    std::vector<std::string> refnames_;
    std::vector<std::string> comments_;
    ska::flat_hash_map<KmerType, lazy::vector<IT1>> map_;
    
    KmerIdx(unsigned k, const char *path=nullptr): k_(k) {
        if(k > sizeof(KmerType) * CHAR_BIT / 2) RUNTIME_ERROR("Error: k must be <= bits per KmerType / 2.");
        if(path) make_idx(path);
    }
    void add_seq(const kseq_t *ks) {
        refnames_.emplace_back(ks->name.s);
        if(ks->comment.s)
            comments_.emplace_back(ks->comment.s);
        else comments_.emplace_back();
        KmerType kmer = 0, mask = ((1ull << k_) - 1);
        const char *s = ks->seq.s, *e = ks->seq.s + ks->seq.l;
        unsigned nfilled = 0;
        while(s < e) {
            kmer <<= 2;
            if((kmer |= cstr_lut[*s]) == BAD_KMER && *s != 'T') {
                s += k_;
                kmer = 0;
                continue;
            }
            if(++nfilled == k_) {
                kmer &= mask;
                auto diff = s - ks->seq.s;
                auto pos = (refnames_.size() << (sizeof(IT1) / 2 * CHAR_BIT)) | diff;
                auto it = map_.find(kmer);
                if(it == map_.end()) {
                    map_.emplace(kmer, lazy::vector<IT1>{pos});
                } else it->second.emplace_back(pos);
                --nfilled;
            }
        }
    }
    void make_idx(const char *path) {
        gzFile fp = gzopen(path, "rb");
        kseq_t *ks = kseq_init(fp);
        while(kseq_read(ks) >= 0)
            add_seq(ks);
        gzclose(fp);
        kseq_destroy(ks);
    }
    void write(const char *path) const {
        gzFile fp = gzopen(path, "w");
        gzwrite(fp, &k_, sizeof(k_));
        uint32_t nnames = refnames_.size();
        gzwrite(fp, &nnames, sizeof(nnames));
        for(size_t i = 0; i < refnames_.size(); ++i)
            gzputs(fp, refnames_[i].data()), gzputc(fp, '\n');
        for(size_t i = 0; i < comments_.size(); ++i)
            gzputs(fp, comments_[i].data()), gzputc(fp, '\n');
        for(const auto &pair: map_) {
            gzwrite(fp, &pair.first, sizeof(pair.first));
            u32 nelem = pair.second.size();
            gzwrite(fp, &nelem, sizeof(nelem));
            for(const auto v: pair.second)
                gzwrite(fp, &v, sizeof(v));
        }
        gzclose(fp);
    }
};

} // namespace bns

#endif /* #ifndef BONSAI_KMER_IDX_H__ */
