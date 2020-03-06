#ifndef BONSAI_KMER_IDX_H__
#define BONSAI_KMER_IDX_H__
#include "encoder.h"
#include "lazy/vector.h"
#include "flat_hash_map/flat_hash_map.hpp"

namespace bns {

#define KI_SWAP_DESTROY(x) do {decltype(x) tmp = std::move(x);} while(0)

template<typename KmerType=u64, typename IT1=u64>
struct KmerIdx {
    static_assert(std::is_integral<KmerType>::value && std::is_integral<IT1>::value, "IT1 and KmerType must both be integral");
    static_assert(std::is_unsigned<KmerType>::value && std::is_unsigned<IT1>::value, "IT1 and KmerType must both be unsigned");
    static constexpr KmerType BAD_KMER = std::numeric_limits<KmerType>::max();

    const unsigned k_;
    std::vector<std::string> refnames_;
    std::vector<std::string> comments_;
    std::vector<uint64_t>     seqlens_;
    std::vector<uint64_t>     cm_seqs_; // Cumulative sequence lengths
    ska::flat_hash_map<KmerType, lazy::vector<IT1>> map_;

    KmerIdx(unsigned k, const char *path=nullptr): k_(k) {
        if(k > sizeof(KmerType) * CHAR_BIT / 2) UNRECOVERABLE_ERROR("Error: k must be <= bits per KmerType / 2.");
        cm_seqs_.emplace_back(0);
        if(path) make_idx(path);
    }
    void add_seq(const kseq_t *ks) {
        refnames_.emplace_back(ks->name.s);
        if(ks->comment.s)
            comments_.emplace_back(ks->comment.s);
        else comments_.emplace_back();
        seqlens_.emplace_back(ks->seq.l);
        IT1 position_increment = cm_seqs_.back();
        cm_seqs_.emplace_back(cm_seqs_.back() + ks->seq.l);
        KmerType kmer = 0;
        const KmerType mask = ((1ull << k_) - 1);
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
                IT1 diff = s - ks->seq.s;
                auto it = map_.find(kmer);
                if(it == map_.end()) {
                    map_.emplace(kmer, lazy::vector<IT1>{diff + position_increment});
                } else it->second.emplace_back(diff + position_increment);
                --nfilled;
            }
            ++s;
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
    void read(const char *path) {
        KI_SWAP_DESTROY(refnames_);
        KI_SWAP_DESTROY(comments_);
        KI_SWAP_DESTROY(seqlens_);
        KI_SWAP_DESTROY(cm_seqs_);
        KI_SWAP_DESTROY(map_);
        gzFile fp = gzopen(path, "rb");
        if(!fp) UNRECOVERABLE_ERROR("ZOMG");
        gzread(fp, &k_, sizeof(k_));
        uint32_t nnames;
        gzread(fp, &nnames, sizeof(nnames));
        for(auto i = 0u; i < nnames; ++i) {
            uint64_t len;
            gzread(fp, &len, sizeof(len));
            seqlens_.emplace_back(len);
        }
        std::string tmp;
        for(auto i = 0u; i < nnames; ++i) {
            char buf[1 << 12];
            if(!gzgets(fp, buf, sizeof(buf))) UNRECOVERABLE_ERROR("ZOMG");
            tmp = buf;
            tmp.pop_back();
            refnames_.push_back(std::move(tmp));
        }
        for(auto i = 0u; i < nnames; ++i) {
            char buf[1 << 12];
            if(!gzgets(fp, buf, sizeof(buf))) UNRECOVERABLE_ERROR("ZOMG");
            tmp = buf;
            tmp.pop_back();
            comments_.push_back(std::move(tmp));
        }
        FOREVER {
            KmerType kmer;
            uint32_t nelem;
            if(unlikely(gzread(fp,  &kmer, sizeof(kmer )) != sizeof(kmer)))  UNRECOVERABLE_ERROR("ZOMG");
            if(unlikely(gzread(fp, &nelem, sizeof(nelem)) != sizeof(nelem))) UNRECOVERABLE_ERROR("ZOMG");
            auto it = map_.emplace(nelem, lazy::vector<IT1>()).first;
            IT1 tmp;
            it->second.reserve(nelem);
            for(size_t i = 0; i < nelem; ++i)
                gzread(fp, &tmp, sizeof(tmp)), it->second.emplace_back(tmp);
            if(gzeof(fp)) break;
        }
    }
    void write(const char *path) const {
        gzFile fp = gzopen(path, "w");
        gzwrite(fp, &k_, sizeof(k_));
        uint32_t nnames = refnames_.size();
        gzwrite(fp, &nnames, sizeof(nnames));
        for(size_t i = 0; i < nnames; ++i)
            gzwrite(fp, &seqlens_[i], sizeof(seqlens_[i]));
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
