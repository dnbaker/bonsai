#pragma once
#include "encoder.h"

namespace bns {

template<typename ScoreType>
class GenomeChunker {
    const std::string       path_;
    const size_t        chunk_sz_;
    struct ContigData {
        std::vector<hll::hll_t>  hlls_;
        const std::string contig_name_;
        ContigData(const kseq_t *ks, size_t chunk_size, unsigned p=16, hll::EstimationMethod estim=hll::ERTL_MLE, hll::JointEstimationMethod jestim=hll::ERTL_JOINT_MLE):
                contig_name_(ks->name.s) {
            hlls_.reserve((ks->seq.l + (chunk_size - 1)) / chunk_size);
            while(hlls_.size() < (ks->seq.l + (chunk_size - 1)) / chunk_size) hlls_.emplace_back(p, estim, jestim);
        }
    };
    std::vector<ContigData> cds_;
public:
    GenomeChunker(const char *path, size_t chunk_size, const Spacer &sp, kseq_t *ks=nullptr,
                  unsigned p=16,
                  hll::EstimationMethod estim=hll::ERTL_MLE, hll::JointEstimationMethod jestim=hll::ERTL_JOINT_MLE)
                  : path_(path), chunk_sz_(chunk_size)
    {
        bool destroy = ks == nullptr;
        double csinv = 1./ chunk_sz_;
        gzFile fp = gzopen(path, "rb");
        if(destroy) ks = kseq_init(fp);
        else        kseq_assign(ks, fp);
        Encoder<ScoreType> enc(sp);
        while(kseq_read(ks) >= 0) {
            cds_.emplace_back(ks, chunk_sz_, p, estim, jestim);
            auto &cd = cds_.back();
            enc.for_each([&](const u64 &kmer) {cd.hlls_[size_t(enc.pos() * csinv)].addh(kmer);},
                         ks->seq.s, ks->seq.l);
        }
        gzclose(fp);
        if(destroy) kseq_destroy(ks);
        for(auto &cd: cds_) for(auto &hll: cd.hlls_) hll.csum();
    }
    template<typename Functor>
    void for_each(const Functor &func) const {
        for(const auto &el: cds_) func(el);
    }
    template<typename Functor>
    void for_each(const Functor &func) {
        for(auto &el: cds_) func(el);
    }
    auto chunk_size()  const {return chunk_sz_;}
    const auto &path() const {return path_;}
};

}
