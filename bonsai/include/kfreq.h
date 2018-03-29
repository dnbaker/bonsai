#pragma once
#include "util.h"
#include "kmerutil.h"

namespace emp {

namespace freq {


// Counts short kmer occurrences using arrays. (Supported: up to 16)
template<typename SizeType, typename=std::enable_if_t<std::is_integral_v<SizeType>>>
class KFreqArray {
    struct SubKFreq {
        const unsigned k_; // Kmer size
        u32            v_; // Value
        uint8_t        f_; // How full?
        std::vector<SizeType> data_;
        SubKFreq(unsigned k): k_(k), v_(0), f_(0), data_(1ull << (k << 1)) {}
        void clear_kmer() {
            f_ = v_ = 0;
        }
        void clear() {
            clear_kmer();
            std::fill(std::begin(data_), std::end(data_), 0);
        }
        void write(gzFile fp) const {
            gzwrite(fp, (void *)data_.data(), data_.size() * sizeof(SizeType));
#if !NDEBUG
            std::fprintf(stderr, "For k = %u:", k_);
            for(const auto &el: data_) {
                std::fprintf(stderr, "%u|", (unsigned)el);
            }
            std::fputc('\n', stderr);
#endif
        }
    };
    const unsigned maxk_;
    std::vector<SubKFreq> freqs_;
public:
    KFreqArray(unsigned k): maxk_(k) {
        while(freqs_.size() < maxk_) {
            freqs_.emplace_back(freqs_.size() + 1);
        }
        //if(maxk_ != 4) throw std::runtime_error("I'm making it for only k == 4 for now because I'm lazy.");
    }
    void clear_kmers() {
        for(auto &freq: freqs_) freq.clear_kmer();
    }
#define __kmask32(k) (UINT32_C(-1) >> (32 - (k << 1)))
    void process_seq(const char *s, size_t l) {
        clear_kmers();
        u32 cc;
        size_t i = 0;
        while(i < l) {
            if((cc = cstr_lut[s[i++]]) == UINT32_C(-1)) {
                clear_kmers();
                //std::fprintf(stderr, "i is now %zu. Kmers are cleared. Continue.\n", i);
                continue;
            }
            //std::fprintf(stderr, "Now incrementing counts\n");
            ++freqs_[0].data_[cc];
            for(auto it(freqs_.begin() + 1); it < freqs_.end(); ++it) {
                //std::fprintf(stderr, "Now doing stuff at pos %zu of %zu\n", std::distance(freqs_.begin(), it), freqs_.size());
                auto &sf(*it);
                sf.v_ <<= 2;
                sf.v_ |= cc;
                if(sf.f_ == sf.k_ - 1) {
                    sf.v_ &= __kmask32(sf.k_);
                    ++sf.data_[sf.v_];
                } else ++sf.f_;
            }
            //std::fprintf(stderr, "Done incrementing counts\n");
        }
    }
#undef __kmask32
    void add(const char *path, kseq_t *ks=nullptr) {
        const bool destroy = (ks == nullptr);
        gzFile fp(gzopen(path, "rb"));
        if(destroy) ks = kseq_init(fp);
        else       kseq_assign(ks, fp);

        while(kseq_read(ks) >= 0) process_seq(ks->seq.s, ks->seq.l);

        if(destroy) kseq_destroy(ks);
        gzclose(fp);
    }
    void clear() {
        for(auto &freq: freqs_)
            freq.clear();
    }
    void write(const char *path) {
        gzFile fp = gzopen(path, "wb");
        if(fp == nullptr) throw std::runtime_error("Could not open file for output.");
        gzwrite(fp, (void *)&maxk_, sizeof(maxk_));
        for(const auto &freq: freqs_) {
            freq.write(fp);
        }
        gzclose(fp);
    }
};

using KFC = KFreqArray<u32>;

}

}
