#ifndef _ENCODER_H_
#define _ENCODER_H_
#include "spacer.h"
#include "htslib/kstring.h"

namespace kpg {

template<int (*is_lt)(uint64_t, uint64_t, void *)>
class Encoder {
    const char *s_;
    const size_t l_;
    const Spacer &sp_;
    unsigned pos_;
    void *data_;
public:
    Encoder(char *s, size_t l, Spacer &sp, void *data=nullptr):
      s_(s),
      l_(l),
      sp_(sp),
      data_(data),
      pos_(0)
    {
    }
    void assign(char *s, size_t l) {
        s_ = s; l_ = l; pos_ = 0;
    }
    void assign(kstring_t *ks) {assign(ks->s, ks->l);}
    uint64_t window(unsigned start) {
        // Encode kmer here.
        assert(sp_.w_ + start <= l_);
        uint64_t best_kmer(0);
        for(unsigned wpos(start), end(sp_.w_ + start); wpos != end; ++wpos) {
            uint64_t new_kmer(cstr_lut[s_[wpos]]);
            unsigned j(0);
            for(auto s: sp_.s_) {
                j += s + 1; ++j;
                assert(j < sp_.c_);
                new_kmer <<= 2;
                new_kmer |= cstr_lut[s_[wpos + j]];
            }
            new_kmer = canonical_representation(new_kmer, sp_.k_);
            if(is_lt(new_kmer, best_kmer, data_)) best_kmer = new_kmer;
        }
        return best_kmer;
    }
    int has_next_window() {return pos_ < l_ - sp_.w_ + 1;}
    uint64_t next_kmer() {
        return window(pos_++);
    }
};

} //namespace kpg
#endif // _ENCODER_H_
