#ifndef _SPACE_UTIL_H
#define _SPACE_UTIL_H

#include <vector>
#include <string>
#include "kmerutil.h"

namespace kpg {

typedef std::vector<uint8_t> spvec_t;

uint32_t comb_size(const spvec_t &spaces) {
    uint32_t ret(spaces.size() + 1); // Since there's 1 fewer entry in spaces
    // We then increment the size of our comb for each space.
    for(auto i: spaces) ret += i;
    return ret;
}


struct Spacer {
    // Instance variables
    spvec_t s_; // Spaces to skip
    const uint32_t k_:8;
    const uint32_t w_:16; // window size
    const uint32_t c_:16; // comb size
    const uint64_t mask_;

public:
    Spacer(unsigned k, uint16_t w, spvec_t *spaces=nullptr):
      s_(spaces ? *spaces: spvec_t(k - 1, 0)),
      k_(k),
      w_(w),
      c_(comb_size(s_)),
      mask_(__kmask_init(k))
    {
        assert(s_.size() + 1 == k);
        fprintf(stderr, "Size of spaces: %zu. k_: %i\n", s_.size(), (int)k);
    }
    void write(uint64_t kmer, FILE *fp=stdout) {
        int offset = ((k_ - 1) << 1);
        fputc(num2nuc((kmer >> offset) & 0x3u), fp);
        for(auto s: s_) {
            assert(offset >= 0);
            offset -= 2;
            while(s--) fputc('-', fp);
            fputc(num2nuc((kmer >> offset) & 0x3u), fp);
        }
        fputc('\n', fp);
    }
    std::string to_string (uint64_t kmer) const {
        std::string ret;
        ret.reserve(c_ - k_ + 1);
        int offset = ((k_ - 1) << 1);
        fprintf(stderr, "offset: %i. vec size: %zu. c_: %u\n", offset, s_.size(), c_);
        ret.push_back(num2nuc((kmer >> offset) & 0x3u));
        for(auto i(s_.crbegin()); i != s_.crend(); ++i) {
            auto s(*i);
            assert(offset >= 0);
            offset -= 2;
            while(s--) ret.push_back('-');
            ret.push_back(num2nuc((kmer >> offset) & 0x3u));
        }
        return ret;
    }
    ~Spacer() {}
};

} // namespace kpg

#endif // #ifndef _SPACE_UTIL_H
