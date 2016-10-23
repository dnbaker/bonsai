#ifndef _SPACE_UTIL_H
#define _SPACE_UTIL_H

#include <vector>
#include "kmerutil.h"

uint32_t comb_size(const std::vector<uint8_t> &spaces) {
    uint32_t ret(spaces.size() + 1); // Since there's 1 fewer entry in spaces
    // We then increment the size of our comb for each space.
    for(auto i: spaces) ret += i;
    return ret;
}

class Spacer {
    // Instance variables
    const std::vector<uint8_t> s_; // Spaces to skip
    const uint32_t k_:8;
    const uint32_t w_:16; // window size
    const uint32_t c_:16; // comb size
    const uint64_t mask_;

public:
    Spacer(unsigned k, uint16_t w, std::vector<uint8_t> *spaces=nullptr):
      k_(k),
      s_(spaces ? *spaces: std::vector<uint8_t>(0, k_ - 1)),
      w_(w),
      c_(comb_size(s_)),
      mask_(__kmask_init(k))
    {}
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
    ~Spacer() {}
};

#endif // #ifndef _SPACE_UTIL_H
