#ifndef _EMP_SPACE_H__
#define _EMP_SPACE_H__

#include <vector>
#include <string>
#include <algorithm>
#include "kmerutil.h"
#include "logutil.h"

namespace emp {

typedef std::vector<u8> spvec_t;

u32 comb_size(const spvec_t &spaces);
spvec_t parse_spacing(const char *space_string, unsigned k);
struct Spacer {
    static const u32 max_k = 31;
    // Instance variables
    spvec_t s_; // Spaces to skip
    const u32 k_:8;  // Kmer size
    const u32 c_:16; // comb size
    const u32 w_:16; // window size

public:
    Spacer(unsigned k, std::uint16_t w, spvec_t spaces=spvec_t{}):
      s_(spaces.size() ? spaces: spvec_t(k - 1, 0)),
      k_(k),
      c_(comb_size(s_)),
      w_(std::max((int)c_, (int)w))
    {
        if(k > max_k) LOG_EXIT("Provided k %u greater than max %u.\n", k_, max_k);
        for(auto &i: s_) ++i; // Convert differences into offsets
        if(s_.size() + 1 != k) {
            LOG_EXIT("Error: input vector must have size 1 less than k. k: %u. size: %zu.\n",
                     k, s_.size());
        }
    }
    Spacer(unsigned k, std::uint16_t w, const char *space_string):
        Spacer(k, w, parse_spacing(space_string, k))
    {
    }
    Spacer(unsigned k): Spacer(k, k) {}
    Spacer(const Spacer &other): s_(other.s_), k_(other.k_), c_(other.c_), w_(other.w_) {}
    void write(u64 kmer, std::FILE *fp=stdout) const {
        kmer ^= XOR_MASK;
        int offset = ((k_ - 1) << 1);
        std::fputc(num2nuc((kmer >> offset) & 0x3u), fp);
        for(auto s: s_) {
            assert(offset >= 0);
            offset -= 2;
            while(s-- > 1) std::fputc('-', fp);
            std::fputc(num2nuc((kmer >> offset) & 0x3u), fp);
        }
        std::fputc('\n', fp);
    }
    std::string to_string(u64 kmer) const {
        kmer ^= XOR_MASK;
        std::string ret;
        ret.reserve(c_ - k_ + 1);
        int offset = ((k_ - 1) << 1);
        ret.push_back(num2nuc((kmer >> offset) & 0x3u));
        for(auto s: s_) {
            assert(offset >= 0);
            offset -= 2;
            while(s-- > 1) ret.push_back('-');
            ret.push_back(num2nuc((kmer >> offset) & 0x3u));
        }
        return ret;
    }
    ~Spacer() {}
};

} // namespace emp

#endif // #ifndef _EMP_SPACE_H__
