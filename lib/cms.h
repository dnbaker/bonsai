#ifndef _CMS_H_
#define _CMS_H_
#include <cstdint>
#include <cinttypes>
#include <cstdio>
#include <cstdlib>
#include <functional>
#include <limits>
#include <vector>

#include "hash.h"
#include "util.h"

#define is_pow2(x) ((x & (x - 1)) == 0)

namespace kpg {

static inline size_t roundup64(size_t x) {
    --x;
    x |= x >> 1;
    x |= x >> 2;
    x |= x >> 4;
    x |= x >> 8;
    x |= x >> 16;
    x |= x >> 32;
    return ++x;
}

// Hash function taken from Murmurhash3, though we only hash a single single-bit element.
template<typename T, size_t ns=4>
struct cms_t {
    const size_t sz_;
    std::vector<T> bits_; // 2-d array of bits. Access position j within array i via bits_[j + i *sz_]
    uint64_t seeds_[ns];
    const size_t mask_;

    cms_t(size_t sz, uint32_t seedseed=1337): // Lucky number
      sz_(roundup64(sz)), bits_(ns * sz_, 0), mask_(sz_ - 1)
    {
        assert(is_pow2(sz_));
        srand(seedseed);
        for(auto &seed: seeds_) seed = ((uint64_t)rand() << 32) | rand();
    }
    INLINE void add(uint64_t hashval, int val) {
        // Multiply a produced hash
        unsigned i(0);
        for(auto seed: seeds_)
            bits_[((hashval ^ seed) & mask_) + sz_ * i++] += val;
    }
    INLINE void add(uint64_t hashval) {add(hashval, 1);}
    T query(uint64_t hashval) {
        T ret(std::numeric_limits<T>::max());
        unsigned i(0);
        for(auto seed: seeds_) {
            const T val((hashval ^ seed) & mask_);
            if(bits_[val + i * sz_] < ret) ret = bits_[val + i * sz_];
            ++i;
        }
        return ret;
    }
    cms_t<T, ns> &operator+=(const cms_t<T, ns> &other) {
        if(other.sz_ != sz_) goto fail;
        for(unsigned i(0); i < ns; ++i) if(seeds_[i] != other.seeds_[i]) goto fail;
        for(size_t i(0), end(bits_.size()); i < end; ++i) bits_[i] += other.bits_[i];
        return *this;
        fail:
        fprintf(stderr, "[%s]: cms_t's have different table sizes or seeds. Abort!\n", __func__);
        exit(1);
        return *this;
    }
};

} //namspace kpg

#endif // #ifndef _CMS_H_
