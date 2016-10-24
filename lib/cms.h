#ifndef _CMS_H_
#ifndef _CMS_H_
#include <cstdint>
#include <cinttypes>
#include <cstdio>
#include <cstdlib>
#include <functional>

#include "hash.h"
#include "util.h"

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
    std::vector<T> bits_; // 2-d array of bits.
    uint32_t seeds_[ns];
    const size_t mask_;

    cms_t(size_t sz, uint32_t seedseed=1337) // Lucky number
      sz_(roundup64(sz)), bits_(ns * sz_, 0), mask_(sz_ - 1)
    {
        srand(seedseed);
        for(auto &seed: seeds_) seed = rand();
    }
    INLINE void add(uint64_t hashval, int val) {
        // Multiply a produced hash
        unsigned i(0);
        for(auto seed: seeds)
            bits_[((hashval ^ seed) & mask_) + sz_ * i++] += val;
    }
    INLINE void add(uint64_t hashval) {add(hashval, 1);}
    T query(uint64_t hashval) {
        T ret(std::numeric_limits<T>::max());
        for(auto seed: seeds) {
            const T val((hashval ^ seed) & mask_);
            if(bits_[val + i * ns] < ret) ret = bits_[val + i * ns];
        }
        return ret;
    }
    cms_t<ns> &operator+=(const cms_t<ns> &other, int fail_silently=0) {
        if(other.sz_ != sz_) {
            if(fail_silently) return *this;
            else goto fail;
        }
        for(unsigned i(0); i < ns; ++i) if(seeds_[i] != other.seeds_[i]) goto fail;
        for(unsigned i(0); i < ns * sz_; ++i) bits_[i] += other.bits_[i];
        return *this;
        fail:
        fprintf(stderr, "[%s]: cms_t's have different table sizes. Abort!\n", __func__);
        exit(1);
        return *this;
    }
    ~cms_t() {free(bits);}
};
}

#endif // #ifndef _CMS_H_
