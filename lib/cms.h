#ifndef _EMP_CMS_H__
#define _EMP_CMS_H__
#include <cstdint>
#include <cinttypes>
#include <cstdio>
#include <cstdlib>
#include <cassert>
#include <functional>
#include <limits>
#include <vector>

#include "hash.h"
#include "util.h"

#define is_pow2(x) ((x & (x - 1)) == 0)

namespace emp {


// Hash function taken from Murmurhash3, though we only hash a single single-bit element.
template<typename T, std::size_t ns=4>
struct cms_t {
    const std::size_t sz_;
    std::vector<T> bits_; // 2-d array of bits. Access position j within array i via bits_[j + i *sz_]
    std::uint64_t seeds_[ns];
    const std::size_t mask_;

    cms_t(std::size_t sz, std::uint32_t seedseed=137): // Lucky number
      sz_(roundup64(sz)), bits_(ns * sz_, 0), mask_(sz_ - 1)
    {
        assert(is_pow2(sz_));
        std::srand(seedseed);
        for(auto &seed: seeds_) seed = ((std::uint64_t)std::rand() << 32) | std::rand();
    }
    INLINE void add(std::uint64_t hashval, int val) {
        // Multiply a produced hash
        unsigned i(0);
        for(auto seed: seeds_)
            bits_[((hashval ^ seed) & mask_) + sz_ * i++] += val;
    }
    INLINE void add(std::uint64_t hashval) {add(hashval, 1);}
    T query(std::uint64_t hashval) {
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
        for(std::size_t i(0), end(bits_.size()); i < end; ++i) bits_[i] += other.bits_[i];
        return *this;
        fail:
        fprintf(stderr, "[%s]: cms_t's have different table sizes or seeds. Abort!\n", __func__);
        std::exit(1);
        return *this;
    }
};

} //namspace emp

#endif // #ifndef _EMP_CMS_H__
