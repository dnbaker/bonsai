#pragma once
#include <cstdint>
#include <cstdlib>
#include <vector>
#include <type_traits>
#include <string>
#include "libpopcnt.h"
#include "hll/hll.h"

namespace pop {
using bitvec_t = std::vector<uint64_t>;

unsigned vec_popcnt(const std::string &vec);
template<typename T>
inline unsigned popcount(T val)    noexcept {return __builtin_popcount(val);}
template<>
inline unsigned popcount(char val) noexcept {
    static const uint8_t lut [] {
        0, 1, 1, 2, 1, 2, 2, 3, 1, 2, 2, 3, 2, 3, 3, 4, 1, 2, 2, 3, 2, 3, 3, 4, 2, 3, 3, 4, 3, 4, 4, 5, 1, 2, 2, 3, 2, 3, 3, 4, 2, 3, 3, 4, 3, 4, 4, 5, 2, 3, 3, 4, 3, 4, 4, 5, 3, 4, 4, 5, 4, 5, 5, 6, 1, 2, 2, 3, 2, 3, 3, 4, 2, 3, 3, 4, 3, 4, 4, 5, 2, 3, 3, 4, 3, 4, 4, 5, 3, 4, 4, 5, 4, 5, 5, 6, 2, 3, 3, 4, 3, 4, 4, 5, 3, 4, 4, 5, 4, 5, 5, 6, 3, 4, 4, 5, 4, 5, 5, 6, 4, 5, 5, 6, 5, 6, 6, 7, 1, 2, 2, 3, 2, 3, 3, 4, 2, 3, 3, 4, 3, 4, 4, 5, 2, 3, 3, 4, 3, 4, 4, 5, 3, 4, 4, 5, 4, 5, 5, 6, 2, 3, 3, 4, 3, 4, 4, 5, 3, 4, 4, 5, 4, 5, 5, 6, 3, 4, 4, 5, 4, 5, 5, 6, 4, 5, 5, 6, 5, 6, 6, 7, 2, 3, 3, 4, 3, 4, 4, 5, 3, 4, 4, 5, 4, 5, 5, 6, 3, 4, 4, 5, 4, 5, 5, 6, 4, 5, 5, 6, 5, 6, 6, 7, 3, 4, 4, 5, 4, 5, 5, 6, 4, 5, 5, 6, 5, 6, 6, 7, 4, 5, 5, 6, 5, 6, 6, 7, 5, 6, 6, 7, 6, 7, 7, 8
    };
    return static_cast<unsigned>(lut[(uint8_t)val]);
}

template<>
inline unsigned popcount(unsigned long long val) noexcept {
#ifndef NO_USE_CQF_ASM
// From cqf https://github.com/splatlab/cqf/
    asm("popcnt %[val], %[val]"
            : [val] "+r" (val)
            :
            : "cc");
    return val;
#else
    // According to GodBolt, gcc7.3 fails to inline this function call even at -Ofast.
    // 
    // 
    return __builtin_popcountll(val);
#endif
}

template<>
inline unsigned popcount(unsigned long val) noexcept {
    return __builtin_popcountl(val);
}

template<typename T>
inline auto vec_popcnt(const T &container) {
    auto i(container.cbegin());
    auto ret(popcount(*i));
    while(++i != container.cend()) ret += popcount(*i);
    return ret;
}

unsigned vec_popcnt(uint64_t *p, size_t l);
    
template<typename T, typename std::enable_if_t<std::is_arithmetic_v<T>>>
inline unsigned bitdiff(T a, T b) {
    // TODO: Modify to use SSE intrinsics to speed up calculation.
    // See https://github.com/WojciechMula/sse-popcount for examples/code.
    // Consider adding #ifndef wrappings based on architecture.
    return popcount(a ^ b);
}

namespace detail {
union SIMDHolder {

public:

#if HAS_AVX_512
    using SType = __m512i;
#define popcnt_fn(x) popcnt512(x)
#elif __AVX2__
#define popcnt_fn(x) popcnt256(x)
    using SType = __m256i;
#else
#define USE_UNROLLED_BITDIFF
#endif
    void add_sum(unsigned &sum) const {
        unroller<0, n64> ur;
        ur(*this, sum);
    }
    template<size_t iternum, size_t niter_left> struct unroller {
        void operator()(const SIMDHolder &ref, unsigned &sum) const {
            sum += ref.val[iternum];
            unroller<iternum+1, niter_left-1>()(ref, sum);
        }
    };
    template<size_t iternum> struct unroller<iternum, 0> {
        void operator()(const SIMDHolder &ref, unsigned &sum) const {}
    };

    static constexpr size_t nels = sizeof(SType) / sizeof(uint8_t);
    static constexpr size_t n64  = sizeof(SType) / sizeof(uint64_t);
    using u8arr = uint8_t[nels];
    SType val;
    u8arr vals;
};

inline unsigned unrolled_bitdiff(const uint64_t *a, const uint64_t *b, size_t nbytes);
inline unsigned byte_bitdiff(const uint8_t *a, const uint8_t *b, size_t nelem) {
    size_t nblocks(nelem / (sizeof(SIMDHolder) / sizeof(uint8_t)));
    const SIMDHolder *pa((const SIMDHolder *)a), *pb((const SIMDHolder *)b);
    SIMDHolder tmp, sum;
    tmp.val = (pa++)->val ^ (pb++)->val; // I'm being lazy here and assuming it's aligned, but I have that guarantee from the aligned vectors.
    sum.val = popcnt_fn(tmp.val);
    while(--nblocks) { // Prefix decrement to account for the fact that I used one block in initialization.
        tmp.val = (pa++)->val ^ (pb++)->val, sum.val += popcnt_fn(tmp.val);
    }
    unsigned ret = unrolled_bitdiff((const uint64_t *)pa, (const uint64_t *)pb, nelem & ((nelem / (sizeof(SIMDHolder) / sizeof(uint8_t))) - 1));
    sum.add_sum(ret);
    return ret;
}

inline unsigned unrolled_bitdiff(const uint64_t *a, const uint64_t *b, size_t nbytes) {
#define ITER ret += popcount(*a++ ^ *b++)
    unsigned ret(popcount(*a++ ^ *b++));
    nbytes -= 8;
    const size_t len(nbytes >> 3); // Num 64-bit integers.
    size_t loop((len + 7) >> 3);
    switch(len & 7) {
        case 0: do {
            ITER;
            case 7: ITER; case 6: ITER; case 5: ITER;
            case 4: ITER; case 3: ITER; case 2: ITER;  case 1: ITER;
        } while (--loop);
    }
    if(__builtin_expect(nbytes &= 0x7u, 0)) {
        uint64_t mask = (0xFFFFFFFFFFFFFFFF >> ((8u - nbytes) << 3));
        const uint64_t ta = *a & mask, tb = *b & mask;
        mask = ta ^ tb; // Re-use mask variable.
        ret += mask = popcount(mask);
    }
    return ret;
#undef ITER
}
} // namespace detail
    
template<typename T>
inline auto vec_bitdiff(const T &a, const T &b) {
#if USE_UNROLLED_BITDIFF
    return detail::unrolled_bitdiff((const uint64_t *)&a[0], (const uint64_t *)&b[0], a.size() * sizeof(a[0]));
#else
    return detail::byte_bitdiff(&a[0], &b[0], a.size() * sizeof(a[0]));
#endif
}

} // namespace popcnt
