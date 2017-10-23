#ifndef _BIT_CMP_H__
#define _BIT_CMP_H__

#include "lib/tx.h"

namespace emp {

enum BitCmp {
    EQUAL,
    FIRST_PARENT,
    SECOND_PARENT,
    INCOMPARABLE
};


inline int veccmp(const void *a, const void *b, size_t nbytes) {
#if __AVX2__
    __m256i aentry, bentry;
    const __m256i *va = reinterpret_cast<const __m256i*>(a),
                  *vb = reinterpret_cast<const __m256i*>(b);
    size_t n255(nbytes >> 5), nlo(nbytes & 31u);
    int aparent(1), bparent(1); // a could be parent of b, b could be parent of a.
    while(n255--) {
        aentry = _mm256_loadu_si256(va++);
        bentry = _mm256_loadu_si256(vb++);
        bparent &= _mm256_testz_si256(aentry, ~bentry);
        aparent &= _mm256_testz_si256(bentry, ~aentry);
        if((aparent & bparent) == 0) return BitCmp::INCOMPARABLE;
    }
    std::uint8_t *eba((std::uint8_t *)va), *ebb((std::uint8_t *)vb);
#elif __SSE2__
    __m128i aentry, bentry;
    const __m128i *va = reinterpret_cast<const __m128i*>(a),
                  *vb = reinterpret_cast<const __m128i*>(b);
    size_t n127(nbytes >> 4), nlo(nbytes & 15u);
    int aparent(1), bparent(1); // a could be parent of b, b could be parent of a.
    while(n127--) {
        aentry = _mm_loadu_si128(va++);
        bentry = _mm_loadu_si128(vb++);
        bparent &= _mm_testz_si128(aentry, ~bentry);
        aparent &= _mm_testz_si128(bentry, ~aentry);
        if((aparent & bparent) == 0) return BitCmp::INCOMPARABLE;
    }
    std::uint8_t *eba((std::uint8_t *)va), *ebb((std::uint8_t *)vb);
#else
    size_t n64(nbytes >> 3), nlo(nbytes & 7u);
    const std::uint64_t *pa((const std::uint64_t *)a), *pb((const std::uint64_t *)b);
    while(n64--) {
        bparent &= !(*pa & (~*pb));
        aparent &= !(*pb & (~*pa));
        if((aparent & bparent) == 0) return BitCmp::INCOMPARABLE;
        ++pa, ++pb;
    }
    std::uint8_t *eba((std::uint8_t *)pa), std::uint8_t *ebb((std::uint8_t *)pb);
#endif
    while(nlo--) {
        bparent &= !(*eba & (~*ebb));
        aparent &= !(*ebb & (~*eba));
        if((aparent & bparent) == 0) return BitCmp::INCOMPARABLE;
        ++eba, ++ebb;
    }
#if VECCMP_USE_SWITCH
    // Variant (should benchmark)
    switch((aparent << 1) | bparent) {
        case 3: return BitCmp::EQUAL;
        case 2: return BitCmp::FIRST_PARENT;
        case 1: return BitCmp::SECOND_PARENT;
        case 0: return BitCmp::INCOMPARABLE;
    }
#else
    static const std::uint8_t retcodes[] {BitCmp::INCOMPARABLE, BitCmp::SECOND_PARENT, BitCmp::FIRST_PARENT, BitCmp::EQUAL};
    return retcodes[(aparent << 1) | bparent];
#endif
}

// Vector comparison function
template<typename Container1, typename Container2>
int veccmp(const Container1 &a, const Container2 &b) {
    using T = std::decay_t<decltype(a[0])>;
    static_assert(std::is_same<T, std::decay_t<decltype(b[0])>>::value, "a and be must be containers of the same type.");
    if(unlikely(a.size() != b.size())) LOG_EXIT("a and b must be the same size! %zu, %zu\n", static_cast<size_t>(a.size()), static_cast<size_t>(b.size()));
    return veccmp(a.data(), b.data(), a.size() * sizeof(T));
}


}

#endif // #ifndef _BIT_CMP_H__
