#pragma once
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
    u8 *eba((u8 *)va), *ebb((u8 *)vb);
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
    u8 *eba((u8 *)va), *ebb((u8 *)vb);
#else
    size_t n64(nbytes >> 3), nlo(nbytes & 7u);
    const u64 *pa((const u64 *)a), *pb((const u64 *)b);
    while(n64--) {
        bparent &= !(*pa & (~*pb));
        aparent &= !(*pb & (~*pa));
        if((aparent & bparent) == 0) return BitCmp::INCOMPARABLE;
        ++pa, ++pb;
    }
    u8 *eba((u8 *)pa), u8 *ebb((u8 *)pb);
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
    static const u8 retcodes[] {BitCmp::INCOMPARABLE, BitCmp::SECOND_PARENT, BitCmp::FIRST_PARENT, BitCmp::EQUAL};
    return retcodes[(aparent << 1) | bparent];
#endif
}

// Vector comparison function
template<typename Container1, typename Container2>
int veccmp(const Container1 &a, const Container2 &b) {
    using T = std::decay_t<decltype(a[0])>;
    static_assert(std::is_same_v<T, std::decay_t<decltype(b[0])>>, "a and be must be containers of the same type.");
    if(unlikely(a.size() != b.size())) LOG_EXIT("a and b must be the same size! %zu, %zu\n", static_cast<size_t>(a.size()), static_cast<size_t>(b.size()));
    return veccmp(a.data(), b.data(), a.size() * sizeof(T));
}


} // namespace emp
