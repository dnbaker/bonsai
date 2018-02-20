#pragma once
#include "tx.h"

namespace emp {

enum BitCmp {
    EQUAL,
    FIRST_PARENT,
    SECOND_PARENT,
    INCOMPARABLE
};

#if !NDEBUG
inline int veccmp_naive(const void *a, const void *b, size_t nbytes);
#endif

inline int veccmp(const void *a, const void *b, size_t nbytes) {
#if __AVX2__
    __m256i aentry, bentry;
    const __m256i *pa = reinterpret_cast<const __m256i*>(a),
                  *pb = reinterpret_cast<const __m256i*>(b);
    size_t n255(nbytes >> 5), nlo(nbytes & 31u);
    int aparent(1), bparent(1); // a could be parent of b, b could be parent of a.
    while(n255--) {
        aentry = _mm256_loadu_si256(pa++);
        bentry = _mm256_loadu_si256(pb++);
        bparent &= _mm256_testz_si256(aentry, ~bentry);
        aparent &= _mm256_testz_si256(bentry, ~aentry);
        if((aparent | bparent) == 0) return BitCmp::INCOMPARABLE;
    }
#elif __SSE2__
    __m128i aentry, bentry;
    const __m128i *pa = reinterpret_cast<const __m128i*>(a),
                  *pb = reinterpret_cast<const __m128i*>(b);
    size_t n127(nbytes >> 4), nlo(nbytes & 15u);
    int aparent(1), bparent(1); // a could be parent of b, b could be parent of a.
    while(n127--) {
        aentry = _mm_loadu_si128(pa++);
        bentry = _mm_loadu_si128(pb++);
        bparent &= _mm_testz_si128(aentry, ~bentry);
        aparent &= _mm_testz_si128(bentry, ~aentry);
        if((aparent == 0) & (bparent == 0)) return BitCmp::INCOMPARABLE;
    }
#else
    size_t n64(nbytes >> 3), nlo(nbytes & 7u);
    const u64 *pa((const u64 *)a), *pb((const u64 *)b);
    while(n64--) {
        bparent &= !(*pa & (~*pb));
        aparent &= !(*pb & (~*pa));
        if((aparent == 0) & (bparent == 0)) return BitCmp::INCOMPARABLE;
        ++pa, ++pb;
    }
#endif
    u8 *eba((u8 *)pa), *ebb((u8 *)pb);
    while(nlo--) {
        bparent &= !(*eba & (~*ebb));
        aparent &= !(*ebb & (~*eba));
        if((aparent == 0) & (bparent == 0)) return BitCmp::INCOMPARABLE;
        ++eba, ++ebb;
    }
    static const u8 retcodes[] {BitCmp::INCOMPARABLE, BitCmp::SECOND_PARENT, BitCmp::FIRST_PARENT, BitCmp::EQUAL};
    assert(retcodes[(aparent << 1) | bparent] == veccmp_naive(a, b, nbytes));
    return retcodes[(aparent << 1) | bparent];
}
#if !NDEBUG
inline int veccmp_naive(const void *a, const void *b, size_t nbytes) {
    const uint64_t *pa((const uint64_t *)a), *pb((const uint64_t *)b);
    int aparent = 1, bparent = 1;
    uint64_t va, vb;
    nbytes >>= 3;
    while(nbytes--) {
        va = *pa++, vb = *pb++;
        if(va & ~vb) {
            bparent = 0;
        }
        if(vb & ~va) {
            aparent = 0;
        }
        if(aparent == bparent && aparent == 0) return BitCmp::INCOMPARABLE;
    }
    static const u8 retcodes[] {BitCmp::INCOMPARABLE, BitCmp::SECOND_PARENT, BitCmp::FIRST_PARENT, BitCmp::EQUAL};
    return retcodes[(aparent << 1) | bparent];
}
#endif

// Vector comparison function
template<typename Container1, typename Container2>
int veccmp(const Container1 &a, const Container2 &b) {
    using T = std::decay_t<decltype(a[0])>;
    static_assert(std::is_same_v<T, std::decay_t<decltype(b[0])>>, "a and be must be containers of the same type.");
    if(unlikely(a.size() != b.size())) LOG_EXIT("a and b must be the same size! %zu, %zu\n", static_cast<size_t>(a.size()), static_cast<size_t>(b.size()));
    return veccmp(a.data(), b.data(), a.size() * sizeof(T));
}


} // namespace emp
