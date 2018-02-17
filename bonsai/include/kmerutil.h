#ifndef _KMER_UTIL_H__
#define _KMER_UTIL_H__
#include <cassert>
#include <cstdio>
#include <unistd.h>
#include <future>
#include "util.h"

// Converting sequences to numeric equivalent
#ifndef num2nuc
# ifndef NUM2NUC_STR
#  define NUM2NUC_STR "ACGTN"
# endif
# define num2nuc(x) NUM2NUC_STR[(uint8_t)x]
#endif


#ifndef BINFINITY
#    define BINFINITY (UINT64_C(-1))
#endif
#ifndef BF
#    define BF BINFINITY
#endif
#define __kmask_init(k) (BF >> (64 - (k << 1)))

#define rc_string  "\0NNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNTNGNNNCNNNNNNNNNNNNANNNNNNNNNNNNtNgNNNcNNNNNNNNNNNNaNNNNNNNNNNN"
#define nuc_cmpl(character) rc_string[(uint8_t)character]

// From bit-twiddling hacks (www.graphics.stanford.edu/~seander/bithacks.html)
#define haszero(v) (((v) - 0x01010101UL) & ~(v) & 0x80808080UL)
#define hasvalue(x,n) (haszero((x) ^ (~0UL/255 * (n))))

namespace emp {


extern const u32 nucpos_arr_acgt[128];

static INLINE uint8_t nuc2num(char c) {return nucpos_arr_acgt[(uint8_t)c];}
extern const int8_t cstr_lut[128];
extern const int8_t cstr_rc_lut[128];
// C++ std lib doesn't actually give you a way to check on the status directly
// without joining the thread. This is a hacky workaroud c/o
// http://stackoverflow.com/questions/10890242/get-the-status-of-a-stdfuture
template<typename R>
static INLINE bool is_ready(std::future<R> const& f) {
    return f.wait_for(std::chrono::seconds(0)) == std::future_status::ready;
}

// Jellyfish/Kraken
static INLINE u64 reverse_complement(u64 kmer, uint8_t n) {
    kmer = ((kmer >> 2)  & 0x3333333333333333UL) | ((kmer & 0x3333333333333333UL) << 2);
    kmer = ((kmer >> 4)  & 0x0F0F0F0F0F0F0F0FUL) | ((kmer & 0x0F0F0F0F0F0F0F0FUL) << 4);
    kmer = ((kmer >> 8)  & 0x00FF00FF00FF00FFUL) | ((kmer & 0x00FF00FF00FF00FFUL) << 8);
    kmer = ((kmer >> 16) & 0x0000FFFF0000FFFFUL) | ((kmer & 0x0000FFFF0000FFFFUL) << 16);
    kmer = ( kmer >> 32                        ) | ( kmer                         << 32);
    return (((u64)-1) - kmer) >> (8 * sizeof(kmer) - (n << 1));
}

// Very slightly modified from bowtie2. [https://github.com/BenLangmead/bowtie2/blob/f2857aae716a3447fa689536c6f9c8b193a05a12/bt2_idx.h#L1847-1895]
inline static void countInU64Ex(uint64_t dw, uint32_t *arrs, unsigned k) {
    const uint64_t COUNT_MASK = (0xFFFFFFFFFFFFFFFF >> (64 - 2 * k));
    static const uint64_t c_table[4] = {
        0xffffffffffffffff,
        0xaaaaaaaaaaaaaaaa,
        0x5555555555555555,
        0x0000000000000000
    };
    uint64_t c0 = c_table[0];
    uint64_t x0 = dw ^ c0;
    uint64_t x1 = (x0 >> 1);
    uint64_t x2 = x1 & UINT64_C(0x5555555555555555);
    uint64_t x3 = x0 & x2;
    x3 &= COUNT_MASK;
    auto tmp = pop::popcount(x3); // because __builtin_popcountll returns a 32-bit element.
    arrs[0] += (uint32_t) tmp;

    c0 = c_table[1];
    x0 = dw ^ c0;
    x1 = (x0 >> 1);
    x2 = x1 & UINT64_C(0x5555555555555555);
    x3 = x0 & x2;
    x3 &= COUNT_MASK;
    tmp = pop::popcount(x3);
    arrs[1] += (uint32_t) tmp;

    c0 = c_table[2];
    x0 = dw ^ c0;
    x1 = (x0 >> 1);
    x2 = x1 & UINT64_C(0x5555555555555555);
    x3 = x0 & x2;
    x3 &= COUNT_MASK;
    tmp = pop::popcount(x3);
    arrs[2] += (uint32_t) tmp;

    c0 = c_table[3];
    x0 = dw ^ c0;
    x1 = (x0 >> 1);
    x2 = x1 & UINT64_C(0x5555555555555555);
    x3 = x0 & x2;
    x3 &= COUNT_MASK;
    arrs[3] += (tmp = pop::popcount(x3));
}

static INLINE u64 canonical_representation(u64 kmer, uint8_t n) {
    const u64 revcom(reverse_complement(kmer, n));
    return kmer < revcom ? kmer : revcom;
}
static INLINE bool canonicalize(u64 &kmer, uint8_t n) {
    const u64 revcom(reverse_complement(kmer, n));
    if(kmer < revcom) return false;
    kmer = revcom;
    return true;
}

} // namespace emp

#endif //ifndef _KMER_UTIL_H__
