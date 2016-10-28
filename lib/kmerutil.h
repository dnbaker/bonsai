#ifndef _KMER_UTIL_H_
#define _KMER_UTIL_H_
#include <cassert>
#include <cstdint>
#include <cstdio>
#include <unistd.h>
#include <future>
#include <functional>

#include "util.h"

// Converting sequences to numeric equivalent
#ifndef num2nuc
# ifndef NUM2NUC_STR
#  define NUM2NUC_STR "ACGTN"
# endif
# define num2nuc(x) NUM2NUC_STR[(uint8_t)x]
#endif


#ifndef BINFINITY
#    define BINFINITY -1ull
#endif
#ifndef BF
#    define BF BINFINITY
#endif
#ifndef XOR_MASK
#    define XOR_MASK 0xe37e28c4271b5a2dULL
#endif
#define __kmask_init(k) (BF >> (64 - (k << 1)))

#define rc_string  "\0NNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNTNGNNNCNNNNNNNNNNNNANNNNNNNNNNNNtNgNNNcNNNNNNNNNNNNaNNNNNNNNNNN"
#define nuc_cmpl(character) rc_string[(uint8_t)character]

// From bit-twiddling hacks (www.graphics.stanford.edu/~seander/bithacks.html)
#define haszero(v) (((v) - 0x01010101UL) & ~(v) & 0x80808080UL)
#define hasvalue(x,n) (haszero((x) ^ (~0UL/255 * (n))))

namespace kpg {


static const uint32_t nucpos_arr_acgt[128] {0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
        0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
        0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0,
        0, 2, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 3, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
        0, 0, 1, 0, 0, 0, 2, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 3, 0, 0, 0, 0, 0, 0, 0,
        0, 0, 0, 0
};

static INLINE uint8_t nuc2num(char c) {return nucpos_arr_acgt[(uint8_t)c];}
static const int8_t cstr_lut[] {
    -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1,
    -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1,
    -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1,
    -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1,
    -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1,
    -1, -1, -1, -1, -1,  0, -1,  1, -1, -1, -1,  2,
    -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1,
     3, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1,
    -1,  0, -1,  1, -1, -1, -1,  2, -1, -1, -1, -1,
    -1, -1, -1, -1, -1, -1, -1, -1,  3, -1, -1, -1,
    -1, -1, -1, -1, -1, -1, -1, -1
};

static const int8_t cstr_rc_lut[] {
    -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1,
    -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1,
    -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1,
    -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1,
    -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1,
    -1, -1, -1, -1, -1,  3, -1,  2, -1, -1, -1,  1,
    -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1,
     0, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1,
    -1,  3, -1,  2, -1, -1, -1,  1, -1, -1, -1, -1,
    -1, -1, -1, -1, -1, -1, -1, -1,  0, -1, -1, -1,
    -1, -1, -1, -1, -1, -1, -1, -1
};

// C++ std lib doesn't actually give you a way to check on the status directly
// without joining the thread. This is a hacky workaroud c/o
// http://stackoverflow.com/questions/10890242/get-the-status-of-a-stdfuture
template<typename R>
static INLINE bool is_ready(std::future<R> const& f) {
    return f.wait_for(std::chrono::seconds(0)) == std::future_status::ready;
}

// Jellyfish/Kraken
static INLINE uint64_t reverse_complement(uint64_t kmer, uint8_t n) {
    kmer = ((kmer >> 2)  & 0x3333333333333333UL) | ((kmer & 0x3333333333333333UL) << 2);
    kmer = ((kmer >> 4)  & 0x0F0F0F0F0F0F0F0FUL) | ((kmer & 0x0F0F0F0F0F0F0F0FUL) << 4);
    kmer = ((kmer >> 8)  & 0x00FF00FF00FF00FFUL) | ((kmer & 0x00FF00FF00FF00FFUL) << 8);
    kmer = ((kmer >> 16) & 0x0000FFFF0000FFFFUL) | ((kmer & 0x0000FFFF0000FFFFUL) << 16);
    kmer = ( kmer >> 32                        ) | ( kmer                         << 32);
    return (((uint64_t)-1) - kmer) >> (8 * sizeof(kmer) - (n << 1));
}

static INLINE uint64_t canonical_representation(uint64_t kmer, uint8_t n) {
    const uint64_t revcom(reverse_complement(kmer, n));
    return kmer < revcom ? kmer : revcom;
}

} // namespace kpg

#endif //ifndef _KMER_UTIL_H_
