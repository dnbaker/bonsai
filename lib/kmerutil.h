#ifndef _KMER_UTIL_H_
#define _KMER_UTIL_H_
#include <cassert>
#include <cstdint>
#include <cstdio>
#include <unistd.h>
#include <future>
#include <functional>

#include "htslib/khash.h"

// Largest odd kmer that can held in 64 bits, since 32 Ts is reserved to signal ambiguity
#define MAX_KMER 31

// Utility
#define MAX2(x, y) (x > y ? x: y)

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
#define __kmask_init(k) (BF >> (64 - (k * 2)))

#define rc_string  "\0NNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNTNGNNNCNNNNNNNNNNNNANNNNNNNNNNNNtNgNNNcNNNNNNNNNNNNaNNNNNNNNNNN"
#define nuc_cmpl(character) rc_string[(uint8_t)character]

// From bit-twiddling hacks (www.graphics.stanford.edu/~seander/bithacks.html)
#define haszero(v) (((v) - 0x01010101UL) & ~(v) & 0x80808080UL)
#define hasvalue(x,n) (haszero((x) ^ (~0UL/255 * (n))))

namespace kpg {

KHASH_MAP_INIT_INT64(kmer, uint64_t)

static const uint32_t nucpos_arr_acgt[128] {0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
        0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
        0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0,
        0, 2, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 3, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
        0, 0, 1, 0, 0, 0, 2, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 3, 0, 0, 0, 0, 0, 0, 0,
        0, 0, 0, 0
};

static inline uint8_t nuc2num(char c) {return nucpos_arr_acgt[(uint8_t)c];}

static const int8_t htseq_lut[] {
    -1,  0,  1, -1,  2, -1, -1, -1,
     3, -1, -1, -1, -1, -1, -1, -1
};
static const int8_t htseq_rc_lut[] {
    -1,  3,  2, -1,  1, -1, -1, -1,
     0, -1, -1, -1, -1, -1, -1, -1
};
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

//Gets the 2 bits for an encoded base at an index in an unsigned 64-bit integer encoded kmer.
#define encoded_base(kmer, k, index) ((kmer >> (2 * (k - 1 - index))) & 0x3)
static inline char decoded_base(uint64_t kmer, const int k, const int index) {
   return num2nuc(encoded_base(kmer, k, index));
}

inline void kmer2cstr(uint64_t kmer, int k, char *buf)
{
    while(k) *buf++ = num2nuc((kmer >> (2 * --k)) & 0x3u);
    *buf = '\0';
}

// Used to determine the direction in which to encode a kmer
inline int cstr_rc_lt(const char *seq, int k, int cpos) {
    const char *_seq1 = cpos + seq, *_seq2 = _seq1 + k - 1;
    for(;k;--k, ++_seq1, --_seq2)
        if(*_seq1 != nuc_cmpl(*_seq2))
            return *_seq1 < nuc_cmpl(*_seq2);
    return 0; // This is reverse-complementarily palindromic. Doesn't matter: it's the same string.
}

//Get unspaced kmer at position.
inline uint64_t cstr_get_kmer(const char *seq, const int cpos,
                              uint64_t kmer_mask, int k) {
    uint64_t ret = 0;
    //LOG_DEBUG("Seq: %s at pointer %p.", seq, (void *)seq);
    if(cstr_rc_lt(seq, k, cpos)) {
        seq += cpos;
        for(;k; --k) {
            ret <<= 2;
            ret |= cstr_lut[static_cast<uint8_t>(*seq++)];
            if(ret == BF) {/*fprintf(stderr, "BINFINITY1. seq: %c, %i.\n", *(seq - 1), (int)*(seq - 1)); */return ret;}
        }
        ret &= kmer_mask;
    } else {
        seq += cpos + k;
        for(;k; --k) {
            ret <<= 2;
            ret |= cstr_rc_lut[static_cast<uint8_t>(*--seq)];
            if(ret == BF) {/*fprintf(stderr, "BINFINITY2: seq: %c.\n", *seq); */ return ret;}
        }
        ret &= kmer_mask;
    }
    return ret;
} /*cstr_set_kmer*/

inline void cstr_set_kmer(uint64_t *ret, const char *seq, const int cpos,
                          uint64_t kmer_mask, int k) {
    *ret = 0;
    //LOG_DEBUG("Seq: %s at pointer %p.", seq, (void *)seq);
    if(cstr_rc_lt(seq, k, cpos)) {
        seq += cpos;
        for(;k; --k) {
            *ret <<= 2;
            *ret |= cstr_lut[static_cast<uint8_t>(*seq++)];
            if(*ret == BF) return;
        }
        *ret &= kmer_mask;
    } else {
        seq += cpos + k;
        for(;k; --k) {
            *ret <<= 2;
            *ret |= cstr_rc_lut[static_cast<uint8_t>(*--seq)];
            if(*ret == BF) return;
        }
        *ret &= kmer_mask;
    }
} /*cstr_set_kmer*/

// C++ std lib doesn't actually give you a way to check on the status directly
// without joining the thread. This is a hacky workaroud c/o
// http://stackoverflow.com/questions/10890242/get-the-status-of-a-stdfuture
template<typename R>
static inline bool is_ready(std::future<R> const& f) {
    return f.wait_for(std::chrono::seconds(0)) == std::future_status::ready;
}

// Jellyfish/Kraken
static inline uint64_t reverse_complement(uint64_t kmer, uint8_t n) {
    kmer = ((kmer >> 2)  & 0x3333333333333333UL) | ((kmer & 0x3333333333333333UL) << 2);
    kmer = ((kmer >> 4)  & 0x0F0F0F0F0F0F0F0FUL) | ((kmer & 0x0F0F0F0F0F0F0F0FUL) << 4);
    kmer = ((kmer >> 8)  & 0x00FF00FF00FF00FFUL) | ((kmer & 0x00FF00FF00FF00FFUL) << 8);
    kmer = ((kmer >> 16) & 0x0000FFFF0000FFFFUL) | ((kmer & 0x0000FFFF0000FFFFUL) << 16);
    kmer = ( kmer >> 32                        ) | ( kmer                         << 32);
    return (((uint64_t)-1) - kmer) >> (8 * sizeof(kmer) - (n << 1));
}

static inline uint64_t canonical_representation(uint64_t kmer, uint8_t n) {
    const uint64_t revcom(reverse_complement(kmer, n));
    return kmer < revcom ? kmer : revcom;
}

// Calculates lexicographic minimizer for unspaced seed.
template<uint8_t k, uint8_t m>
uint64_t kmer_minimizer(uint64_t kmer) {
    static const uint64_t mask((1 << (m * 2)) - 1);
    static const uint8_t key_bits(k * 2);
    static const uint64_t xor_mask(XOR_MASK & mask);
    uint64_t min_bin_key(~0);
    for (unsigned i(0); i < key_bits / 2 - m + 1; ++i) {
        const uint64_t temp_bin_key(xor_mask ^ canonical_representation(kmer & mask, m));
        if (temp_bin_key < min_bin_key)
            min_bin_key = temp_bin_key;
        kmer >>= 2;
    }
    return min_bin_key;
}

// Minimizer is now the floor of k / 2.
template<uint8_t k>
uint64_t kmer_minimizer(uint64_t kmer) {
    static const uint64_t mask((1 << (k / 2 * 2)) - 1);
    static const uint8_t key_bits(k * 2);
    static const uint64_t xor_mask(XOR_MASK & mask);
    uint64_t min_bin_key(~0);
    for (unsigned i(0); i < key_bits / 2 - k / 2 + 1; ++i) {
        const uint64_t temp_bin_key(xor_mask ^ canonical_representation(kmer & mask, k / 2));
        if (temp_bin_key < min_bin_key)
            min_bin_key = temp_bin_key;
        kmer >>= 2;
    }
    return min_bin_key;
}

// Taken from Kraken
static inline uint64_t kmer_minimizer(uint64_t kmer, const uint8_t k, const uint8_t m) {
    const uint64_t mask((1 << (m * 2)) - 1);
    const uint8_t key_bits(k * 2);
    const uint64_t xor_mask(XOR_MASK & mask);
    uint64_t min_bin_key(~0);
    for (uint64_t i = 0; i < ((key_bits / 2) - m + 1uL); ++i) {
        const uint64_t temp_bin_key(xor_mask ^ canonical_representation(kmer & mask, m));
        if (temp_bin_key < min_bin_key)
            min_bin_key = temp_bin_key;
        kmer >>= 2;
    }
    return min_bin_key;
}

static inline int max_hp(uint64_t kmer)
{
    int run(1), max_run(0);
    uint8_t last(kmer&0x3);
    for(kmer >>= 2; kmer; kmer>>=2)
        if((kmer&0x3) == last) ++run;
        else max_run = MAX2(run, max_run), last = kmer&0x3, run = 1;
    return max_run;
}

static inline int fails_hp(uint64_t kmer, int threshold)
{
    char run(1), last(kmer&0x3);
    for(kmer >>= 2; kmer; kmer>>=2) {
        if((kmer&0x3) == last) {
            if(++run == threshold) return 1;
        } else last = kmer&0x3, run = 0;
    }
    return 0;
}
static inline int fails_hp(char *str, int threshold)
{
    char last(*str), run(1);
    for(++str; *str; ++str) {
        if(*str == last) {
            if(++run == threshold) return 1;
        } else last = *str, run = 0;
    }
    return 0;
}

} // namespace kpg

#endif //ifndef _KMER_UTIL_H_
