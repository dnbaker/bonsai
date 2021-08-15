#ifndef _HASH_H_
#define _HASH_H_
#include <climits>
#include <cassert>
#include <type_traits>
#include "kseq_declare.h"

#ifndef DO_DUFF
#define DO_DUFF(len, ITER) \
    do { \
        if(len) {\
            std::uint64_t loop = (len + 7) >> 3;\
            switch(len & 7) {\
                case 0: do {\
                    ITER; [[fallthrough]];\
                    case 7: ITER; [[fallthrough]]; case 6: ITER; [[fallthrough]]; case 5: ITER; [[fallthrough]];\
                    case 4: ITER; [[fallthrough]]; case 3: ITER; [[fallthrough]]; case 2: ITER; [[fallthrough]]; case 1: ITER;\
                } while (--loop);\
            }\
        }\
    } while(0)
#endif

namespace bns {
using u64 = uint64_t;

// Thomas Wang hash
// Original site down, available at https://naml.us/blog/tag/thomas-wang
// This is our core 64-bit hash.
// It has a 1-1 mapping from any one 64-bit integer to another
// and can be inverted with irving_inv_hash.
constexpr INLINE u64 wang_hash(u64 key) {
  key = (~key) + (key << 21); // key = (key << 21) - key - 1;
  key = key ^ (key >> 24);
  key = (key + (key << 3)) + (key << 8); // key * 265
  key = key ^ (key >> 14);
  key = (key + (key << 2)) + (key << 4); // key * 21
  key = key ^ (key >> 28);
  key = key + (key << 31);
  return key;
}

static inline constexpr int log2_64(uint64_t value)
{
    // https://stackoverflow.com/questions/11376288/fast-computing-of-log2-for-64-bit-integers
    const int tab64[64] {
        63,  0, 58,  1, 59, 47, 53,  2,
        60, 39, 48, 27, 54, 33, 42,  3,
        61, 51, 37, 40, 49, 18, 28, 20,
        55, 30, 34, 11, 43, 14, 22,  4,
        62, 57, 46, 52, 38, 26, 32, 41,
        50, 36, 17, 19, 29, 10, 13, 21,
        56, 45, 25, 31, 35, 16,  9, 12,
        44, 24, 15,  8, 23,  7,  6,  5
    };
    value |= value >> 1;
    value |= value >> 2;
    value |= value >> 4;
    value |= value >> 8;
    value |= value >> 16;
    value |= value >> 32;
    // This could be replaced with a __builtin_clz
    return tab64[((uint64_t)((value - (value >> 1))*0x07EDD5E59A4E28C2)) >> 58];
}

struct wang_hash_struct {

    // But reduces the number of collisions in the hash table because the values will
    // always be divisible by pointer alignment.
    static constexpr size_t OFFSET = log2_64(alignof(void *));

    template<typename U, typename=typename std::enable_if<std::is_pointer<U>::value || std::is_arithmetic<U>::value>::type>
    constexpr u64 operator()(U key) const {
        static_assert(std::is_pointer<U>::value || std::is_arithmetic<U>::value, "Must be arithmetic  or a pointer.");
#if __cplusplus < 201703L
        if (std::is_pointer<U>::value)
#else
        if constexpr(std::is_pointer<U>::value)
#endif
            return wang_hash((u64)(key) >> OFFSET);
        else
            return wang_hash((u64)(key));
    }

};

template<typename T>
struct idt_struct {
    constexpr T operator()(const T key) const {return key;}
};

// Geoffrey Irving
// https://naml.us/blog/tag/thomas-wang
INLINE u64 irving_inv_hash(u64 key) {
  u64 tmp;
  // Invert key = key + (key << 31)
  tmp = key-(key<<31);
  key = key-(tmp<<31);
  // Invert key = key ^ (key >> 28)
  tmp = key^key>>28;
  key = key^tmp>>28;
  // Invert key *= 21
  key *= 14933078535860113213u;
  // Invert key = key ^ (key >> 14)
  tmp = key^key>>14;
  tmp = key^tmp>>14;
  tmp = key^tmp>>14;
  key = key^tmp>>14;
  // Invert key *= 265
  key *= 15244667743933553977u;
  // Invert key = key ^ (key >> 24)
  tmp = key^key>>24;
  key = key^tmp>>24;
  // Invert key = (~key) + (key << 21)
  tmp = ~key;
  tmp = ~(key-(tmp<<21));
  tmp = ~(key-(tmp<<21));
  key = ~(key-(tmp<<21));
  return key;
}


// Taken from khash. https://github.com/attractivechaos/klib
static INLINE int X31_hash_string(const char *s)
{
    int h = *s++;
    if (h) while(*s) h = (h << 5) - h + *s++;
    return h;
}

// SBDM http://www.cse.yorku.ca/~oz/sdbm.bun/. Mirrored https://github.com/davidar/sdbm.
// Slightly modified.
// This is a 32-bit unsigned hash function for strings.
static INLINE unsigned
dbm_hash(const char *str, size_t len)
{
    unsigned n = *str++;
    if(n) {
#define HASHC  do {(n = *str++ + (n << 16) + (n << 6) - n);} while(0)
        --len;
        DO_DUFF(len, HASHC);
#undef HASHC
    }
    return n;
}

static INLINE unsigned dbm_hash(const char *str)
{
    unsigned n = *str;
    if(n) for(++str; *str; n = *str++ + ((n << 16) + (n << 6) - n));
    return n;
}

static INLINE unsigned dbm_hash(const std::string str) {return dbm_hash(str.data(), str.size());}

// rotate "v" to the left 1 position
template<size_t n>
INLINE constexpr uint64_t lrot(uint64_t x) {
    return (x << n) | (x >> (64 - n));
}
template<size_t n>
INLINE constexpr uint64_t rrot(uint64_t x) {
    return (x >> n) | (x << (64 - n));
}

// rotate 31-left bits of "v" to the left by "s" positions
// Precondition: s must be less than lbits
template<unsigned lbits>
INLINE constexpr uint64_t rol(uint64_t x, unsigned s) {
    s = std::min(s, uint32_t(s - lbits));
    assert(s < lbits);
    return ((x << s) | (x >> (lbits - s))) & ((1ul << lbits) - 1);
}


// From bit twiddling hacks
// Swaps ranges of bits starting at b1 and b2
template<unsigned b1, unsigned b2, unsigned range=1>
INLINE constexpr uint64_t swapbits(uint64_t x) {
    static_assert(b1 <= 64 && b2 <= 64, "Can't swap bits which don't exist.");
    const uint64_t tmp = ((x >> b1) ^ (x >> b2)) & ((1ul << range) - 1);
    x ^= (tmp << b1) | (tmp << b2);
    return x;
}
inline uint64_t swapbits033(const uint64_t v) {return swapbits<0, 33>(v);}
inline uint64_t swapbits3263(const uint64_t v) {return swapbits<32, 63>(v);}

// rotate 33-right bits of "v" to the left by "s" positions
inline uint64_t rol33(const uint64_t v, unsigned s) {
    return rol<33>(v, s);
}

} // namespace bns

#endif // #ifdef _HASH_H_
