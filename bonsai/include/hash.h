#ifndef _HASH_H_
#define _HASH_H_
#include "util.h"
#include <climits>

namespace emp {

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

struct wang_hash_struct {

    // But reduces the number of collisions in the hash table because the values will
    // always be divisible by pointer alignment.
    static constexpr size_t OFFSET = log2_64(alignof(void *));

    template<typename U, typename=std::enable_if_t<std::is_pointer_v<U> || std::is_arithmetic_v<U>>>
    constexpr u64 operator()(U key) const {
        static_assert(std::is_pointer_v<U> || std::is_arithmetic_v<U>, "Must be arithmetic  or a pointer.");
        if constexpr(std::is_pointer_v<U>)
            return wang_hash(reinterpret_cast<u64>(key) >> OFFSET);
        else
            return wang_hash(key);
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
#define HASHC  do {(n = *str++ + (n << 16) + (n << 6) - n)} while(0)
    DO_DUFF(len, HASHC);
#undef HASHC
    return n;
}

static INLINE unsigned dbm_hash(const char *str)
{
    unsigned n = *str;
    if(n) for(++str; *str; n = *str++ + ((n << 16) + (n << 6) - n));
    return n;
}

static INLINE unsigned dbm_hash(const std::string str) {return dbm_hash(str.data(), str.size());}

} // namespace emp

#endif // #ifdef _HASH_H_
