#ifndef _HASH_H_
#define _HASH_H_
#include "util.h"
#define REGISTER


namespace emp {

// Thomas Wang hash
// Original site down, available at https://naml.us/blog/tag/thomas-wang
// This is our core 64-bit hash.
// It has a 1-1 mapping from any one 64-bit integer to another
// and can be inverted with irving_inv_hash.
constexpr INLINE std::uint64_t wang_hash(std::uint64_t key) {
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
    constexpr std::uint64_t operator()(std::uint64_t key) const {return wang_hash(key);}
};
template<typename T>
struct ptr_wang_hash_struct {
    template<typename = typename std::enable_if<std::is_pointer<T>::value || sizeof(T) == 8>::type>
    constexpr std::uint64_t operator()(T key) const {return wang_hash(reinterpret_cast<std::uint64_t>(key));}
};

template<typename T>
struct idt_struct {
    constexpr T operator()(const T key) const {return key;}
};

// Geoffrey Irving
// https://naml.us/blog/tag/thomas-wang
INLINE std::uint64_t irving_inv_hash(std::uint64_t key) {
  std::uint64_t tmp;

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
    int h = *s;
    if (h) for (++s ; *s; ++s) h = (h << 5) - h + *s;
    return h;
}

// SBDM http://www.cse.yorku.ca/~oz/sdbm.bun/. Mirrored https://github.com/davidar/sdbm.
// Slightly modified.
// This is a 32-bit unsigned hash function for strings.
static INLINE unsigned
dbm_hash(const char *str, std::size_t len)
{
    unsigned n = 0;

#define HASHC  (n = *str++ + (n << 16) + (n << 6) - n)
#ifdef DUFF
    if(len) {
        int loop = (len + 7) >> 3;
        switch(len & 7) {
            case 0: do {
                HASHC;
                case 7: HASHC; case 6: HASHC; case 5: HASHC;
                case 4: HASHC; case 3: HASHC; case 2: HASHC;  case 1: HASHC;
                } while (--loop);
        }
    }
#else
    while (len--) HASHC;
#endif
    return n;
}
#undef HASHC
static INLINE unsigned dbm_hash(const char *str)
{
    unsigned n = *str;
    if(n) for(++str; *str; n = *str++ + ((n << 16) + (n << 6) - n));
    return n;
}

static INLINE unsigned dbm_hash(const std::string str) {return dbm_hash(str.data(), str.size());}

template<typename T>
struct dbm_hash_struct {
    constexpr std::uint64_t operator()(const T &el, size_t nb=sizeof(T)) const {return dbm_hash(reinterpret_cast<const char *>(&el), nb);}
};

} // namespace emp

#endif // #ifdef _HASH_H_
