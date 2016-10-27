#ifndef _HASH_H_
#define _HASH_H_
#include "util.h"

namespace kpg {

// Murmurhash3's Finalizer
static INLINE uint64_t u64hash(uint64_t key) {
    key ^= key >> 33;
    key *= 0xff51afd7ed558ccd;
    key ^= key >> 33;
    key *= 0xc4ceb9fe1a85ec53;
    key ^= key >> 33;
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
#define DUFF
static INLINE unsigned
dbm_hash(register const char *str, register size_t len)
{
    register unsigned n = 0;

#define __hashc  (n = *str++ + 65599 * n)
#ifdef DUFF
    if(len) {

        register int loop = (len + 7) >> 3;

        switch(len & (7)) {
        case 0: do {
            __hashc;  case 7: __hashc;
        case 6: __hashc;  case 5: __hashc;
        case 4: __hashc;  case 3: __hashc;
        case 2: __hashc;  case 1: __hashc;
            } while (--loop);
        }
    }
#else
    while (len--) n = *str++ + 65599 * n;
#endif
    return n;
}
#undef DUFF
#undef __hashc
static INLINE unsigned dbm_hash(register const char *str)
{
    register unsigned n = 0;
    while(*str) n = *str++ + ((n << 16) + (n << 6) - n);
    return n;
}

} // namespace kpg

#endif // #ifdef _HASH_H_
