#ifndef _KPG_UTIL_H_
#define _KPG_UTIL_H_
#include <cstdlib>
#include <cstdio>
#include <cinttypes>
#include <cstdint>

#if __GNUC__ || __clang__
#define INLINE __attribute__((always_inline))
#else
#define INLINE inline
#endif

#endif // #ifdef _KPG_UTIL_H_
