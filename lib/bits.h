#ifndef _BITS_H__
#define _BITS_H__

#include "lib/util.h"

namespace emp {

namespace popcnt {


std::uint64_t vec_popcnt(const std::string &vec);

template<typename T>
INLINE unsigned popcount(T val) noexcept {
    return __builtin_popcount(val);
}

template<>
INLINE unsigned popcount(char val) noexcept {
    return __builtin_popcount((int)val);
}

template<>
INLINE unsigned popcount(unsigned long long val) noexcept {
// From cqf https://github.com/splatlab/cqf/
    asm("popcnt %[val], %[val]"
            : [val] "+r" (val)
            :
            : "cc");
    return val;
}


template<>
INLINE unsigned popcount(unsigned long val) noexcept {
    return __builtin_popcountl(val);
}

template<typename T>
INLINE std::uint64_t vec_popcnt(T &container) {
    auto i(container.cbegin());
    std::uint64_t ret(popcount(*i));
    while(++i != container.cend()) ret += popcount(*i);
    return ret;
}

std::uint64_t vec_popcnt(std::uint64_t *p, std::size_t l);

// Note: This could be further optimized 

} //namespace popcnt


} // namespace emp

#endif
