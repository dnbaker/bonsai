#ifndef _BITS_H__
#define _BITS_H__

#include "lib/util.h"

namespace emp {

namespace popcnt {


std::uint64_t vec_popcnt(const std::string &vec);

template<typename T>
constexpr unsigned popcount(T val) noexcept {
    return __builtin_popcount(val);
}

template<>
constexpr unsigned popcount(char val) noexcept {
    return __builtin_popcount((int)val);
}

template<>
constexpr unsigned popcount(unsigned long long val) noexcept {
    return __builtin_popcountll(val);
}

template<>
constexpr unsigned popcount(unsigned long val) noexcept {
    return __builtin_popcountl(val);
}

template<typename T>
std::uint64_t vec_popcnt(T &container) {
    auto i(container.cbegin());
    std::uint64_t ret(popcount(*i));
    while(++i != container.cend()) ret += popcount(*i);
    return ret;
}

std::uint64_t vec_popcnt(std::uint64_t *p, std::size_t l);

} //namespace popcnt


} // namespace emp

#endif
