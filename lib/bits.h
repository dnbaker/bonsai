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
std::uint64_t vec_popcnt(T &container) {
    auto i(container.cbegin());
    std::uint64_t ret(popcount(*i));
    while(++i != container.cend()) ret += popcount(*i);
    return ret;
}

std::uint64_t vec_popcnt(std::uint64_t *p, std::size_t l);
    
    
INLINE unsigned bitdiff(std::uint64_t a, std::uint64_t b) {
    return popcount(a ^ b);
}
    
template<typename T>
std::uint64_t vec_bitdiff(T &a, T &b) {
    auto ai(a.cbegin()), bi(b.cbegin());
    std::uint64_t ret(bitdiff(*ai, *bi));
    while(++ai != a.cend() && ++bi != b.cend()) ret += bitdiff(*ai, *bi);
    return ret;
}

} //namespace popcnt


} // namespace emp

#endif
