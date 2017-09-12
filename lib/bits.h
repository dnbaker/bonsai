#ifndef _BITS_H__
#define _BITS_H__
#include <type_traits>
#include "lib/util.h"

namespace emp {

namespace popcnt {


unsigned vec_popcnt(const std::string &vec);

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
INLINE auto vec_popcnt(T &container) {
    auto i(container.cbegin());
    auto ret(popcount(*i));
    while(i != container.cend()) ret += popcount(*i++);
    return ret;
}

unsigned vec_popcnt(std::uint64_t *p, std::size_t l);
    
template<typename T, typename std::enable_if_t<std::is_arithmetic<T>::value>>
INLINE unsigned bitdiff(T a, T b) {
    return popcount(a ^ b);
}
    
template<typename T>
INLINE auto vec_bitdiff(const T &a, const T &b) {
    assert(a.size() == b.size());
    auto ai(a.cbegin()), bi(b.cbegin());
    auto ret(bitdiff(*ai, *bi));
    while(ai != a.cend()) ret += bitdiff(*ai++, *bi++);
    return ret;
}

} //namespace popcnt


} // namespace emp

#endif
