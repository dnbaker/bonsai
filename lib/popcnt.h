#pragma once
#include <cstdint>
#include <cstdlib>
#include <vector>
#include <type_traits>
#include <string>

namespace popcnt {
using bitvec_t = std::vector<uint64_t>;

unsigned vec_popcnt(const std::string &vec);
template<typename T>
inline unsigned popcount(T val)    noexcept {return __builtin_popcount(val);}
template<>
inline unsigned popcount(char val) noexcept {
    static const uint8_t lut [] {
        0, 1, 1, 2, 1, 2, 2, 3, 1, 2, 2, 3, 2, 3, 3, 4, 1, 2, 2, 3, 2, 3, 3, 4, 2, 3, 3, 4, 3, 4, 4, 5, 1, 2, 2, 3, 2, 3, 3, 4, 2, 3, 3, 4, 3, 4, 4, 5, 2, 3, 3, 4, 3, 4, 4, 5, 3, 4, 4, 5, 4, 5, 5, 6, 1, 2, 2, 3, 2, 3, 3, 4, 2, 3, 3, 4, 3, 4, 4, 5, 2, 3, 3, 4, 3, 4, 4, 5, 3, 4, 4, 5, 4, 5, 5, 6, 2, 3, 3, 4, 3, 4, 4, 5, 3, 4, 4, 5, 4, 5, 5, 6, 3, 4, 4, 5, 4, 5, 5, 6, 4, 5, 5, 6, 5, 6, 6, 7, 1, 2, 2, 3, 2, 3, 3, 4, 2, 3, 3, 4, 3, 4, 4, 5, 2, 3, 3, 4, 3, 4, 4, 5, 3, 4, 4, 5, 4, 5, 5, 6, 2, 3, 3, 4, 3, 4, 4, 5, 3, 4, 4, 5, 4, 5, 5, 6, 3, 4, 4, 5, 4, 5, 5, 6, 4, 5, 5, 6, 5, 6, 6, 7, 2, 3, 3, 4, 3, 4, 4, 5, 3, 4, 4, 5, 4, 5, 5, 6, 3, 4, 4, 5, 4, 5, 5, 6, 4, 5, 5, 6, 5, 6, 6, 7, 3, 4, 4, 5, 4, 5, 5, 6, 4, 5, 5, 6, 5, 6, 6, 7, 4, 5, 5, 6, 5, 6, 6, 7, 5, 6, 6, 7, 6, 7, 7, 8
    };
    return static_cast<unsigned>(lut[(uint8_t)val]);
}

template<>
inline unsigned popcount(unsigned long long val) noexcept {
#ifndef NO_USE_CQF_ASM
// From cqf https://github.com/splatlab/cqf/
    asm("popcnt %[val], %[val]"
            : [val] "+r" (val)
            :
            : "cc");
    return val;
#else
    return __builtin_popcountll(val);
#endif
}

template<>
inline unsigned popcount(unsigned long val) noexcept {
    return __builtin_popcountl(val);
}

template<typename T>
inline auto vec_popcnt(const T &container) {
    auto i(container.cbegin());
    auto ret(popcount(*i));
    while(++i != container.cend()) ret += popcount(*i);
    return ret;
}

unsigned vec_popcnt(uint64_t *p, size_t l);
    
template<typename T, typename std::enable_if_t<std::is_arithmetic_v<T>>>
inline unsigned bitdiff(T a, T b) {
    // TODO: Modify to use SSE intrinsics to speed up calculation.
    // See https://github.com/WojciechMula/sse-popcount for examples/code.
    // Consider adding #ifndef wrappings based on architecture.
    return popcount(a ^ b);
}
    
template<typename T>
inline auto vec_bitdiff(const T &a, const T &b) {
    assert(a.size() == b.size());
    auto ai(a.cbegin()), bi(b.cbegin());
    auto ret(bitdiff(*ai, *bi));
    while(ai != a.cend()) ret += bitdiff(*ai++, *bi++);
    return ret;
}

} // namespace popcnt
