#include "bits.h"

namespace emp {

namespace popcnt {


std::uint64_t vec_popcnt(const char *p, const std::size_t l) {
    const std::uint64_t *arr((std::uint64_t *)p);
    std::uint64_t ret(0);
    const std::uint64_t nloops(l >> 3);
    for(std::size_t i(0); i < nloops; ++i)  ret += popcount(arr[i]);
    for(const char *q(reinterpret_cast<const char *>(arr + nloops)), *e(p + l);
        q != e; ++q) ret += popcount(*q);
    return ret;
}

std::uint64_t vec_popcnt(std::uint64_t *p, std::size_t l) {
    std::uint64_t ret(popcount(*p));
    while(l--) ret += popcount(*++p);
    return ret;
}

std::uint64_t vec_popcnt(const std::string &vec) {
    return vec_popcnt(vec.data(), vec.size());
}

} // namespace popcnt

} // namespace emp
