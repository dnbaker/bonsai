#include "bits.h"
#include <vector>

namespace emp {

namespace popcnt {


unsigned vec_popcnt(const char *p, const std::size_t l) {
    const std::uint64_t *arr(reinterpret_cast<const std::uint64_t *>(p));
    std::uint64_t ret(0);
    const std::uint64_t nloops(l >> 3);
    for(std::size_t i(0); i < nloops; ++i)  ret += popcount(arr[i]);
    for(const char *q(reinterpret_cast<const char *>(arr + nloops)), *e(p + l);
        q != e; ++q) ret += popcount(*q);
    return ret;
}

unsigned vec_popcnt(std::uint64_t *p, std::size_t l) {
    if(unlikely(l == 0)) return 0;
    std::uint64_t ret(popcount(*p));
    while(--l) ret += popcount(*++p);
    return ret;
}

unsigned vec_popcnt(const std::string &vec) {
    return vec_popcnt(vec.data(), vec.size());
}
template<>
auto vec_bitdiff(const bitvec_t &a, const bitvec_t &b);

} // namespace popcnt

} // namespace emp
