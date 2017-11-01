#include "bits.h"
#include <vector>

namespace emp {

namespace popcnt {


unsigned vec_popcnt(const char *p, const std::size_t l) {
    const u64 *arr(reinterpret_cast<const u64 *>(p));
    u64 ret(0);
    const u64 nloops(l >> 3);
    for(std::size_t i(0); i < nloops; ++i)  ret += popcount(arr[i]);
    for(const char *q(reinterpret_cast<const char *>(arr + nloops)), *e(p + l);
        q != e; ++q) ret += popcount(*q);
    return ret;
}

unsigned vec_popcnt(u64 *p, std::size_t l) {
    if(unlikely(l == 0)) return 0;
    u64 ret(popcount(*p));
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
