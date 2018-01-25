#include "popcnt.h"

namespace popcnt {

unsigned vec_popcnt(const char *p, const size_t l) {
    const uint64_t *arr(reinterpret_cast<const uint64_t *>(p));
    uint64_t ret(0);
    const uint64_t nloops(l >> 3);
    for(size_t i(0); i < nloops; ++i)  ret += popcount(arr[i]);
    for(const char *q(reinterpret_cast<const char *>(arr + nloops)), *e(p + l);
        q != e; ++q) ret += popcount(*q);
    return ret;
}

unsigned vec_popcnt(uint64_t *p, size_t l) {
    if(__builtin_expect(l == 0, 0)) return 0;
    uint64_t ret(popcount(*p));
    while(--l) ret += popcount(*++p);
    return ret;
}

unsigned vec_popcnt(const std::string &vec) {
    return vec_popcnt(vec.data(), vec.size());
}
template<>
auto vec_bitdiff(const bitvec_t &a, const bitvec_t &b);

}
