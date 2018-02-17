#include "popcnt.h"

namespace pop {

unsigned vec_popcnt(const char *p, const size_t l) {
    return ::popcnt((void *)p, l);
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
