#ifndef _SETCMP_H__
#define _SETCMP_H__
#include "bitmap.h"
#include "khash64.h"
#include "feature_min.h"

namespace bns {


template<typename KhashType>
size_t intersection_size(const KhashType *a, const KhashType *b) {
    size_t ret(0);
    const KhashType *lhs, *rhs;
    if(a->n_buckets <= b->n_buckets) {
        lhs = a, rhs = b;
    } else lhs = b, rhs = a;
    for(khiter_t ki(0); ki != kh_end(lhs); ++ki)
        if(kh_exist(lhs, ki))
            ret += (khash_get(rhs, kh_key(lhs, ki)) != kh_end(rhs));
    return ret;
}

template<typename KhashType>
double jaccard_index(const KhashType *a, const KhashType *b) {
    const auto is(intersection_size(a, b));
    return static_cast<double>(is) / (kh_size(a) + kh_size(b) - is);
}

template<typename KhashType>
size_t union_size(const KhashType *a, const KhashType *b) {
    size_t ret(0);
    for(khiter_t ki(0); ki != kh_end(a); ++ki)
        if(kh_exist(a, ki))
            ret += (khash_get(b, kh_key(a, ki)) == kh_end(b));
    ret += kh_size(b);
    return ret;
}

} // namespace bns


#endif // #ifndef _SETCMP_H__
