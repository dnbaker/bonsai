#ifndef _SETCMP_H__
#define _SETCMP_H__
#include "bitmap.h"
#include "khash64.h"
#include "feature_min.h"

namespace bns {


template<typename KhashType>
size_t intersection_size(const KhashType *a, const KhashType *b) {
    size_t ret(0);
    for(khiter_t ki(0); ki != kh_end(a); ++ki)
        if(kh_exist(a, ki))
            ret += (khash_get(b, kh_key(a, ki)) != kh_end(b));
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
