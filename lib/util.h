#ifndef _KPG_UTIL_H_
#define _KPG_UTIL_H_
#include <cstdlib>
#include <cstdio>
#include <cinttypes>
#include <cstdint>
#include <type_traits>
#include "htslib/khash.h"
#include "lib/logutil.h"

#ifdef __GNUC__
#  define likely(x) __builtin_expect((x),1)
#  define unlikely(x) __builtin_expect((x),0)
#  define UNUSED(x) __attribute__((unused)) x
#else
#  define likely(x) (x)
#  define unlikely(x) (x)
#  define UNUSED(x) (x)
#endif

#if __GNUC__ || __clang__
#define INLINE __attribute__((always_inline)) inline
#else
#define INLINE inline
#endif

namespace kpg {

KHASH_SET_INIT_INT64(all)
KHASH_MAP_INIT_INT64(c, uint32_t)
KHASH_MAP_INIT_INT64(64, uint64_t)
KHASH_MAP_INIT_INT(p, uint32_t)
KHASH_MAP_INIT_STR(name, uint32_t)

size_t count_lines(const char *fn);
khash_t(name) *build_name_hash(const char *fn);
void destroy_name_hash(khash_t(name) *hash);
khash_t(p) *build_parent_map(const char *fn);

static INLINE size_t roundup64(size_t x) {
    --x;
    x |= x >> 1;
    x |= x >> 2;
    x |= x >> 4;
    x |= x >> 8;
    x |= x >> 16;
    x |= x >> 32;
    return ++x;
}

template <typename T>
void print_khash(T *rex) {
    fprintf(stderr, "n buckets %u, nocc %u, size %u, upper bound %u.\n",
            rex->n_buckets, rex->n_occupied, rex->size, rex->upper_bound);
}

#define __fw(item, fp) \
  fwrite(&(item), 1, sizeof(item), fp)
template <typename T>
void write_khash_map(T *map, const char *path) {
    FILE *fp(fopen(path, "wb"));
    __fw(map->n_buckets, fp);
    __fw(map->n_occupied, fp);
    __fw(map->size, fp);
    __fw(map->upper_bound, fp);
    fwrite(map->flags, __ac_fsize(map->n_buckets), sizeof(*map->flags), fp);
    fwrite(map->keys, map->n_buckets, sizeof(*map->keys), fp);
    fwrite(map->vals, map->n_buckets, sizeof(*map->vals), fp);
    fclose(fp);
}
#undef __fw

template <typename T>
T *load_khash_map(const char *path) {
    T *rex = (T *)calloc(1, sizeof(T));
    FILE *fp(fopen(path, "rb"));
    typedef typename std::remove_pointer<decltype(rex->keys)>::type keytype_t;
    typedef typename std::remove_pointer<decltype(rex->vals)>::type valtype_t;
    fread(&rex->n_buckets, 1, sizeof(rex->n_buckets), fp);
    fread(&rex->n_occupied, 1, sizeof(rex->n_occupied), fp);
    fread(&rex->size, 1, sizeof(rex->size), fp);
    fread(&rex->upper_bound, 1, sizeof(rex->upper_bound), fp);
    rex->flags = (uint32_t *)malloc(sizeof(*rex->flags) * __ac_fsize(rex->n_buckets));
    rex->keys = (keytype_t *)malloc(sizeof(*rex->keys) * rex->n_buckets);
    rex->vals = (valtype_t *)malloc(sizeof(*rex->vals) * rex->n_buckets);
    fread(rex->flags, __ac_fsize(rex->n_buckets), sizeof(*rex->flags), fp);
    fread(rex->keys, 1, rex->n_buckets * sizeof(*rex->keys), fp);
    fread(rex->vals, 1, rex->n_buckets * sizeof(*rex->vals), fp);
    fclose(fp);
    return rex;
}

inline void kset_union(khash_t(all) *a, khash_t(all) *b) {
    khint_t ki2;
    int khr;
    for(ki2 = kh_begin(b); ki2 != kh_end(b); ++ki2)
        if(kh_exist(b, ki2))
            if(kh_get(all, a, kh_key(b, ki2)) == kh_end(a))
                kh_put(all, a, kh_key(b, ki2), &khr);
}

uint32_t lca(khash_t(p) *map, uint32_t a, uint32_t b);
unsigned node_depth(khash_t(p) *map, uint32_t a);

} // namespace kpg

#endif // #ifdef _KPG_UTIL_H_
