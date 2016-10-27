#ifndef _KPG_UTIL_H_
#define _KPG_UTIL_H_
#include <cstdlib>
#include <cstdio>
#include <cinttypes>
#include <cstdint>
#include <type_traits>
#include "htslib/khash.h"

#ifdef __GNUC__
#  define likely(x) __builtin_expect((x),1)
#  define unlikely(x) __builtin_expect((x),0)
#else
#  define likely(x)
#  define unlikely(x)
#endif

#if __GNUC__ || __clang__
#define INLINE __attribute__((always_inline)) inline
#else
#define INLINE inline
#endif

namespace kpg {

KHASH_SET_INIT_INT64(all)
KHASH_MAP_INIT_INT64(c, uint32_t)
KHASH_MAP_INIT_INT(p, uint32_t)
KHASH_MAP_INIT_STR(name, uint32_t)

size_t count_lines(const char *fn);
khash_t(name) *build_name_hash(const char *fn);
khash_t(name) *destroy_name_hash(khash_t(name) *hash);

template <typename T>
void write_khash_map(T *map, const char *path) {
    FILE *fp(fopen(path, "wb"));
    fwrite("KHASH_DB", 1, 8, fp);
    fwrite(map, 1, sizeof(*map), fp);
    fwrite(map->flags, 1, sizeof(*map->flags) * map->n_buckets, fp);
    fwrite(map->keys, 1, sizeof(*map->keys) * map->n_buckets, fp);
    fwrite(map->vals, 1, sizeof(*map->vals) * map->n_buckets, fp);
    fclose(fp);
}

template <typename T>
T *load_khash_map(const char *path) {
    char magic[9]{0};
    T *rex = (T *)malloc(sizeof(T));
    FILE *fp(fopen(path, "rb"));
    fread(magic, 1, sizeof magic, fp);
    if(memcmp(magic, "KHASH_DB", 8)) {
        fprintf(stderr,  "Missing magic number. (%s found).\n", magic);
        exit(1);
    }
    typedef typename std::remove_pointer<decltype(rex->keys)>::type keytype_t;
    typedef typename std::remove_pointer<decltype(rex->vals)>::type valtype_t;
    fread(rex, 1, sizeof(*rex), fp);
    rex->flags = (uint32_t *)malloc(rex->n_buckets * sizeof(uint32_t));
    fread(rex->flags, sizeof((*rex->flags)) * rex->n_buckets, 1uL, fp);
    rex->keys = (keytype_t *)malloc(rex->n_buckets * sizeof(keytype_t));
    fread(rex->keys, sizeof((*rex->keys)) * rex->n_buckets, 1uL, fp);
    rex->vals = (valtype_t *)malloc(rex->n_buckets * sizeof(valtype_t));
    fread(rex->vals, sizeof(*rex->vals) * rex->n_buckets, 1uL, fp);
    fclose(fp);
    return rex;
}

uint32_t lca(khash_t(p) *map, uint32_t a, uint32_t b);
unsigned node_depth(khash_t(p) *map, uint32_t a);

} // namespace kpg

#endif // #ifdef _KPG_UTIL_H_
