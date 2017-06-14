#ifndef _EMP_UTIL_H__
#define _EMP_UTIL_H__
#include <cstdlib>
#include <cstdio>
#include <cinttypes>
#include <cstdint>
#include <type_traits>
#include <map>
#include <forward_list>
#include <vector>
#include <unordered_map>
#include <limits>
#include <chrono>
#include "khash64.h"
#include "lib/logutil.h"
#include "lib/sample_gen.h"

#ifdef __GNUC__
#  ifndef likely
#    define likely(x) __builtin_expect((x),1)
#  endif
#  ifndef unlikely
#    define unlikely(x) __builtin_expect((x),0)
#  endif
#  ifndef UNUSED
#    define UNUSED(x) __attribute__((unused)) x
#  endif
#else
#  ifndef likely
#    define likely(x) (x)
#  endif
#  ifndef unlikely
#    define unlikely(x) (x)
#  endif
#  ifndef UNUSED
#    define UNUSED(x) (x)
#  endif
#endif

#ifndef INLINE
#  if __GNUC__ || __clang__
#  define INLINE __attribute__((always_inline)) inline
#  define PACKED __attribute__((packed))
#  else
#  define INLINE inline
#  define PACKED
#  endif
#endif

#ifndef TAXID_TYPEDEF
#define TAXID_TYPEDEF
using tax_t = std::uint32_t;
#endif

#define TIME_CODE(code, name) do             \
{                                            \
    auto i(std::chrono::system_clock::now());\
    { code }                                 \
    auto j(std::chrono::system_clock::now());\
    fprintf(stderr, "Task %s took %lf\ns", name, std::chrono::duration<double>(j - i).count());\
} while(0)

template<typename T> class TD; // For debugging

namespace emp {

KHASH_SET_INIT_INT64(all)
KHASH_MAP_INIT_INT64(c, tax_t)
KHASH_MAP_INIT_INT64(64, std::uint64_t)
KHASH_MAP_INIT_INT(p, tax_t)
KHASH_MAP_INIT_STR(name, tax_t)
using strlist = std::forward_list<std::string>;

std::string rand_string(std::size_t n);
std::size_t count_lines(const char *fn) noexcept;
khash_t(name) *build_name_hash(const char *fn) noexcept;
void destroy_name_hash(khash_t(name) *hash) noexcept;
khash_t(p) *build_parent_map(const char *fn) noexcept;
std::unordered_map<tax_t, std::vector<tax_t>> invert_parent_map(khash_t(p) *) noexcept;
tax_t get_taxid(const char *fn, khash_t(name) *name_hash);
std::string get_firstline(const char *fn);
std::vector<tax_t> get_all_descendents(const std::unordered_map<tax_t, std::vector<tax_t>> &map, tax_t tax);
// Resolve_tree is modified from Kraken 1 source code, which
// is MIT-licensed. https://github.com/derrickwood/kraken
tax_t resolve_tree(std::map<tax_t, tax_t> &hit_counts,
                   khash_t(p) *parent_map) noexcept;

const char *bool2str(bool val);

// Modified from bit-twiddling hacks to work with 64-bit integers.
template<typename T>
static INLINE auto roundup64(T x) noexcept {
    kroundup64(x);
    return x;
}

INLINE std::uint64_t rand64() noexcept {
    return (((std::uint64_t)std::rand()) << 32) | std::rand();
}
template <typename T> 
void khash_destroy(T *map) noexcept {
    LOG_EXIT("NotImplementedError");
}

template<>
void khash_destroy(khash_t(64) *map) noexcept;
template<>
void khash_destroy(khash_t(all) *map) noexcept;
template<>
void khash_destroy(khash_t(c) *map) noexcept;
template<>
void khash_destroy(khash_t(p) *map) noexcept;

template <typename T>
void print_khash(T *rex) noexcept {
    std::fprintf(stderr, "n buckets %zu, nocc %zu, size %zu, upper bound %zu.\n",
                 rex->n_buckets, rex->n_occupied, rex->size, rex->upper_bound);
}

#define __fw(item, fp) \
  std::fwrite(&(item), 1, sizeof(item), fp)

template<typename T>
void khash_write_impl(T *map, std::FILE *fp) noexcept {
    for(khiter_t ki(0); ki != kh_end(map); ++ki)
        if(!kh_exist(map, ki))
            kh_key(map, ki) = 0, kh_val(map, ki) = 0;
    __fw(map->n_buckets, fp);
    __fw(map->n_occupied, fp);
    __fw(map->size, fp);
    __fw(map->upper_bound, fp);
    std::fwrite(map->flags, __ac_fsize(map->n_buckets), sizeof(*map->flags), fp);
    std::fwrite(map->keys, map->n_buckets, sizeof(*map->keys), fp);
    std::fwrite(map->vals, map->n_buckets, sizeof(*map->vals), fp);
}

template <typename T>
std::size_t khash_write(T *map, const char *path) noexcept {
    std::FILE *fp(std::fopen(path, "wb"));
    khash_write_impl(map, fp);
    std::size_t ret(ftell(fp));
    std::fclose(fp);
    return ret;
}

#undef __fw

template <typename T>
T *khash_load_impl(std::FILE *fp) noexcept {
    T *rex((T *)std::calloc(1, sizeof(T)));
    typedef typename std::remove_pointer<decltype(rex->keys)>::type keytype_t;
    typedef typename std::remove_pointer<decltype(rex->vals)>::type valtype_t;
    std::fread(&rex->n_buckets, 1, sizeof(rex->n_buckets), fp);
    std::fread(&rex->n_occupied, 1, sizeof(rex->n_occupied), fp);
    std::fread(&rex->size, 1, sizeof(rex->size), fp);
    std::fread(&rex->upper_bound, 1, sizeof(rex->upper_bound), fp);
    LOG_DEBUG("buckets: %zu. nocc: %zu. size: %zu. ub: %zu\n", (size_t)rex->n_buckets, size_t(rex->n_occupied), size_t(rex->size), size_t(rex->upper_bound));
    rex->flags = (std::uint32_t *)std::malloc(sizeof(*rex->flags) * __ac_fsize(rex->n_buckets));
    rex->keys = (keytype_t *)std::malloc(sizeof(*rex->keys) * rex->n_buckets);
    if(!rex->keys) std::fprintf(stderr, "Could not allocate %zu bytes of memory (%zu GB)\n", sizeof(*rex->keys) * rex->n_buckets, sizeof(*rex->keys) * rex->n_buckets >> 30), exit(1);
    rex->vals = (valtype_t *)std::malloc(sizeof(*rex->vals) * rex->n_buckets);
    if(!rex->vals) std::fprintf(stderr, "Could not allocate %zu bytes of memory (%zu GB)\n", sizeof(*rex->vals) * rex->n_buckets, sizeof(*rex->vals) * rex->n_buckets >> 30), exit(1);
    LOG_DEBUG("About to read into flags at %p\n", (void *)rex->flags);
    std::fread(rex->flags, __ac_fsize(rex->n_buckets), sizeof(*rex->flags), fp);
    LOG_DEBUG("About to read into keys at %p\n", (void *)rex->keys);
    std::fread(rex->keys, 1, rex->n_buckets * sizeof(*rex->keys), fp);
    LOG_DEBUG("About to read into vals at %p\n", (void *)rex->vals);
    std::fread(rex->vals, 1, rex->n_buckets * sizeof(*rex->vals), fp);
    return rex;
}

template <typename T>
T *khash_load(const char *path) noexcept {
    std::FILE *fp(std::fopen(path, "rb"));
    T *rex(khash_load_impl<T>(fp));
    std::fclose(fp);
    return rex;
}

void kset_union(khash_t(all) *a, khash_t(all) *b) noexcept;

tax_t lca(const khash_t(p) *map, tax_t a, tax_t b) noexcept;
unsigned node_depth(khash_t(p) *map, tax_t a) noexcept;
unsigned node_dist(khash_t(p) *map, tax_t leaf, tax_t root) noexcept;
std::unordered_map<tax_t, strlist> tax2genome_map(khash_t(name) *name_map, const std::vector<std::string> &paths);
INLINE tax_t get_parent(khash_t(p) *taxmap, tax_t key) noexcept {
    // Returns maximum value if not found.
    khiter_t ki;
    return ((ki = kh_get(p, taxmap, key)) != kh_end(taxmap)) ? kh_val(taxmap, ki)
                                                             : std::numeric_limits<tax_t>::max();
}

std::map<tax_t, tax_t> build_kraken_tax(const std::string &fname);
uint32_t lca(std::map<uint32_t, uint32_t> &parent_map, uint32_t a, uint32_t b);
std::map<uint32_t, uint32_t> kh2kr(khash_t(p) *map);

bool isfile(const char *path) noexcept;
template<typename T>
std::size_t fllen(std::forward_list<T> &list) {
    size_t ret(0);
    for(auto i(list.begin()); i != list.end(); ++i) {
        ++ret;
    }
    return ret;
}

template<typename Key, typename Map>
inline bool has_key(Key &key, Map &map) {
    return map.find(key) != map.end();
}

#if 0
superkingdom
kingdom
subkingdom
superphylum
phylum
subphylum
superclass
class
subclass
infraclass
cohort
superorder
order
suborder
infraorder
parvorder
superfamily
family
subfamily
tribe
subtribe
genus
subgenus
species
subspecies
group
subgroup
varitas
forma
#endif

INLINE const char *get_lvlname(ClassLevel lvl) {return classlvl_arr[static_cast<int>(lvl) + 2];}
ClassLevel get_linelvl(const char *line);


} // namespace emp

#endif // #ifdef _EMP_UTIL_H__
