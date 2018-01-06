#ifndef _EMP_UTIL_H__
#define _EMP_UTIL_H__
#include <algorithm>
#include <chrono>
#include <cinttypes>
#include <cstdint>
#include <cstdio>
#include <cstdlib>
#include <forward_list>
#include <limits>
#include <map>
#include <type_traits>
#include <unordered_map>
#include <unordered_set>
#include <vector>
#include "kspp/ks.h"
#include "klib/kstring.h"
#include "khash64.h"
#include "lib/logutil.h"
#include "lib/sample_gen.h"
#include "lib/rand.h"
#include "lib/lazyvec.h"

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

#ifdef USE_PDQSORT
# include "pdqsort/pdqsort.h"
# ifndef SORT
#  define SORT ::pdq::sort
# endif
# define SORT_BRANCHLESS ::pdq::sort_branchless
#else
# ifndef SORT
#  define SORT ::std::sort
# endif
# define SORT_BRANCHLESS ::std::sort
#endif
// SORT_BRANCHLESS is a lie for std::sort.

#ifndef HAS_KPUTUW__
#define kputuw_ kputuw
#endif

#define TIME_CODE(code, name) do             \
{                                            \
    auto i_##_name(::std::chrono::system_clock::now());\
    { code }                                 \
    auto j_##_name(::std::chrono::system_clock::now());\
    fprintf(stderr, "Task %s took %lfs\n", name, std::chrono::duration<double>(j_##_name - i_##_name).count());\
} while(0)

namespace emp {

using u32 = std::uint32_t;
using tax_t = u32;
using u64 = std::uint64_t;
using i32 = std::int32_t;
using i64 = std::int64_t;
using std::size_t;
using u8 = std::uint8_t;
using std::forward_list;
using std::unordered_map;
using strlist = forward_list<std::string>;
using cpslist = forward_list<std::string*>;
using bitvec_t = std::vector<u64>;

class Timer {
    using TpType = std::chrono::system_clock::time_point;
    std::string name_;
    TpType start_, stop_;
public:
    Timer(std::string &&name=""): name_{name}, start_(std::chrono::system_clock::now()) {}
    void stop() {stop_ = std::chrono::system_clock::now();}
    void restart() {start_ = std::chrono::system_clock::now();}
    void report() {std::cerr << "Took " << std::chrono::duration<double>(stop_ - start_).count() << "s for task '" << name_ << "'\n";}
    ~Timer() {stop(); /* hammertime */ report();}
    void rename(const char *name) {name_ = name;}
};

KHASH_SET_INIT_INT64(all)
KHASH_MAP_INIT_INT64(c, tax_t)
KHASH_MAP_INIT_INT64(64, u64)
KHASH_MAP_INIT_INT(p, tax_t)
KHASH_MAP_INIT_STR(name, tax_t)

std::string rand_string(size_t n);
size_t count_lines(const char *fn) noexcept;
khash_t(name) *build_name_hash(const char *fn) noexcept;
void destroy_name_hash(khash_t(name) *hash) noexcept;
void print_name_hash(khash_t(name) *hash) noexcept;
khash_t(p) *build_parent_map(const char *fn) noexcept;

lazy::vector<u64, size_t> load_binary_kmers(const char *path);
khash_t(all) *load_binary_kmerset(const char *path);

tax_t get_max_val(const khash_t(p) *hash) noexcept;
std::unordered_map<tax_t, std::vector<tax_t>> invert_parent_map(khash_t(p) *) noexcept;
tax_t get_taxid(const char *fn, khash_t(name) *name_hash);
std::string get_firstline(const char *fn);
std::vector<std::string> get_paths(const char *path);
std::vector<tax_t> get_all_descendents(const std::unordered_map<tax_t, std::vector<tax_t>> &map, tax_t tax);
std::vector<tax_t> get_desc_lca(tax_t a, tax_t, const std::unordered_map<tax_t, std::vector<tax_t>> &parent_map);
// Resolve_tree is modified from Kraken 1 source code, which
// is MIT-licensed. https://github.com/derrickwood/kraken
tax_t resolve_tree(std::map<tax_t, tax_t> &hit_counts,
                   khash_t(p) *parent_map) noexcept;

const char *bool2str(bool val);

#ifdef roundup64
#undef roundup64
#endif

// Modified from bit-twiddling hacks to work with 64-bit integers.
template<typename T>
static INLINE auto roundup64(T x) noexcept {
    kroundup64(x);
    return x;
}

INLINE u64 rand64() noexcept {
    return rng::random_twist();
}

template <typename T>
void khash_destroy(T *map) noexcept {
    LOG_EXIT("NotImplementedError");
}

template<typename T, typename KType>
inline khint_t khash_put(T *map, KType key, int *ret);
template<> inline khint_t khash_put(khash_t(64) *map, uint64_t key, int *ret) {
    return kh_put(64, map, key, ret);
}
template<> inline khint_t khash_put(khash_t(all) *map, uint64_t key, int *ret) {
    return kh_put(all, map, key, ret);
}

template<typename T> khint_t khash_get(T *map, uint64_t key) {
    if constexpr(std::is_same_v<std::decay_t<decltype(map->keys[0])>, uint64_t>) {
        return kh_get(64, (khash_t(64) *)map, key);
    } else {
        return kh_get(c, (khash_t(c) *)map, key);
    }
}

template<>
void khash_destroy(khash_t(64) *map) noexcept;
template<>
void khash_destroy(khash_t(all) *map) noexcept;
template<>
void khash_destroy(khash_t(c) *map) noexcept;
template<>
void khash_destroy(khash_t(p) *map) noexcept;
template<>
void khash_destroy(khash_t(name) *map) noexcept;

template <typename T>
void print_khash(T *rex) noexcept {
    fprintf(stderr, "n buckets %zu, nocc %zu, size %zu, upper bound %zu.\n",
                 rex->n_buckets, rex->n_occupied, rex->size, rex->upper_bound);
}

#define __fw(item, fp) \
  fwrite(&(item), 1, sizeof(item), fp)

template<typename T>
void khash_write_impl(T *map, std::FILE *fp) noexcept {
    for(khiter_t ki(0); ki != kh_end(map); ++ki)
        if(!kh_exist(map, ki))
            kh_key(map, ki) = 0, kh_val(map, ki) = 0;
    __fw(map->n_buckets, fp);
    __fw(map->n_occupied, fp);
    __fw(map->size, fp);
    __fw(map->upper_bound, fp);
    fwrite(map->flags, __ac_fsize(map->n_buckets), sizeof(*map->flags), fp);
    fwrite(map->keys, map->n_buckets, sizeof(*map->keys), fp);
    fwrite(map->vals, map->n_buckets, sizeof(*map->vals), fp);
}

template <typename T>
size_t khash_write(T *map, const char *path) noexcept {
    std::FILE *fp(fopen(path, "wb"));
    khash_write_impl(map, fp);
    size_t ret(ftell(fp));
    fclose(fp);
    return ret;
}

#undef __fw

template <typename T>
T *khash_load_impl(std::FILE *fp) noexcept {
    T *rex((T *)std::calloc(1, sizeof(T)));
    typedef typename std::remove_pointer<decltype(rex->keys)>::type keytype_t;
    typedef typename std::remove_pointer<decltype(rex->vals)>::type valtype_t;
    fread(&rex->n_buckets, 1, sizeof(rex->n_buckets), fp);
    fread(&rex->n_occupied, 1, sizeof(rex->n_occupied), fp);
    fread(&rex->size, 1, sizeof(rex->size), fp);
    fread(&rex->upper_bound, 1, sizeof(rex->upper_bound), fp);
    LOG_DEBUG("buckets: %zu. nocc: %zu. size: %zu. ub: %zu\n", (size_t)rex->n_buckets, size_t(rex->n_occupied), size_t(rex->size), size_t(rex->upper_bound));
    rex->flags = (u32 *)std::malloc(sizeof(*rex->flags) * __ac_fsize(rex->n_buckets));
    rex->keys = (keytype_t *)std::malloc(sizeof(*rex->keys) * rex->n_buckets);
    if(!rex->keys) fprintf(stderr, "Could not allocate %zu bytes of memory (%zu GB)\n", sizeof(*rex->keys) * rex->n_buckets, sizeof(*rex->keys) * rex->n_buckets >> 30), exit(1);
    rex->vals = (valtype_t *)std::malloc(sizeof(*rex->vals) * rex->n_buckets);
    if(!rex->vals) fprintf(stderr, "Could not allocate %zu bytes of memory (%zu GB)\n", sizeof(*rex->vals) * rex->n_buckets, sizeof(*rex->vals) * rex->n_buckets >> 30), exit(1);
    LOG_DEBUG("About to read into flags at %p\n", (void *)rex->flags);
    fread(rex->flags, __ac_fsize(rex->n_buckets), sizeof(*rex->flags), fp);
    LOG_DEBUG("About to read into keys at %p\n", (void *)rex->keys);
    fread(rex->keys, 1, rex->n_buckets * sizeof(*rex->keys), fp);
    LOG_DEBUG("About to read into vals at %p\n", (void *)rex->vals);
    fread(rex->vals, 1, rex->n_buckets * sizeof(*rex->vals), fp);
    return rex;
}

template <typename T>
T *khash_load(const char *path) noexcept {
    std::FILE *fp(fopen(path, "rb"));
    T *rex(khash_load_impl<T>(fp));
    fclose(fp);
    return rex;
}

void kset_union(khash_t(all) *a, khash_t(all) *b) noexcept;

tax_t lca(const khash_t(p) *map, tax_t a, tax_t b) noexcept;
unsigned node_depth(const khash_t(p) *map, tax_t a) noexcept;
unsigned node_dist(khash_t(p) *map, tax_t leaf, tax_t root) noexcept;
std::unordered_map<tax_t, strlist> tax2genome_map(khash_t(name) *name_map, const std::vector<std::string> &paths);
INLINE tax_t get_parent(const khash_t(p) *taxmap, tax_t key) noexcept {
    // Returns maximum value if not found.
    khiter_t ki;
    return ((ki = kh_get(p, taxmap, key)) != kh_end(taxmap)) ? kh_val(taxmap, ki)
                                                             : std::numeric_limits<tax_t>::max();
}

std::map<tax_t, tax_t> build_kraken_tax(const std::string &fname);
uint32_t lca(std::map<uint32_t, uint32_t> &parent_map, uint32_t a, uint32_t b);
std::map<uint32_t, uint32_t> kh2kr(khash_t(p) *map);

inline bool isfile(const char *path) noexcept {return access(path, F_OK) != -1;}
inline bool isfile(const std::string &path) noexcept {return isfile(path.data());}

template<typename T>
size_t size(const T &container) {
    return container.size();
}
template<typename T>
size_t size(const forward_list<T> &list) {
    size_t ret(0);
    for(auto i(list.begin()); i != list.end(); ++i) {
        ++ret;
    }
    return ret;
}

template<typename Key, typename Map>
inline bool has_key(const Key &key, const Map &map) {
    return map.find(key) != map.end();
}

#if !NDEBUG
template<typename T, class Compare=std::less<typename T::value_type>>
void assert_sorted_impl(const T &container, Compare cmp=Compare{}) {
    using std::begin;
    auto it(begin(container)), it2(begin(container));
    if(it2 != end(container)) ++it2;
    while(it2 != end(container)) {
        assert(cmp(*it, *it2));
        ++it, ++it2;
    }
}
#define assert_sorted(container, ...) ::emp::assert_sorted_impl(container, ##__VA_ARGS__)
#else
#define assert_sorted(container, ...)
#endif


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

INLINE const char *get_lvlname(ClassLevel lvl) {return classlvl_arr[static_cast<int>(lvl) + LINE_LVL_OFFSET];}
ClassLevel get_linelvl(const char *line, std::string &buffer, const std::unordered_map<std::string, ClassLevel> &map);
std::vector<tax_t> get_sorted_taxes(const khash_t(p) *taxmap, const char *path);
std::unordered_map<tax_t, ClassLevel> get_tax_depths(const khash_t(p) *taxmap, const char *path);
std::unordered_map<tax_t, strlist> tax2desc_genome_map(
        const std::unordered_map<tax_t, strlist> &tx2g,
        const khash_t(p) *taxmap, const std::vector<tax_t> &taxes,
        const std::unordered_map<tax_t, ClassLevel> &lvl_map);

template<typename Container, typename=std::enable_if_t<std::is_same_v<typename Container::value_type, tax_t>>>
tax_t lca(khash_t(p) *taxmap, const Container& v) noexcept {
    if(v.size() == 0) {
        fprintf(stderr, "Warning: no elements provided. Returning 0 for lca.\n");
        return 0;
    }
    auto it(v.begin());
    auto ret(*it);
    while(++it != v.end()) {
        //LOG_DEBUG("lca before: %u. lca to put with: %u.\n", ret, *it);
        ret = lca(taxmap, ret, *it);
    }
    return ret;
}

template<typename ValType, typename SetType>
std::vector<ValType> vector_set_filter(const std::vector<ValType> &vec, const SetType &set) {
    std::vector<ValType> ret;
    typename SetType::const_iterator it;
    for(const auto el: vec) if((it = set.find(el)) != set.end()) ret.emplace_back(el);
    return ret;
}

template<typename T>
std::string bitvec2str(const T &a) {
    std::string ret;
    for(auto it(a.cbegin()), eit(a.cend()); it != eit; ++it)
        for(int j(63); j >= 0; j--)
            ret += ((*it & (1ull << j)) != 0) + '0';
    return ret;
}

static inline kstring_t *kspp2ks(ks::string &ks) {
    static_assert(sizeof(kstring_t) == sizeof(ks::string), "ks::string must have the same size.");
    return reinterpret_cast<kstring_t *>(std::addressof(ks));
}

static inline const kstring_t *kspp2ks(const ks::string &ks) {
    static_assert(sizeof(kstring_t) == sizeof(ks::string), "ks::string must have the same size.");
    return reinterpret_cast<const kstring_t *>(std::addressof(ks));
}

inline constexpr int log2_64(uint64_t value)
{
    // https://stackoverflow.com/questions/11376288/fast-computing-of-log2-for-64-bit-integers
    const int tab64[64] {
        63,  0, 58,  1, 59, 47, 53,  2,
        60, 39, 48, 27, 54, 33, 42,  3,
        61, 51, 37, 40, 49, 18, 28, 20,
        55, 30, 34, 11, 43, 14, 22,  4,
        62, 57, 46, 52, 38, 26, 32, 41,
        50, 36, 17, 19, 29, 10, 13, 21,
        56, 45, 25, 31, 35, 16,  9, 12,
        44, 24, 15,  8, 23,  7,  6,  5
    };
    value |= value >> 1;
    value |= value >> 2;
    value |= value >> 4;
    value |= value >> 8;
    value |= value >> 16;
    value |= value >> 32;
    return tab64[((uint64_t)((value - (value >> 1))*0x07EDD5E59A4E28C2)) >> 58];
}
INLINE double kmer_entropy(uint64_t kmer, unsigned k) {
    if(kmer == UINT64_C(-1) || kmer == UINT64_C(0)) return 0.; // Whooooooo
    // A better solution would be avoiding redundant work, but that only works
    // for contiguous seeds.
    unsigned cts[4]{0};
    while(k--) ++cts[kmer & 0x3], kmer >>= 2;
    const double div(1./k);
    double tmp(div * cts[0]), sum(tmp * std::log2(tmp));
    for(unsigned i(1); i < 4; ++i)
        tmp = div * cts[i], sum += tmp * std::log2(tmp);
    return tmp;
}

template<typename T>
INLINE const char *get_str(const T &str) {
    return str.data();
}
template<typename T>
INLINE char *get_str(T &str) {
    return str.data();
}
INLINE char *get_str(char *str)             {return str;}
INLINE const char *get_str(const char *str) {return str;}

} // namespace emp

#endif // #ifdef _EMP_UTIL_H__
