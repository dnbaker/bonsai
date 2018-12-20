#pragma once
#ifdef NDEBUG
#  if NDEBUG == 0
#    undef NDEBUG
#  endif
#endif
#include <algorithm>
#include <cassert>
#include <chrono>
#include <cinttypes>
#include <cstdio>
#include <cstring>
#include <forward_list>
#include <fstream>
#include <functional>
#include <limits>
#include <map>
#include <random>
#include <set>
#include <sstream>
#include <type_traits>
#include <unordered_map>
#include <unordered_set>
#include <sys/stat.h>
#include <zlib.h>

#include "kspp/ks.h"
#include "clhash/include/clhash.h"
#include "khash64.h"
#include "klib/kstring.h"
#include "kseq_declare.h"
#include "lazy/vector.h"
#include "linear/linear.h"
#include "logutil.h"
#include "popcnt.h"
#include "sample_gen.h"

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
#ifndef gzbuffer
#  if ZLIB_VERNUM < 0x1235
#    define gzbuffer(fp, size)
#  else
#    define gzbuffer(fp, size) gzbuffer(fp, size)
#  endif
#endif


#ifdef USE_PDQSORT
# include "pdqsort/pdqsort.h"
# ifndef SORT
#  define SORT pdqsort
# endif
# define SORT_BRANCHLESS ::pdqsort_branchless
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

#define STLFREE(x) do {decltype(x) tmp; using std::swap; swap(tmp, x);} while(0)

#ifndef XOR_MASK
#    define XOR_MASK 0xe37e28c4271b5a2dULL
#endif

#define TIME_CODE(code, name) do             \
{                                            \
    auto i_##_name(::std::chrono::system_clock::now());\
    { code }                                 \
    auto j_##_name(::std::chrono::system_clock::now());\
    fprintf(stderr, "Task %s took %lfs\n", name, std::chrono::duration<double>(j_##_name - i_##_name).count());\
} while(0)

namespace bns {
using namespace std::literals;

struct ForPool {
    void *fp_;
    ForPool(int nthreads): fp_(kt_forpool_init(nthreads)) {}
    void forpool(void (*func)(void*,long,int), void *data, long n) {
        kt_forpool(fp_, func, data, n);
    }
    ~ForPool() {
        kt_forpool_destroy(fp_);
    }
};

using MainFnPtr = int (*) (int, char **);

using i16 = std::uint16_t;
using i32 = std::int32_t;
using i64 = std::int64_t;
using i8  = std::uint8_t;
using u16 = std::uint16_t;
using u32 = std::uint32_t;
using u64 = std::uint64_t;
using u8  = std::uint8_t;
using std::size_t;
using tax_t = u32;
using pop::bitvec_t;

template<typename MutexType>
struct LockSmith {
    // Simple lock-holder to avoid writing to the same file twice.
    MutexType &m_;
    LockSmith(MutexType &m): m_(m) {m_.lock();}
    ~LockSmith() {m_.unlock();}
};

class Timer {
    using TpType = std::chrono::system_clock::time_point;
    std::string name_;
    TpType start_, stop_;
public:
    Timer(std::string &&name=""): name_{std::move(name)}, start_(std::chrono::system_clock::now()) {}
    void stop() {stop_ = std::chrono::system_clock::now();}
    void restart() {start_ = std::chrono::system_clock::now();}
    void report() {std::cerr << "Took " << std::chrono::duration_cast<std::chrono::nanoseconds>(stop_ - start_).count() << "ns for task '" << name_ << "'\n";}
    ~Timer() {stop(); /* hammertime */ report();}
    void rename(const char *name) {name_ = name;}
};

namespace {
template<typename T> class TD;
}

KHASH_SET_INIT_INT64(all)
KHASH_MAP_INIT_INT64(c, tax_t)
KHASH_MAP_INIT_INT64(64, u64)
KHASH_MAP_INIT_INT(p, tax_t)
KHASH_MAP_INIT_STR(name, tax_t)

// Resolve_tree is modified from Kraken 1 source code, which
// is MIT-licensed. https://github.com/derrickwood/kraken

#ifdef roundup64
#undef roundup64
#endif

// Modified from bit-twiddling hacks to work with 64-bit integers.
template<typename T>
static INLINE auto roundup64(T x) noexcept {
    --x;
    x |= x>>1;
    x |= x>>2;
    x |= x>>4;
    if constexpr(sizeof(x) >= 2)
        x |= x>>8;
    if constexpr(sizeof(x) >= 4)
        x |= x>>16;
    if constexpr(sizeof(x) >= 8)
        x |= x>>32;
    ++x;
    return x;
}

template <typename T>
static void khash_destroy(T *map) noexcept {
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

template<typename T, typename KType> khint_t khash_get(T *map, KType key) {
#if __cplusplus < 201703L
    if (std::is_same<std::decay_t<KType>, uint64_t>::value) {
        return kh_get(64, (khash_t(64) *)map, (uint64_t)key);
    } else if (std::is_same<KType, const char *>::value) {
        return kh_get(name, (khash_t(name) *)map, (const char *)key);
    } else {
        return kh_get(c, (khash_t(c) *)map, (uint32_t)key);
    }
#else
    if constexpr(std::is_same_v<std::decay_t<KType>, uint64_t>) {
        return kh_get(64, (khash_t(64) *)map, (uint64_t)key);
    } else if constexpr(std::is_same_v<KType, const char *>) {
        return kh_get(name, (khash_t(name) *)map, (const char *)key);
    } else {
        return kh_get(c, (khash_t(c) *)map, (uint32_t)key);
    }
#endif
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

#define KHR(x) khashraii_##x##_t

#define KHRAII_DEC(x)\
\
class KHR(x) {\
    khash_t(x) h_;\
public:\
    KHR(x)() {std::memset(this, 0, sizeof(*this));}\
    void resize(size_t n) {kh_resize(x, &h_, n);} \
    KHR(x) (size_t n): KHR(x)() {resize(n);} \
    operator khash_t(x) *() {return &h_;}\
    operator const khash_t(x) *() const {return &h_;}\
    khash_t(x) *operator->() {return &h_;}\
    const khash_t(x) *operator->() const {return &h_;}\
    KHR(x)(const KHR(x) &other): KHR(x)(other->n_buckets) {\
        std::memcpy(h_.keys, other->keys, sizeof(*h_.keys) * other->n_buckets);\
        std::memcpy(h_.vals, other->vals, sizeof(*h_.vals) * other->n_buckets);\
        std::memcpy(h_.flags, other->flags, sizeof(*h_.flags) * __ac_fsize(other->n_buckets));\
    }\
    KHR(x)(KHR(x) &&other) {\
        std::memcpy(this, &other, sizeof(*this));\
        other.h_.keys = nullptr;\
        other.h_.vals = nullptr;\
        other.h_.flags = nullptr;\
        other.h_.n_buckets = other.h_.n_occupied = 0;\
    }\
    ~KHR(x)() {\
        std::free(h_.keys);\
        std::free(h_.vals);\
        std::free(h_.flags);\
    }\
};

KHRAII_DEC(64)
KHRAII_DEC(all)
KHRAII_DEC(p)
KHRAII_DEC(c)

template <typename T>
void print_khash(T *rex) noexcept {
    fprintf(stderr, "n buckets %zu, nocc %zu, size %zu, upper bound %zu.\n",
                 rex->n_buckets, rex->n_occupied, rex->size, rex->upper_bound);
}


#define __fw(item, fn) ::write(fn, static_cast<const void *>(std::addressof(item)), sizeof(item))
template<typename T>
size_t khash_write_impl(const T *map, const int fn) noexcept {
    for(khiter_t ki(0); ki != kh_end(map); ++ki)
        if(!kh_exist(map, ki))
            kh_key(map, ki) = 0, kh_val(map, ki) = 0;
    size_t ret = __fw(map->n_buckets, fn);
    ret += __fw(map->n_occupied, fn);
    ret += __fw(map->size, fn);
    ret += __fw(map->upper_bound, fn);
    ret += ::write(fn, map->flags, __ac_fsize(map->n_buckets) * sizeof(*map->flags));
    ret += ::write(fn, map->keys, map->n_buckets * sizeof(*map->keys));
    ret += ::write(fn, map->vals, map->n_buckets * sizeof(*map->vals));
    return ret;
}
#undef __fw
template<typename T>
size_t khash_write_impl(const T *map, std::FILE *fp) noexcept {
    std::fflush(fp);
    const auto ret = khash_write_impl<T>(map, fileno(fp));
    std::fflush(fp);
    return ret;
}
#define __gz(item, fp) gzwrite(fp, static_cast<const void *>(&(item)), sizeof(item))
template<typename T>
size_t khash_write_impl(const T *map, gzFile fp) noexcept {
    for(khiter_t ki(0); ki != kh_end(map); ++ki)
        if(!kh_exist(map, ki))
            kh_key(map, ki) = 0, kh_val(map, ki) = 0;
    size_t ret = __gz(map->n_buckets, fp);
    ret += __gz(map->n_occupied, fp);
    ret += __gz(map->size, fp);
    ret += __gz(map->upper_bound, fp);
    ret += gzwrite(fp, static_cast<const void *>(map->flags), __ac_fsize(map->n_buckets) * sizeof(*map->flags));
    ret += gzwrite(fp, static_cast<const void *>(map->keys), map->n_buckets * sizeof(*map->keys));
    ret += gzwrite(fp, static_cast<const void *>(map->vals), map->n_buckets * sizeof(*map->vals));
    return ret;
}
#undef __gz

template <typename T>
size_t khash_write(const T *map, const char *path, bool write_gz=false) noexcept {
    if(write_gz) {
        gzFile fp = gzopen(path, "wb");
        const auto ret = khash_write_impl(map, fp);
        gzclose(fp);
        return ret;
    }
    std::FILE *fp(fopen(path, "wb"));
    const auto ret = khash_write_impl(map, fp);
    fclose(fp);
    return ret;
}


template <typename T>
T *khash_load_impl(const int fn) noexcept {
    T *rex((T *)std::calloc(1, sizeof(T)));
    using keytype_t = std::remove_pointer_t<decltype(rex->keys)>;
    using valtype_t = std::remove_pointer_t<decltype(rex->vals)>;
    ::read(fn, &rex->n_buckets, sizeof(rex->n_buckets));
    ::read(fn, &rex->n_occupied, sizeof(rex->n_occupied));
    ::read(fn, &rex->size, sizeof(rex->size));
    ::read(fn, &rex->upper_bound, sizeof(rex->upper_bound));
    rex->flags = (u32 *)std::malloc(sizeof(*rex->flags) * __ac_fsize(rex->n_buckets));
    if(!rex->flags) fprintf(stderr, "Could not allocate %zu bytes of memory (%zu GB)\n", (sizeof(*rex->flags) * __ac_fsize(rex->n_buckets)), (sizeof(*rex->flags) * __ac_fsize(rex->n_buckets)) >> 30), exit(1);
    rex->keys = (keytype_t *)std::malloc(sizeof(*rex->keys) * rex->n_buckets);
    if(!rex->keys) fprintf(stderr, "Could not allocate %zu bytes of memory (%zu GB)\n", sizeof(*rex->keys) * rex->n_buckets, sizeof(*rex->keys) * rex->n_buckets >> 30), exit(1);
    rex->vals = (valtype_t *)std::malloc(sizeof(*rex->vals) * rex->n_buckets);
    if(!rex->vals) fprintf(stderr, "Could not allocate %zu bytes of memory (%zu GB)\n", sizeof(*rex->vals) * rex->n_buckets, sizeof(*rex->vals) * rex->n_buckets >> 30), exit(1);
    ::read(fn, rex->flags, __ac_fsize(rex->n_buckets) * sizeof(*rex->flags));
    ::read(fn, rex->keys, rex->n_buckets * sizeof(*rex->keys));
    ::read(fn, rex->vals, rex->n_buckets * sizeof(*rex->vals));
    return rex;
}

template <typename T>
T *khash_load_impl(std::FILE *fp) noexcept {
    std::fflush(fp);
    auto rex = khash_load_impl<T>(fileno(fp));
    std::fflush(fp);
    return rex;
}

template <typename T>
T *khash_load(const char *path) noexcept {
    std::FILE *fp(fopen(path, "rb"));
    T *rex(khash_load_impl<T>(fp));
    fclose(fp);
    return rex;
}

INLINE tax_t get_parent(const khash_t(p) *taxmap, tax_t key) noexcept {
    // Returns maximum value if not found.
    khiter_t ki;
    return ((ki = kh_get(p, taxmap, key)) != kh_end(taxmap)) ? kh_val(taxmap, ki)
                                                             : std::numeric_limits<tax_t>::max();
}

inline bool isfile(const char *path) noexcept {return access(path, F_OK) != -1;}
inline bool isfile(const std::string &path) noexcept {return isfile(path.data());}

template<typename T>
size_t size(const T &container) {return container.size();}
template<typename T>
size_t size(const std::forward_list<T> &list) {
    size_t ret(0);
    for(auto i(list.begin()); i != list.end(); ++i) ++ret;
    return ret;
}

template<typename Key, typename Map>
inline bool has_key(const Key &key, const Map &map) {
    return map.find(key) != map.end();
}

INLINE u32 nuccount(u64 kmer, unsigned k) {
    // Modified from Bowtie2.
    // Returns counts in the 4 different 8-bit registers of a 32-bit integer.
    const uint64_t COUNT_MASK = (0xFFFFFFFFFFFFFFFF >> (64 - 2 * k));
    static const uint64_t c_table[4] = {
        0xffffffffffffffff, // b1111111111111111111111111111111111111111111111111111111111111111
        0xaaaaaaaaaaaaaaaa, // b1010101010101010101010101010101010101010101010101010101010101010
        0x5555555555555555, // b0101010101010101010101010101010101010101010101010101010101010101
        0x0000000000000000  // b0000000000000000000000000000000000000000000000000000000000000000
    };
    uint64_t c0 = c_table[0];
    uint64_t x0 = kmer ^ c0;
    uint64_t x1 = (x0 >> 1);
    uint64_t x2 = x1 & UINT64_C(0x5555555555555555);
    uint64_t x3 = x0 & x2;
    x3 &= COUNT_MASK;
    auto tmp = pop::popcount(x3); // because __builtin_popcountll returns a 32-bit element.
    u32 ret = (tmp <<= 24);

    c0 = c_table[1];
    x0 = kmer ^ c0;
    x1 = (x0 >> 1);
    x2 = x1 & UINT64_C(0x5555555555555555);
    x3 = x0 & x2;
    x3 &= COUNT_MASK;
    tmp = pop::popcount(x3);
    ret |= (tmp <<= 16);

    c0 = c_table[2];
    x0 = kmer ^ c0;
    x1 = (x0 >> 1);
    x2 = x1 & UINT64_C(0x5555555555555555);
    x3 = x0 & x2;
    x3 &= COUNT_MASK;
    tmp = pop::popcount(x3);
    ret |= (tmp <<= 8);

    c0 = c_table[3];
    x0 = kmer ^ c0;
    x1 = (x0 >> 1);
    x2 = x1 & UINT64_C(0x5555555555555555);
    x3 = x0 & x2;
    x3 &= COUNT_MASK;
    return ret |= (tmp = pop::popcount(x3));
}

INLINE const char *get_lvlname(ClassLevel lvl) {return classlvl_arr[static_cast<int>(lvl) + LINE_LVL_OFFSET];}

template<typename Container, typename=typename std::enable_if<std::is_same<typename Container::value_type, tax_t>::value>::type>
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
        for(int j(63); j >= 0; ret += ((*it & (1ull << j--)) != 0) + '0');
    return ret;
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
    // This could be replaced with a __builtin_clz
    return tab64[((uint64_t)((value - (value >> 1))*0x07EDD5E59A4E28C2)) >> 58];
}

static inline void kseq_assign(kseq_t *ks, gzFile fp) {
    if(!ks->f) {
        ks->f = (kstream_t*)calloc(1, sizeof(kstream_t));
        ks->f->buf = (unsigned char*)malloc(KSTREAM_SIZE);
    } else ks->f->is_eof = ks->f->begin = ks->f->end = 0;
    ks->f->f = fp;
}

static inline kseq_t kseq_init_stack() {
    kseq_t ret;
    std::memset(&ret, 0, sizeof(ret));
    return ret;
}

static INLINE void kseq_destroy_stack(kseq_t &ks) {
    free(ks.name.s); free(ks.comment.s); free(ks.seq.s); free(ks.qual.s);
    ks_destroy(ks.f);
    std::memset(&ks, 0, sizeof(ks));
}


INLINE double kmer_entropy(uint64_t kmer, unsigned k) {
    const u32 counts(nuccount(kmer, k));
    const double div(1./k);
    double tmp(div * (counts >> 24)),    sum(tmp * std::log2(tmp));
    tmp = div * ((counts >> 16) & 0xFF), sum += tmp * std::log2(tmp);
    tmp = div * ((counts >> 8) & 0xFF),  sum += tmp * std::log2(tmp);
    return tmp = div * (counts & 0xFF),  sum += tmp * std::log2(tmp);
}

template<typename T> INLINE const char *get_cstr(const T &str) {return str.data();}
template<typename T> INLINE char *get_cstr(T &str)             {return str.data();}
INLINE char *get_cstr(char *str)             {return str;}
INLINE const char *get_cstr(const char *str) {return str;}
INLINE char *strchrnul(char *str, int c) {
    while(*str && *str != c) ++str;
    return str;
}

namespace detail {
    static const std::unordered_map<int, const char *> zerr2str {
        {Z_OK, "Z_OK"},
        {Z_STREAM_END, "Z_STREAM_END"},
        {Z_STREAM_END, "Z_STREAM_END"},
        {Z_NEED_DICT, "Z_NEED_DICT"},
        {Z_ERRNO, "Z_ERRNO"},
        {Z_STREAM_ERROR, "Z_STREAM_ERROR"},
        {Z_DATA_ERROR, "Z_DATA_ERROR"},
        {Z_MEM_ERROR, "Z_MEM_ERROR"},
        {Z_BUF_ERROR, "Z_BUF_ERROR"},
        {Z_VERSION_ERROR, "Z_VERSION_ERROR"}
    };
}

#ifndef DO_DUFF
#define DO_DUFF(len, ITER) \
    do { \
        if(len) {\
            std::uint64_t loop = (len + 7) >> 3;\
            switch(len & 7) {\
                case 0: do {\
                    ITER; [[fallthrough]];\
                    case 7: ITER; [[fallthrough]]; case 6: ITER; [[fallthrough]]; case 5: ITER; [[fallthrough]];\
                    case 4: ITER; [[fallthrough]]; case 3: ITER; [[fallthrough]]; case 2: ITER; [[fallthrough]]; case 1: ITER;\
                } while (--loop);\
            }\
        }\
    } while(0)
#endif

class zlib_error: public std::runtime_error {
public:
    zlib_error(int c, const char *fn=nullptr):
        std::runtime_error(ks::sprintf("zlib error code %u (%s) accessing file %s.", c,
                                       detail::zerr2str.at(c), fn ? fn: "<no file provided>").data())
    {}
};

#ifndef RUNTIME_ERROR
#define RUNTIME_ERROR(msg) \
        throw std::runtime_error(std::string("[") + __FILE__ + ':' + __PRETTY_FUNCTION__ + std::to_string(__LINE__) + "] " + msg)
#endif

struct KSeqBufferHolder {
    std::vector<kseq_t> kseqs_;
    KSeqBufferHolder(size_t n) {
        while(kseqs_.size() < n) kseqs_.emplace_back(kseq_init_stack());
    }
    ~KSeqBufferHolder() {
        this->free();
    }
    kseq_t &operator[](size_t index) {
        return kseqs_[index];
    }
    const kseq_t &operator[](size_t index) const {
        return kseqs_[index];
    }
    auto data() {return kseqs_.data();}
    void free() {
        for(auto &ks: kseqs_) kseq_destroy_stack(ks);
    }
};

static size_t count_lines(const char *fn) noexcept {
    std::ifstream is(fn);
    if(!is.good()) LOG_EXIT("Could not open file at %s\n", fn);
    std::string line;
    size_t n(0);
    while(std::getline(is, line)) ++n;
    return n;
}

static std::unordered_map<tax_t, std::vector<tax_t>> invert_parent_map(const khash_t(p) *tax) noexcept {
    std::unordered_map<tax_t, std::vector<tax_t>> ret;
    typename std::unordered_map<tax_t, std::vector<tax_t>>::iterator m;
    for(khiter_t ki(0); ki < kh_end(tax); ++ki) {
        if(!kh_exist(tax, ki)) continue;
        if((m = ret.find(kh_val(tax, ki))) == ret.end()) m = ret.emplace(kh_val(tax, ki), std::vector<tax_t>{kh_key(tax, ki)}).first;
        else                                             m->second.emplace_back(kh_key(tax, ki));
    }
#if !NDEBUG
    for(auto &pair: ret) {
        for(auto child: pair.second) assert(get_parent(tax, child) == pair.first);
    }
#endif
    return ret;
}

static std::vector<tax_t> get_all_descendents(const std::unordered_map<tax_t, std::vector<tax_t>> &map, tax_t tax) {
    std::vector<tax_t> ret;
    auto m(map.find(tax));
    if(m != map.end()) {
        LOG_DEBUG("Found a match for parent %u with %zu children.\n", tax, m->second.size());
        for(tax_t newtax: m->second) {
            assert(newtax != tax);
            ret.push_back(newtax);
            auto desc(get_all_descendents(map, newtax));
            for(auto d: desc) std::fprintf(stderr, "Descendent %u of total %zu descendences of %u\n", d, desc.size(), tax);
            ret.insert(ret.end(), desc.begin(), desc.end());
        }
    }
    LOG_DEBUG("Returning descendents. Size: %zu\n", ret.size());
    return ret;
}

// Rewritten from Kraken's source code.
// Consider rewriting to use kbtree instead of std::map.
// lh3's benchmarks indicate that his is only as fast as std::map,
// though it uses less than half the memory. The thing is, taxonomic trees
// aren't that deep.
// I don't think we'd gain anything practical with that.
static tax_t lca(const khash_t(p) *map, tax_t a, tax_t b) noexcept {
    // Use linear::set to use a flat map trees for this (very) small set rather than hash table.
    // linear sets will be faster up to ~100 elements, and the taxonomic tree isn't that deep.
    if(unlikely(map == nullptr)) {
        std::fprintf(stderr, "null taxonomy.\n");
        std::exit(EXIT_FAILURE);
    }
    linear::set<tax_t> nodes;
    if(a == b) return a;
    if(b == 0) return a;
    if(a == 0) return b;
    khint_t ki;
    while(a) {
        nodes.insert(a);
        if((ki = kh_get(p, map, a)) == kh_end(map)) {
            std::fprintf(stderr, "Missing taxid %u. Returning node b (%u)!\n", a, b);
            return (tax_t)-1;
        }
        a = kh_val(map, ki);
    }
    while(b) {
        if(nodes.find(b) != nodes.end()) return b;
        if((ki = kh_get(p, map, b)) == kh_end(map)) {
            std::fprintf(stderr, "Missing taxid %u. Returning node b (%u)!\n", b, a);
            return (tax_t)-1;
        }
        b = kh_val(map, ki);
    }
    return 1;
}

static unsigned node_dist(const khash_t(p) *map, tax_t leaf, tax_t root) noexcept {
    unsigned ret(0);
    khint_t ki;
    while(leaf) {
        if((ki = kh_get(p, map, leaf)) == kh_end(map)) {
            std::fprintf(stderr, "Tax ID %u missing. Abort!\n", leaf);
            std::exit(1);
        }
        ++ret;
        if((leaf = kh_val(map, ki)) == root) return ret;
    }
    LOG_EXIT("leaf %u is not a child of root %u\n", leaf, root);
    return ret;
}
static unsigned node_depth(const khash_t(p) *map, tax_t a) noexcept {
    unsigned ret(0);
    khint_t ki;
    while(a) {
        if((ki = kh_get(p, map, a)) == kh_end(map)) {
            std::fprintf(stderr, "Tax ID %u missing. Abort!\n", a);
            std::exit(1);
        }
        a = kh_val(map, ki);
        ++ret;
    }
    return ret;
}

static khash_t(name) *build_name_hash(const char *fn) noexcept {
    size_t bufsz(2048), namelen;
    char *buf((char *)std::malloc(bufsz));
    ssize_t len;
    std::FILE *fp(std::fopen(fn, "r"));
    khash_t(name) *ret(kh_init(name));
    kh_resize(name, ret, count_lines(fn));
    char *p;
    int khr;
    khint_t ki;
    while((len = getline(&buf, &bufsz, fp)) >= 0) {
        switch(*buf) case '\0': case '\n': case '#': continue;
        p = ::bns::strchrnul(buf, '\t');
        ki = kh_put(name, ret, buf, &khr);
        if(khr == 0) { // Key already present.
            LOG_INFO("Key %s already present. Updating value from "
                     "%i to %s.\tNote: if you have performed TaxonomyReformation, this is an error.\n", kh_key(ret, ki), kh_val(ret, ki), p + 1);
        } else {
            namelen = p - buf;
            char *tmp = static_cast<char *>(std::malloc(namelen + 1));
            std::memcpy(tmp, buf, namelen);
            tmp[namelen] = '\0';
            kh_key(ret, ki) = tmp;
        }
        kh_val(ret, ki) = std::atoi(p + 1);
    }
    std::fclose(fp);
    std::free(buf);
    return ret;
}

static void destroy_name_hash(khash_t(name) *hash) noexcept {
    if(hash == nullptr) return;
    for(khint_t ki(kh_begin(hash)); ki != kh_end(hash); ++ki)
        if(kh_exist(hash, ki))
            std::free((void *)kh_key(hash, ki));
    kh_destroy(name, hash);
}

static std::map<tax_t, tax_t> build_kraken_tax(const std::string &fname) {
    const char *fn(fname.data());
    std::FILE *fp(std::fopen(fn, "r"));
    std::map<tax_t, tax_t> ret;
    size_t bufsz = 4096;
    char *buf((char *)std::malloc(bufsz));
    ssize_t len;
    tax_t t1, t2;
    while((len = getline(&buf, &bufsz, fp)) >= 0) {
        switch(*buf) case '\n': case '0': case '#': continue;
        t1 = atoi(buf);
        t2 = atoi(strchr(buf, '|') + 2);
        ret[t1] = t2;
    }
    ret[1] = 0;
    std::fclose(fp);
    std::free(buf);
    return ret;
}

static tax_t lca(const std::map<tax_t, tax_t> &parent_map, tax_t a, tax_t b)
{
  if(a == 0 || b == 0) return a ? a : b;

  linear::set<uint32_t> a_path;
  while (a) a_path.insert(a), a = parent_map.find(a)->second;
  while (b) {
    if(a_path.find(b) != a_path.end())
      return b;
    b = parent_map.find(b)->second;
  }
  return 1;
}

static khash_t(p) *build_parent_map(const char *fn) {
    std::ifstream is(fn);
    khash_t(p) *ret(kh_init(p));
    khint_t ki;
    int khr;
    std::string line;
    char *p;
    while(std::getline(is, line)) {
        switch(line[0]) case '\n': case '\0': case '#': continue;
        ki = kh_put(p, ret, std::atoi(line.data()), &khr);
        kh_val(ret, ki) = (p = std::strchr(line.data(), '|')) ? std::atoi(p + 2)
                                                              : tax_t(-1);
        if(kh_val(ret, ki) == tax_t(-1)) LOG_WARNING("Malformed line: %s", line.data());
    ki = kh_put(p, ret, 1, &khr);
    kh_val(ret, ki) = 0; // Root of the tree.
    if(kh_size(ret) < 2) RUNTIME_ERROR(std::string("Failed to create taxmap from ") + fn);
    LOG_DEBUG("Built parent map of size %zu from path %s\n", kh_size(ret), fn);
    return ret;
}

static void kset_union(khash_t(all) *a, khash_t(all) *b) noexcept {
    khint_t ki2;
    int khr;
    for(ki2 = kh_begin(b); ki2 != kh_end(b); ++ki2)
        if(kh_exist(b, ki2))
            kh_put(all, a, kh_key(b, ki2), &khr);
}

 // Tree resolution: take all hit taxa (plus ancestors), then
 // return leaf of highest weighted leaf-to-root path.
 // Taken, modified from Kraken with new, faster containers.
#if !NDEBUG
static tax_t resolve_tree(const std::map<tax_t, tax_t> &hit_counts,
                          const khash_t(p) *parent_map) noexcept
{
  linear::set<tax_t> max_taxa;
  tax_t max_taxon(0), max_score(0);

  // Sum each taxon's LTR path
  for(auto it(hit_counts.cbegin()), e(hit_counts.cend()); it != e; ++it) {
    tax_t taxon(it->first), node(taxon), score(0);
    // Instead of while node > 0
    while(node)
        score += hit_counts.at(node), node = kh_val(parent_map, kh_get(p, parent_map, node));
    if(score > max_score) {
      max_taxa.clear();
      max_score = score;
      max_taxon = taxon;
      //LOG_DEBUG("max score: %u. max taxon: %u\n", max_score, max_taxon);
    } else if(score == max_score) {
      if(max_taxa.empty()) // Is this check needed?
        max_taxa.insert(max_taxon);
      max_taxa.insert(taxon);
    }
  }
  // If two LTR paths are tied for max, return LCA of all
  if(max_taxa.size()) {
    auto sit(max_taxa.begin());
    for(max_taxon = *sit++;sit != max_taxa.end(); max_taxon = lca(parent_map, max_taxon, *sit++));
  }
  return max_taxon;
}
#endif

static tax_t resolve_tree(const linear::counter<tax_t, u16> &hit_counts,
                          const khash_t(p) *parent_map) noexcept
{
  linear::set<tax_t> max_taxa;
  tax_t max_taxon(0), max_score(0);

  // Sum each taxon's LTR path
  for(unsigned i(0); i < hit_counts.size(); ++i) {
    tax_t taxon(hit_counts.keys()[i]), node(taxon), score(0);
    // Instead of while node > 0
    while(node) {
        score += hit_counts.count(node);
        node = kh_val(parent_map, kh_get(p, parent_map, node));
    }
    // In C++20, we can use the spaceship operator!
    if(score > max_score) {
      max_taxa.clear();
      max_score = score;
      max_taxon = taxon;
    } else if(score == max_score) {
      if(max_taxa.empty()) // Is this check needed?
        max_taxa.insert(max_taxon);
      max_taxa.insert(taxon);
    }
  }
  // If two LTR paths are tied for max, return LCA of all
  if(max_taxa.size()) {
    auto sit(max_taxa.begin());
    for(max_taxon = *sit++;sit != max_taxa.end(); max_taxon = lca(parent_map, max_taxon, *sit++));
  }
#if !NDEBUG
  std::map<tax_t, tax_t> cpy;
  for(size_t i(0); i < hit_counts.keys().size(); ++i) {
    cpy.emplace(hit_counts.keys()[i], hit_counts.vals()[i]);
  }
  assert(max_taxon == resolve_tree(cpy, parent_map));
#endif
  return max_taxon;
}


static std::string rand_string(size_t n) {
    std::string ret;
    ret.reserve(n);
    static const char set[] = "abcdefghijklmnopqrstuvwxyz123456";
    while(ret.size() < n) ret += set[std::rand() & 31];
    assert(ret.size() == n);
    return ret;
}

static std::string get_firstline(const char *fn) {
    gzFile fp(gzopen(fn, "rb"));
    if(fp == nullptr) LOG_EXIT("Could not read from file %s\n", fn);
    static const size_t bufsz(2048);
    char buf[bufsz];
    std::string ret(gzgets(fp, buf, bufsz));
    if(ret.back() != '\n')
        RUNTIME_ERROR(std::string("[E:get_firstline] line from ") +
                                  fn + " did not fit in buffer of size " +
                                  std::to_string(bufsz) +
                                  ". Try recompiling with a larger "
                                  "buffer or rewrite get_firstline.");
    ret.pop_back();
    gzclose(fp);
    return ret;
}

static tax_t get_taxid(const char *fn, const khash_t(name) *name_hash) {
    gzFile fp(gzopen(fn, "rb"));
    if(fp == nullptr) LOG_EXIT("Could not read from file %s\n", fn);
    static const size_t bufsz(2048);
    khint_t ki;
    char buf[bufsz];
    char *line(gzgets(fp, buf, bufsz));
    if(line == nullptr) {
        int err;
        LOG_INFO("zlib error: %s\n", gzerror(fp, &err));
        throw zlib_error(err, fn);
    }
    char *p(++line);
#ifdef SYNTHETIC_GENOME_EXPERIMENTS
    const tax_t ret(std::atoi(p));
#else
    if(std::strchr(p, '|')) {
        p = std::strrchr(p, '|');
        while(*--p != '|');
        char *q(std::strchr(++p, '|'));
        *q = 0;
        line = p;
        //if(strchr(p, '.')) *strchr(p, '.') = 0;
    } else {
        while(!std::isspace(*p)) ++p;
        *p = 0;
    }
    const tax_t ret(unlikely((ki = kh_get(name, name_hash, line)) == kh_end(name_hash)) ? UINT32_C(1) : kh_val(name_hash, ki));
#endif
    gzclose(fp);
    return ret;
}

static std::map<uint32_t, uint32_t> kh2kr(khash_t(p) *map) {
    std::map<uint32_t, uint32_t> ret;
    if(map)
        for(khiter_t ki(0); ki != kh_end(map); ++ki)
            if(kh_exist(map, ki))
                ret[kh_key(map, ki)] = kh_val(map, ki);
    return ret;
}


static std::unordered_map<tax_t, std::forward_list<std::string>> tax2genome_map(khash_t(name) *name_map, const std::vector<std::string> &paths) {
    tax_t taxid;
    std::unordered_map<tax_t, std::forward_list<std::string>> ret;
    typename std::unordered_map<tax_t, std::forward_list<std::string>>::iterator m;
    ret.reserve(paths.size());
#if !NDEBUG
    ks::string ks;
#endif
    for(const auto &path: paths) {
        if((taxid = get_taxid(path.data(), name_map)) == UINT32_C(-1)) continue;
        if((m = ret.find(taxid)) == ret.end()) m = ret.emplace(taxid, std::forward_list<std::string>{path}).first;
        else if(std::find(m->second.begin(), m->second.end(), path) == m->second.end()) m->second.push_front(path);
#if !NDEBUG
        ks.clear();
        ks.sprintf("Path %s has taxid %u.\n", path.data(), taxid);
        for(const std::string &p: m->second) {
            ks.sprintf("Element: \"%s\",", p.data());
        }
        ks.back() = '\n';
        std::cerr << ks.data();
        if((ret.size() & (ret.size() - 1)) == 0) std::cerr << "Size of tax2genome_map is now " << ret.size() << '\n';
#endif
    }
    return ret;
}


static std::unordered_map<tax_t, std::set<tax_t>> make_ptc_map(
        const khash_t(p) *taxmap, const std::vector<tax_t> &taxes,
        const std::unordered_map<tax_t, ClassLevel> &lvl_map) {
    std::unordered_map<tax_t, std::set<tax_t>> ret;
    ret.reserve(kh_size(taxmap));
    std::vector<tax_t> sorted_taxes = taxes;
#if !NDEBUG
    bool fail(false);
    for(const auto tax: sorted_taxes)
        if(lvl_map.find(tax) == lvl_map.end()) std::cerr << "Missing tax level for " << tax << '\n', fail = true;
    if(fail) RUNTIME_ERROR("Failed for missin tax levels.");
#endif
    pdqsort(sorted_taxes.begin(), sorted_taxes.end(), [&lvl_map](const tax_t a, const tax_t b) {
            try {
               return lvl_map.at(a) < lvl_map.at(b);
            } catch(std::out_of_range &ex) {
               std::cerr << "Out of range: " << ex.what() << '\n';
               if(lvl_map.find(a) == lvl_map.end()) {std::cerr << "Missing tax " << a << '\n'; throw;}
               if(lvl_map.find(b) == lvl_map.end()) {std::cerr << "Missing tax " << b << '\n'; throw;}
               RUNTIME_ERROR("ZOMG");
            }
    });
#if !NDEBUG
    for(auto it1(sorted_taxes.begin()), it2(it1 + 1); it2 != sorted_taxes.end(); ++it1, ++it2) {
        assert(lvl_map.at(*it1) <= lvl_map.at(*it2));
    }
#endif
    typename std::unordered_map<tax_t, std::set<tax_t>>::iterator it;
    khiter_t ki;
    tax_t tmptax;
    for(auto sit(sorted_taxes.begin()), eit(sorted_taxes.end()); sit != eit; ++sit) {
        const tax_t tax(*sit);
#if 0
        std::cerr << "Processing tax " << tax << '#' << static_cast<size_t>(sit - sorted_taxes.begin()) << " of " << sorted_taxes.size() << '\n';
#endif
        if((it = ret.find(tax)) == ret.end()) ret.emplace(tax, std::set<tax_t>{});
        tmptax = tax;
        while((ki = kh_get(p, taxmap, tmptax)) != kh_end(taxmap)) {
            if(kh_val(taxmap, ki) == 0) goto loop_end;
            tmptax = kh_val(taxmap, ki);
            if((it = ret.find(tmptax)) == ret.end()) ret.emplace(tmptax, std::set<tax_t>{tax});
            else                                     it->second.insert(tax);
        }
        // Only reach if taxid is missing from taxonomy file.
        RUNTIME_ERROR(std::string("Missing parent for taxid ") + std::to_string(tax));
        loop_end:
#if 0
        std::cerr << "Now finishing up for tax " << tax << '\n';
#else
        continue; // Syntactic boilerplate.
#endif
    }
    return ret;
}

static std::unordered_map<tax_t, std::forward_list<std::string>> tax2desc_genome_map(
        const std::unordered_map<tax_t, std::forward_list<std::string>> &tx2g,
        const khash_t(p) *taxmap, const std::vector<tax_t> &taxes,
        const std::unordered_map<tax_t, ClassLevel> &lvl_map) {
    std::unordered_map<tax_t, std::forward_list<std::string>> ret;
    for(const auto &pair: make_ptc_map(taxmap, taxes, lvl_map)) {
        typename std::unordered_map<tax_t, std::forward_list<std::string>>::const_iterator pit;
        std::forward_list<std::string> list;
        if((pit = tx2g.find(pair.first)) != tx2g.end()) {
            for(const auto &str: pit->second) list.push_front(str);
        }
        else {
            if(kh_get(p, taxmap, pair.first) == kh_end(taxmap)) {
                std::cerr << "No parent for node " << (int)pair.first << '\n';
                RUNTIME_ERROR(std::string("Invalid taxid ") + std::to_string(pair.first));
            } else {
#if 0
                std::cerr << "Valid tax id " << kh_key(taxmap, kh_get(p, taxmap, pair.first))
                          << " with parent " << kh_val(taxmap, kh_get(p, taxmap, pair.first))
                          << ". Do nothing.\n";
#endif
            }
        }
        for(const auto child: pair.second)
            if((pit = tx2g.find(child)) != tx2g.end())
                for(const auto &el: pit->second)
                    list.push_front(el);
        ret.emplace(pair.first, std::move(list));
    }
    return ret;
}


static inline const char *bool2str(bool val) {
    static const char * vals [] {
        "false",
        "true"
    };
    return vals[val];
}

#define _KHD(x) template<> void khash_destroy(khash_t(x) *map) noexcept {kh_destroy(x, map);}

_KHD(all)
_KHD(c)
_KHD(64)
_KHD(p)

template<> void khash_destroy(khash_t(name) *map) noexcept {destroy_name_hash(map);}

static ClassLevel get_linelvl(const char *line, std::string &buffer, const std::unordered_map<std::string, ClassLevel> &map) {
    const char *p(strchr(line, '|'));
    if(!p || (p = strchr(p + 1, '|')) == nullptr)
        RUNTIME_ERROR("Improperly formatted line");
    p = p + 2;
    const char *q(p);
    while(*q != '\t' && *q) ++q;
    if(!*q) RUNTIME_ERROR("Improperly formatted line");
    buffer = std::string(p, q);
    auto m(map.find(buffer));
    if(m == map.end()) {
        for(const auto &pair: map) {
            std::cerr << "Key: " << pair.first << ". Value: " << static_cast<int>(pair.second) << ".\n";
        }
        RUNTIME_ERROR(std::string("Unexpected field entry '") + buffer + "' for line " + line);
    }
    return m->second;
}

static std::unordered_map<tax_t, ClassLevel> get_tax_depths(const khash_t(p) *taxmap, const char *path) {
    std::unordered_map<tax_t, ClassLevel> ret;
    std::ifstream ifs(path);
    std::string buffer;
    tax_t t;
    if(!ifs.good()) RUNTIME_ERROR(std::string("could not open file at ") + path);
    for(std::string line; std::getline(ifs, line);) {
        t = atoi(line.data());
        if(kh_get(p, taxmap, t) == kh_end(taxmap)) continue;
        ret.emplace(t, get_linelvl(line.data(), buffer, classlvl_map));
    }
    return ret;
}

static std::vector<tax_t> get_sorted_taxes(const khash_t(p) *taxmap, const char *path) {
    std::vector<tax_t> taxes;
    {
        std::unordered_set<tax_t> taxset;
        for(khiter_t ki(0); ki != kh_end(taxmap); ++ki) if(kh_exist(taxmap, ki)) taxset.insert(kh_key(taxmap, ki));
        taxes = std::vector<tax_t>(taxset.begin(), taxset.end());
    }
    std::unordered_map<tax_t, ClassLevel> taxclassmap(get_tax_depths(taxmap, path));
    typename std::unordered_map<tax_t, ClassLevel>::iterator ma, mb;
    SORT(taxes.begin(), taxes.end(), [&tcm=taxclassmap,tm=taxmap,&ma,&mb](const tax_t a, const tax_t b) {
        return (ma = tcm.find(a)) == tcm.end() ? (mb = tcm.find(b)) == tcm.end() ? a < b: false
                                               : (mb = tcm.find(b)) == tcm.end() ? true
                                                                                 : ma->second == mb->second ? get_parent(tm, a) < get_parent(tm, b)
                                                                                                            : ma->second > mb->second;
    });
    return taxes;
}

static std::vector<tax_t> get_desc_lca(tax_t a, tax_t b, const std::unordered_map<tax_t, std::vector<tax_t>> &parent_map, const khash_t(p) *taxmap) {
    return get_all_descendents(parent_map, lca(taxmap, a, b));
}

static void print_name_hash(khash_t(name) *hash) noexcept {
    for(khiter_t ki(0); ki < kh_size(hash); ++ki)
        if(kh_exist(hash, ki))
            std::fprintf(stderr, "Key: %s. value: %u\n", kh_key(hash, ki), kh_val(hash, ki));
}

static tax_t get_max_val(const khash_t(p) *hash) noexcept {
    tax_t mx(std::numeric_limits<tax_t>::min());
    for(khiter_t ki(0); ki < kh_size(hash); ++ki)
        if(kh_exist(hash, ki))
            mx = std::max(std::max(kh_key(hash, ki), kh_val(hash, ki)), mx);
    return mx;
}

static khash_t(all) *load_binary_kmerset(const char *path) {
    std::FILE *fp(std::fopen(path, "rb"));
    if(fp == nullptr) throw std::system_error(std::error_code(2, std::system_category()), std::string("Cannot open path at ") + path + ".\n");
    khash_t(all) *ret(kh_init(all));
    u64 n;
    std::fread(&n, 1, sizeof(n), fp);
    if(kh_resize(all, ret, n) < 0) LOG_EXIT("Could not resize hash table to next power of 2 above %zu. New size: %zu\n", n, kh_n_buckets(ret));
    LOG_DEBUG("About to place %zu elements into a hash table of max size %zu\n", n, kh_n_buckets(ret));
    for(int khr; std::fread(&n, 1, sizeof(u64), fp) == sizeof(u64); kh_put(all, ret, n, &khr));
    std::fclose(fp);
#if defined(VERIFY_KHASH_LOAD)
#if VERIFY_KHASH_LOAD
    // Just make sure it all worked.
    for(khiter_t ki(0); ki < kh_end(ret); ++ki) {
        if(kh_exist(ret, ki)) assert(kh_get(all, ret, kh_key(ret, ki)) != kh_end(ret));
    }
    fp = std::fopen(path, "rb");
    std::fread(&n, 1, sizeof(u64), fp); // Skip first number.
    while(std::fread(&n, 1, sizeof(u64), fp) == sizeof(u64)) assert(kh_get(all, ret, n) != kh_end(ret));
    std::fclose(fp);
#endif // If the value is 1
#endif // If the value is defined
    return ret;
}


static lazy::vector<u64, size_t> load_binary_kmers(const char *path) {
    std::FILE *fp(std::fopen(path, "rb"));
    if(fp == nullptr) throw std::system_error(std::error_code(2, std::system_category()), std::string("Cannot open path at ") + path + ".\n");
    lazy::vector<u64, size_t> ret;
    u64 n;
    ::read(fileno(fp), static_cast<void *>(&n), sizeof(n));
    ret.resize(n, lazy::LAZY_VEC_NOINIT);
    ssize_t nread;
    if((nread = ::read(fileno(fp), static_cast<void *>(&ret[0]), n * sizeof(u64))) != ssize_t(n * sizeof(u64)))
        RUNTIME_ERROR(ks::sprintf("Only read %zd bytes from file when expected %zu.", nread, size_t(n * sizeof(u64))).data());
    std::fclose(fp);
    return ret;
}

static std::vector<std::string> get_paths(const char *path) {
    std::ifstream is(path);
    std::vector<std::string> ret;
    for(std::string line;std::getline(is, line);ret.emplace_back(std::move(line)));
    return ret;
}

static ssize_t filesize(const char* filename)
{
    struct stat st;
    stat(filename, &st);
    return st.st_size;
}
static int filesize(const int fd)
{
    struct stat st;
    fstat(fd, &st);
    return st.st_size;
}
static int filesize(std::FILE *fp)
{
    return filesize(fileno(fp));
}


enum WRITE {
    UNCOMPRESSED = 0,
    ZLIB = 1,
    ZSTD = 2
};

template<bool Condition>
using disable_if_t = typename std::enable_if_t<!Condition>;

namespace detail {

struct HasValueTypeImpl
{
    template<class T>
    static auto test(T&&) -> decltype( std::declval<typename T::value_type>(), std::true_type() );
    static auto test(...) -> std::false_type;
};

template<class T>
using HasValueType = decltype( HasValueTypeImpl::test(std::declval<T>()) );

template<class T, class U>
using IsConstructible = typename std::is_constructible<T, U>::type;

template<class Container>
class back_emplace_iterator
    : public std::iterator<std::output_iterator_tag, void, void, void, void>
{
    template<class T>
    using IsSelf = typename std::is_same<std::decay_t<T>, back_emplace_iterator>::type;

    Container &container;

public:
    typedef Container container_type;

    explicit back_emplace_iterator(Container& x): container(x){}
    template<class T, typename=disable_if_t<IsSelf<T>::value>>
    back_emplace_iterator& operator =(T&& t)
    {
        container.emplace_back(std::forward<T>(t));
        return *this;
    }
    template<class T=typename Container::value_type, class=std::enable_if_t<HasValueType<T>::value && IsConstructible<T, std::initializer_list<typename T::value_type>>::value>>
    back_emplace_iterator& operator =(std::initializer_list<typename T::value_type> ilist)
    {
        container.emplace_back(ilist);
        return *this;
    }
    back_emplace_iterator& operator =(typename Container::value_type&& t)
    {
        container.emplace_back(std::move(t));
        return *this;
    }
    back_emplace_iterator& operator *() { return *this; }
    back_emplace_iterator& operator ++() { return *this; }
    back_emplace_iterator& operator ++(int) { return *this; }
};

template<class Container>
inline back_emplace_iterator<Container> back_emplacer(Container& c) {
    return back_emplace_iterator<Container>(c);
}

} // namespace detail



} // namespace bns
namespace std {
    using bns::detail::back_emplacer;
    using bns::detail::back_emplace_iterator;
} // na
