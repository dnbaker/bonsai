#ifndef _KS_WRAPPER_H__
#define _KS_WRAPPER_H__

#ifndef _GNU_SOURCE
#  define _GNU_SOURCE 1
#endif

#if ZWRAP_USE_ZSTD
#  include "zstd_zlibwrapper.h"
#else
#  include <zlib.h>
#endif

#include <cstdint>
#include <cassert>
#include <cinttypes>
#include <cstdarg>
#include <cstdio>
#include <cstring>
#include <iostream>
#include <iostream>
#include <stdexcept>
#include <string>
#include <vector>
#include <unistd.h>
//#include <experimental/functional>
// If this fails to be located and you have a new compiler, you may need to remove "experimental/" from this include.
#include <algorithm>


#ifndef roundup64__
#define roundup64__(x) (--(x), (x)|=(x)>>1, (x)|=(x)>>2, (x)|=(x)>>4, (x)|=(x)>>8, (x)|=(x)>>16, (x)|=(x)>>32, ++(x))
#endif

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
#  ifndef INLINE
#    define INLINE __attribute__((always_inline)) inline
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
#  ifndef INLINE
#    define INLINE
#  endif
#endif

#ifndef ENABLE_BM_SEARCH
#undef __cpp_lib_boyer_moore_searcher
#endif


namespace ks {

#define realloc_unique(ptr, new_size) do {\
    std::decay_t<decltype(ptr.get())> __tmp_ptr(std::realloc(ptr.get(), new_size));\
    if(__tmp_ptr == nullptr) throw std::bad_alloc("Could not reallocate unique_ptr\n");\
    ptr.release(); ptr.reset(__tmp_ptr);\
    } while(0)


using std::uint64_t;
using namespace std::literals;

static std::vector<int> ksBM_prep(const ::std::uint8_t *pat, int m)
{
    int i, *bmGs, *bmBc;
    std::vector<int> prep(m + 256);
    bmGs = prep.data(); bmBc = prep.data() + m;
    { // preBmBc()
        for (i = 0; i < 256; ++i) bmBc[i] = m;
        for (i = 0; i < m - 1; ++i) bmBc[pat[i]] = m - i - 1;
    }
    std::vector<int> suff(m);
    { // suffixes()
        int f = 0, g;
        suff[m - 1] = m;
        g = m - 1;
        for (i = m - 2; i >= 0; --i) {
            if (i > g && suff[i + m - 1 - f] < i - g)
                suff[i] = suff[i + m - 1 - f];
            else {
                if (i < g) g = i;
                f = i;
                while (g >= 0 && pat[g] == pat[g + m - 1 - f]) --g;
                suff[i] = f - g;
            }
        }
    }
    { // preBmGs()
        int j = 0;
        for (i = 0; i < m; ++i) bmGs[i] = m;
        for (i = m - 1; i >= 0; --i)
            if (suff[i] == i + 1)
                for (; j < m - 1 - i; ++j)
                    if (bmGs[j] == m)
                        bmGs[j] = m - 1 - i;
        for (i = 0; i <= m - 2; ++i)
            bmGs[m - 1 - suff[i]] = m - 1 - i;
    }
    return prep;
}

static inline auto ksBM_prep(const std::string &str) {
    return ksBM_prep((const ::std::uint8_t *)str.data(), str.size());
}

static inline void *kmemmem(const void *_str, int n, const void *_pat, int m, const std::vector<int> &prep)
{
    int i, j;
    const int *bmGs, *bmBc;
    const ::std::uint8_t *str, *pat;
    str = (const ::std::uint8_t*)_str; pat = (const ::std::uint8_t*)_pat;
    bmGs = prep.data(); bmBc = prep.data() + m;
    j = 0;
    while (j <= n - m) {
        for (i = m - 1; i >= 0 && pat[i] == str[i+j]; --i);
        if (i >= 0) {
            int max = bmBc[str[i+j]] - m + 1 + i;
            if (max < bmGs[i]) max = bmGs[i];
            j += max;
        } else return (void*)(str + j);
    }
    return 0;
}



class string {
    uint64_t l, m;
    char     *s;
public:

    static const uint64_t DEFAULT_SIZE = 4;
    using value_type = char;
    using size_type  = uint64_t;
/*
TODO: Add SSO to avoid allocating for small strings, which we currently do
      defensively in order to avoid segfaults.
*/
    void default_allocate() {
        if(likely(s != nullptr)) return;
        if((s = static_cast<char *>(std::malloc(DEFAULT_SIZE))) == nullptr) throw std::bad_alloc();
        if(m < DEFAULT_SIZE) {
            m = DEFAULT_SIZE;
            l = 0;
            *s = 0;
        }
    }

    INLINE explicit string(uint64_t size): l(0), m(size), s(m ? static_cast<char *>(std::malloc(m * sizeof(char))): nullptr) {
        default_allocate();
    }
    operator const char *() const {return s;}

    inline string(uint64_t used, uint64_t max, const char *str):
        l(used), m(max)  {
        if((s = static_cast<char *>(std::malloc(m * sizeof(char)))) == nullptr) throw std::bad_alloc();
        std::memcpy(s, str, (l + 1) * sizeof(char));
        default_allocate();
    }
    inline string(char *str, size_t len): l(len), m(len), s(str) { // Stealing the other thing.
#if !NDEBUG
        std::fprintf(stderr, "[%s:%s:%d] Acquired ownership of string at %p with len %zu has been taken.", __PRETTY_FUNCTION__, __FILE__, __LINE__, static_cast<const void *>(str), len);
#endif
        terminate();
    }
    inline string(const char *str, uint64_t used): string(used, used, str) {}

    INLINE explicit string(const char *str) {
        if(str == nullptr) {
            m = l = 0;
            s = nullptr;
            default_allocate();
        } else {
            m = 0;
            s = nullptr;
            l = std::strlen(str);
            resize(l + 1);
            assert(s);
            std::memcpy(s, str, (l + 1) * sizeof(char));
        }
    }

    INLINE string(): l{0}, m{0}, s{0} {
        default_allocate();
    }
    INLINE ~string() {std::free(s);}

#ifdef KSTRING_H
    // Access kstring
    INLINE kstring_t *ks()             {return reinterpret_cast<kstring_t *>(this);}
    INLINE const kstring_t *ks() const {return reinterpret_cast<const kstring_t *>(this);}
#endif
    INLINE void free() {
        std::free(s);
        std::memset(this, 0, sizeof(*this));
    }
    template<typename T>
    string &append(const T &to_append) {
        *this += to_append;
        return *this;
    }
    string &append(const char *str, size_t len) {
        this->resize(this->l + len);
        std::memcpy(this->s + this->l, str, len);
        this->l += len;
        terminate();
        return *this;
    }
    string &append(size_t n, char c) {
        this->resize(this->l + n);
        for(size_t final_len = this->l + n; this->l != final_len; this->s[this->l++] = c);
        terminate();
        return *this;
    }

    // Copy
    INLINE string(const string &other): l(other.l), m(other.m), s(static_cast<char *>(std::malloc(other.m))) {
        std::memcpy(s, other.s, l + 1);
    }

    INLINE string(const std::string &str): l(str.size()), m(l), s(static_cast<char *>(std::malloc(m))) {
        roundup64__(m);
        std::memcpy(s, str.data(), (l + 1) * sizeof(char));
    }

    INLINE string &operator=(const string &other)    {return *this = string(other);}
    INLINE string &operator=(const char *str)        {return *this = string(str);}
    INLINE string &operator=(const std::string &str) {return *this = string(str);}
    INLINE string &operator=(string &&other)    {
        this->free();
        std::memcpy(this, &other, sizeof(*this));
        other.zero();
        return *this;
    }

    // Move
    INLINE string(string &&other) {
        std::memcpy(this, &other, sizeof(other));
        other.zero();
    }
    INLINE auto       &len()       {return l;}
    INLINE const auto &len() const {return l;}

    // Comparison functions
    INLINE int cmp(const char *str)     const {return std::strcmp(s, str);}
    INLINE int cmp(const string &other) const {return cmp(other.s);}

    INLINE bool operator==(const string &other) const {
        return l == other.l && std::memcmp(this->s, other.s, l) == 0;
    }
    INLINE bool operator==(const ::std::string &other) const {
        return l == other.size() && std::memcmp(this->s, other.data(), l) == 0;
    }

    void zero() {
        std::memset(this, 0, sizeof(*this));
    }

    INLINE bool operator==(const char *str) const {
        return s ? str ? std::strcmp(str, s) == 0: 0: 1;
    }

    template<typename T> INLINE bool operator!=(const T &o) {return !(this->operator==(o));}

    INLINE bool palindrome() const {
        return std::equal(s, s + (l + 1) / 2, std::reverse_iterator<const char *>(s + l));
    }
    INLINE void reverse() {
        auto revit = std::reverse_iterator<char *>(s + l);
        for(auto fit = s, efit = fit + (l + 1) / 2; fit < efit; std::swap(*fit++, *revit++));
    }
    string reversed() const {
        string cpy(*this);
        cpy.reverse();
        return cpy;
    }
    bool startswith(const char *str, uint64_t slen) const {return std::memcmp(s, str, slen) == 0;}
    bool startswith(const char *str) const {return startswith(str, std::strlen(str));}
    template<typename T> bool startswith(const T &str) const {return startswith(str.data(), str.size());}
    bool endswith(const char *str, uint64_t slen) const {
        return std::memcmp(str, s + l - slen, slen) == 0;
    }
    bool endswith(const char *str) const {return endswith(str, std::strlen(str));}
    template<typename T> bool endswith(const T &str) const {return endswith(str.data(), str.size());}

    // Appending:
    INLINE int putc_(int c) {
        if (unlikely(l + 1 >= m)) {
            char *tmp;
            m = l + 2;
            roundup64__(m);
            if ((tmp = static_cast<char *>(std::realloc(s, m * sizeof(char)))))
                s = tmp;
            else
                return EOF;
        }
        s[l++] = (char)c;
        return 0;
    }
    INLINE int putw_(int c)  {
        char buf[16];
        int i, len = 0;
        unsigned int x = c;
        if (c < 0) x = -x;
        do { buf[len++] = (char)(x%10 + '0'); x /= 10; } while (x > 0);
        if (c < 0) buf[len++] = '-';
        if (unlikely(len + l + 1 >= m)) {
            char *tmp;
            m = len + l + 2;
            roundup64__(m);
            if ((tmp = static_cast<char*>(std::realloc(s, m * sizeof(char)))))
                s = tmp;
            else
                return EOF;
        }
        for (i = len - 1; i >= 0; --i) s[l++] = buf[i];
        return 0;
    }
    INLINE int putuw_(int c) {
        char buf[16];
        int len, i;
        unsigned x;
        if (c == 0) return putc('0');
        for (len = 0, x = c; x > 0; x /= 10) buf[len++] = (char)(x%10 + '0');
        if (unlikely(len + l + 1 >= m)) {
            char *tmp;
            m = len + l + 2;
            roundup64__(m);
            if ((tmp = static_cast<char *>(std::realloc(s, m * sizeof(char)))))
                s = tmp;
            else
                return EOF;
        }
        for (i = len - 1; i >= 0; --i) s[l++] = buf[i];
        return 0;
    }
    INLINE int putl_(long c)  {
        char buf[32];
        int i, len = 0;
        unsigned long x = c;
        if (c < 0) x = -x;
        do { buf[len++] = (char)(x%10 + '0'); x /= 10; } while (x > 0);
        if (c < 0) buf[len++] = '-';
        if (unlikely(len + l + 1 >= m)) {
            char *tmp;
            m = len + l + 2;
            roundup64__(m);
            if ((tmp = static_cast<char *>(std::realloc(s, m * sizeof(char)))))
                s = tmp;
            else
                return EOF;
        }
        for (i = len - 1; i >= 0; --i) s[l++] = buf[i];
        return 0;
    }
    INLINE long putsn_(const char *str, long len) {
        if (unlikely(len + l + 1 >= m)) {
            char *tmp;
            m = len + l + 2;
            roundup64__(m);
            if ((tmp = static_cast<char *>(std::realloc(s, m * sizeof(char)))))
                s = tmp;
            else
                return EOF;
        }
        std::memcpy(s + l, str, len * sizeof(char));
        l += len;
        return len;
    }
    std::string to_std_string() const {
        return std::string(s, s + l);
    }
    INLINE char       &back()       {return s[l - 1];}
    INLINE const char &back() const {return s[l - 1];}

    INLINE char       &terminus()       {return s[l];}
    INLINE const char &terminus() const {return s[l];}
    INLINE void       terminate()       {terminus() = '\0';}
    INLINE int putuw(int c) {
        c = putuw_(c); s[l] = 0; return c;
    }
    INLINE int putc(int c) {
        c = putc_(c); s[l] = 0; return c;
    }
    INLINE int putw(int c)  {
        c = putw_(c); s[l] = 0; return c;
    }
    INLINE long putl(long c)  {
        c = putl_(c), s[l] = 0; return c;
    }
    INLINE long putsn(const char *str, long len)  {
        len = putsn_(str, len); s[l] = 0; return l;
    }
    INLINE long int puts(const char *s) {return putsn_(s, std::strlen(s));}
    int vsprintf(const char *fmt, va_list ap)
    {
        va_list args;
        va_copy(args, ap);
        int len(vsnprintf(s + l, m - l, fmt, args)); // This line does not work with glibc 2.0. See `man snprintf'.
        va_end(args);
        if (unlikely((unsigned)len + 1 > m - l)) {
            resize(l + len + 2);
            va_copy(args, ap);
            len = vsnprintf(s + l, m - l, fmt, args);
            va_end(args);
        }
        l += len;
        return len;
    }

    int sprintf(const char *fmt, ...) {
        va_list ap;
        va_start(ap, fmt);
        const int l(this->vsprintf(fmt, ap));
        va_end(ap);
        return l;
    }

    // Transfer ownership
    char  *release() {auto ret(s); l = m = 0; s = nullptr; return ret;}

    // STL imitation
    INLINE uint64_t size()     const {return l;}
    INLINE uint64_t capacity() const {return m;}
    INLINE auto  begin()     const {return s;}
    INLINE auto    end()     const {return s + l;}
    INLINE auto cbegin()     const {return const_cast<const char *>(s);}
    INLINE auto   cend()     const {return const_cast<const char *>(s + l);}
    INLINE char pop() {const char ret(s[--l]); s[l] = 0; return ret;}
    INLINE void pop(uint64_t n) {
        l = l > n ? l - n: 0;
        s[l] = 0;
    }

    void clear() {l = 0; if(s) *s = '\0';}

    INLINE const char     *data() const {return s;}
    INLINE char           *data()       {return s;}
    // In-place modify std::string.
    std::string &set(std::string &ret) const {
        ret.reserve(l);
        ret.resize(0);
        std::copy(begin(), end(), std::back_inserter(ret));
        return ret;
    }
    std::string str() const {
        std::string ret;
        set(ret);
        return ret;
    }
    INLINE void set_size(uint64_t newlen) {
        l = newlen;
    }
    INLINE int resize(uint64_t size) {
        if (m < size) {
            char *tmp;
            m = std::max(size, UINT64_C(4));
            roundup64__(m);
            if ((tmp = static_cast<char*>(std::realloc(s, m * sizeof(char)))) == nullptr) {
                std::cerr << ("Could not allocate sufficient memory for "s + std::to_string(m) + " bytes.\n");
                throw std::bad_alloc();
            }
            s = tmp;
        }
        return 0;
    }

    // std::string imitation: Appending
    // Append char
    INLINE auto &operator+=(const char c) {putc(c);  return *this;}

    // Append formatted int/unsigned/long
    INLINE auto &operator+=(int c)        {putw(c);  return *this;}
    INLINE auto &operator+=(unsigned c)   {putuw(c); return *this;}
    INLINE auto &operator+=(long c)       {putl(c);  return *this;}
    char *locate(const char *str, uint64_t len) {
        return (char *)memmem(s, l, str, len);
    }
    const char *locate(const char *str, uint64_t len) const {return static_cast<const char *>(const_cast<string *>(this)->locate(str, len));}
    const char *locate(const char *str) const {return locate(str, std::strlen(str));}
    char *locate(const char *str) {return locate(str, std::strlen(str));}

    char *bmlocate(const char *str, uint64_t len) {
#if __cpp_lib_boyer_moore_searcher
        std::boyer_moore_searcher searcher(str, str + len);
        return std::search<const char *, decltype(searcher)>(s, s + l, searcher);
#else
        auto prep = ksBM_prep((const ::std::uint8_t *)str, len);
        return (char *)kmemmem((const ::std::uint8_t *)s, l, (const ::std::uint8_t *)str, len, prep);
#endif
    }
    char *bmhlocate(const char *str, uint64_t len) {
#if __cpp_lib_boyer_moore_searcher
        std::boyer_moore_horspool_searcher searcher(str, str + len);
        return std::search(s, s + l, searcher);
#else
        auto prep = ksBM_prep((const ::std::uint8_t *)str, len);
        return (char *)kmemmem(s, l, str, len, prep);
#endif
    }
#if __cpp_lib_boyer_moore_searcher
    auto make_bm() const {
        return std::boyer_moore_searcher(s, s + l);
    }
    auto make_bmh() const {
        return std::boyer_moore_horspool_searcher(s, s + l);
    }
#endif
    const char *bmlocate(const char *str, uint64_t len) const {return static_cast<const char *>(const_cast<string *>(this)->bmlocate(str, len));}
    const char *bmlocate(const char *str) const {return bmlocate(str, std::strlen(str));}
    char *bmlocate(const char *str) {return bmlocate(str, std::strlen(str));}
    const char *bmhlocate(const char *str, uint64_t len) const {return static_cast<const char *>(const_cast<string *>(this)->bmhlocate(str, len));}
    const char *bmhlocate(const char *str) const {return bmhlocate(str, std::strlen(str));}
    char *bmhlocate(const char *str) {return bmhlocate(str, std::strlen(str));}
    bool contains(const char *str, uint64_t len) const {return locate(str, len) != nullptr;}
    bool contains(const char *str) const {return contains(str, std::strlen(str));}
    template<typename T> bool contains(const T &str) const {return contains(str.data(), str.size());}
    bool bmcontains(const char *str, uint64_t len) const {
        return bmlocate(str, len) != end();
    }
    bool bmcontains(const char *str) const {return bmcontains(str, std::strlen(str));}
    template<typename T> bool bmcontains(const T &str) const {return bmcontains(str.data(), str.size());}

    // Append string forms
#ifdef KSTRING_H
    INLINE auto &operator+=(const kstring_t *ks) {
        putsn(ks->s, ks->l);
        return *this;
    }
    INLINE auto &operator+=(const kstring_t &ks) {
        return operator+=(&ks);
    }
#endif
    INLINE auto &operator+=(const std::string &s) {
        putsn(s.data(), s.size());
        return *this;
    }
    INLINE auto &operator+=(const string &other) {putsn(other.s, other.l); return *this;}
    INLINE auto &operator+=(const char *s)       {puts(s); return *this;}

    // Access
    INLINE const char &operator[](uint64_t index) const {return s[index];}
    INLINE char       &operator[](uint64_t index)       {return s[index];}

    INLINE size_t write(FILE *fp) const   {return std::fwrite(s, sizeof(char), l, fp);}
    INLINE auto write(const char *path) const {
        std::FILE *fp(std::fopen(path, "r"));
        if(!fp) throw 1;
        const auto ret(write(fp));
        std::fclose(fp);
        return ret;
    }
    INLINE ssize_t write(int fd) const {return    ::write(fd, s, l * sizeof(char));}
    INLINE auto write(gzFile fp) const {return    gzwrite(fp, s, l * sizeof(char));}
    template<typename T> auto flush(T target) {
        auto ret = this->write(target);
        this->clear();
        return ret;
    }
};

// s MUST BE a null terminated string; [l = strlen(s)]
template<typename T=std::uint64_t, typename Alloc=std::allocator<T>, typename=typename std::enable_if<std::is_arithmetic<T>::value>::type>
void split(char *s, int delimiter, uint64_t l, std::vector<T, Alloc> &offsets)
{
    unsigned i, last_char, last_start;
    offsets.clear();

#define _split_aux_ do {s[i] = 0, offsets.push_back(last_start);} while(0)

    for (i = 0, last_char = last_start = 0; i <= l; ++i) {
        if (delimiter == 0) {
            if (std::isspace(s[i]) || s[i] == 0) {
                if (std::isgraph(last_char))                  _split_aux_; // the end of a field
            } else {
                if (std::isspace(last_char) || last_char == 0) last_start = i;
            }
        } else {
            if (s[i] == delimiter || s[i] == 0) {
                if (last_char != 0 && last_char != static_cast<unsigned>(delimiter)) _split_aux_; // the end of a field
            } else {
                if (last_char == static_cast<unsigned>(delimiter) || last_char == 0) last_start = i;
            }
        }
        last_char = s[i];
    }

#undef _split_aux_

}

template<typename T=std::uint64_t, typename=typename std::enable_if<std::is_arithmetic<T>::value>::type>
inline void split(char *s, int delimiter, std::vector<T> &offsets) {
    split(s, delimiter, std::strlen(s), offsets);
}


template<typename T=std::uint64_t, typename Alloc=std::allocator<T>, typename=typename std::enable_if<std::is_arithmetic<T>::value>::type>
inline std::vector<T, Alloc> split(char *s, uint64_t l, int delimiter=0)
{
    std::vector<T, Alloc> ret;
    ks::split(s, delimiter, l, ret);
    return ret;
}

static string sprintf(const char *fmt, ...) {
    string ret;
    va_list ap;
    va_start(ap, fmt);
    ret.vsprintf(fmt, ap);
    va_end(ap);
    return ret;
}

template<typename T=std::uint64_t, typename Alloc=std::allocator<T>, typename=typename std::enable_if<std::is_arithmetic<T>::value>::type>
inline std::vector<string> toksplit(char *s, uint64_t l, int delimiter=0) {
    auto vec(split<T, Alloc>(s, l, delimiter));
    std::vector<string> ret;
    ret.reserve(vec.size());
    for(const auto i: vec) ret.emplace_back(s + i);
    return ret;
}

template<typename T=std::uint64_t, typename Alloc=std::allocator<T>, typename=typename std::enable_if<std::is_arithmetic<T>::value>::type>
std::vector<T, Alloc> split(string &s, int delimiter=0) {return split<T, Alloc>(s.data(), s.size(), delimiter);}
template<typename T=std::uint64_t, typename Alloc=std::allocator<T>, typename=typename std::enable_if<std::is_arithmetic<T>::value>::type>
std::vector<T, Alloc> split(std::string &s, int delimiter=0) {return split<T, Alloc>(&s[0], s.size(), delimiter);}
template<typename T=std::uint64_t, typename Alloc=std::allocator<T>, typename=typename std::enable_if<std::is_arithmetic<T>::value>::type>
std::vector<T, Alloc> split(char *s, int delimiter=0) {return split<T, Alloc>(s, std::strlen(s), delimiter);}

using KString = ks::string;

} // namespace ks

using ks::KString;

#undef roundup64__

#endif // #ifndef _KS_WRAPPER_H__
