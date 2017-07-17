#ifndef _KS_WRAPPER_H__
#define _KS_WRAPPER_H__
#include <cstdint>
#include <cstring>
#include <cstdarg>
#include <cstdio>
#include <iostream>
#include <unistd.h>

#ifndef roundup64
#define roundup64(x) (--(x), (x)|=(x)>>1, (x)|=(x)>>2, (x)|=(x)>>4, (x)|=(x)>>8, (x)|=(x)>>16, (x)|=(x)>>32, ++(x))
#endif

namespace ks {

using std::size_t;

class KString {
    size_t l, m;
    char     *s;
public:

    explicit KString(size_t size): l(size), m(roundup64(size)), s(size ? static_cast<char *>(std::malloc(size * sizeof(char))): nullptr) {}

    explicit KString(size_t used, size_t max, char *str, bool assume_ownership=false):
        l(used), m(max), s(str) {
        if(assume_ownership == false) {
            s = static_cast<char *>(std::malloc(m * sizeof(char)));
            std::memcpy(s, str, (l + 1) * sizeof(char));
        }
    }

    explicit KString(const char *str) {
        if(str == nullptr) {
            std::memset(this, 0, sizeof *this);
        } else {
            l = std::strlen(str);
            m = roundup64(l);
            s = static_cast<char *>(std::malloc(m * sizeof(char)));
            std::memcpy(s, str, (l + 1) * sizeof(char));
        }
    }

    KString(): KString(nullptr) {}
    ~KString() {std::free(s);}

#ifdef KSTRING_H
    // Access kstring
    kstring_t *ks()             {return reinterpret_cast<kstring_t *>(this);}
    const kstring_t *ks() const {return reinterpret_cast<const kstring_t *>(this);}
#endif

    // Copy
    KString(const KString &other): l(other.l), m(other.m), s(static_cast<char *>(std::malloc(other.m))) {
        std::memcpy(s, other.s, l + 1);
    }

    KString(const std::string &str): l(str.size()), m(l), s(static_cast<char *>(std::malloc(m))) {
        roundup64(m);
        std::memcpy(s, str.data(), (l + 1) * sizeof(char));
    }

    // Stealing ownership in a very mean way.
    KString(std::string &&str): l(str.size()), m(l), s(const_cast<char *>(str.data())) {
        roundup64(m);
        std::memset(&str, 0, sizeof(str));
    }

    KString operator=(const KString &other)   {return KString(other);}
    KString operator=(const char *str)        {return KString(str);}
    KString operator=(const std::string &str) {return KString(str);}

    // Move
    KString(KString &&other) {
        std::memcpy(this, &other, sizeof(other));
        std::memset(&other, 0, sizeof(other));
    }

    // Comparison functions
    int cmp(const char *str)      const {return std::strcmp(s, str);}
    int cmp(const KString &other) const {return cmp(other.s);}

    bool operator==(const KString &other) const {
        if(other.l != l) return 0;
        if(l) for(size_t i(0); i < l; ++i) if(s[i] != other.s[i]) return 0;
        return 1;
    }

    bool operator==(const char *str) const {
        return s ? str ? std::strcmp(str, s) == 0: 0: 1;
    }

    bool operator==(const std::string &str) const {
        return str.size() == l ? l ? std::strcmp(str.data(), s) == 0: 1: 0;
    }

    bool palindrome() const {
        for(size_t i(0), e(l >> 1); i < e; ++i)
            if(s[i] != s[l - i - 1])
                return 0;
        return 1;
    }

    // Appending:
    int putc_(int c) {
        if (l + 1 >= m) {
            char *tmp;
            m = l + 2;
            roundup64(m);
            if ((tmp = static_cast<char *>(std::realloc(s, m * sizeof(char)))))
                s = tmp;
            else
                return EOF;
        }
        s[l++] = c;
        return 0;
    }
    int putw_(int c)  {
        char buf[16];
        int i, len = 0;
        unsigned int x = c;
        if (c < 0) x = -x;
        do { buf[len++] = x%10 + '0'; x /= 10; } while (x > 0);
        if (c < 0) buf[len++] = '-';
        if (len + l + 1 >= m) {
            char *tmp;
            m = len + l + 2;
            roundup64(m);
            if ((tmp = static_cast<char*>(std::realloc(s, m * sizeof(char)))))
                s = tmp;
            else
                return EOF;
        }
        for (i = len - 1; i >= 0; --i) s[l++] = buf[i];
        return 0;
    }
    int putuw_(int c) {
        char buf[16];
        int len, i;
        unsigned x;
        if (c == 0) return putc('0');
        for (len = 0, x = c; x > 0; x /= 10) buf[len++] = x%10 + '0';
        if (len + l + 1 >= m) {
            char *tmp;
            m = len + l + 2;
            roundup64(m);
            if ((tmp = static_cast<char *>(std::realloc(s, m * sizeof(char)))))
                s = tmp;
            else
                return EOF;
        }
        for (i = len - 1; i >= 0; --i) s[l++] = buf[i];
        return 0;
    }
    int putl_(long c)  {
        char buf[32];
        int i, len = 0;
        unsigned long x = c;
        if (c < 0) x = -x;
        do { buf[len++] = x%10 + '0'; x /= 10; } while (x > 0);
        if (c < 0) buf[len++] = '-';
        if (len + l + 1 >= m) {
            char *tmp;
            m = len + l + 2;
            roundup64(m);
            if ((tmp = static_cast<char *>(std::realloc(s, m * sizeof(char)))))
                s = tmp;
            else
                return EOF;
        }
        for (i = len - 1; i >= 0; --i) s[l++] = buf[i];
        return 0;
    }
    int putuw(int c) {
        c = putuw_(c);
        s[l] = 0;
        return c;
    }
    long putsn_(const char *str, long len) {
        if (len + l + 1 >= m) {
            char *tmp;
            m = len + l + 2;
            roundup64(m);
            if ((tmp = static_cast<char *>(std::realloc(s, m * sizeof(char)))))
                s = tmp;
            else
                return EOF;
        }
        std::memcpy(s + l, str, len * sizeof(char));
        l += len;
        return len;
    }
    int putc(int c) {
        c = putc_(c), s[l] = 0;
        return c;
    }
    char       &back()       {return s[l - 1];}
    const char &back() const {return s[l - 1];}
    int putw(int c)  {
        c = putw_(c), s[l] = 0;
        return c;
    }
    int putl(long c)  {
        c = putl_(c), s[l] = 0;
        return c;
    }
    template<typename FMT=char, typename=std::enable_if_t<std::is_same<char, FMT>::value>>
    int puts(const char *s) {return putsn_(s, std::strlen(s) + 1);}
    template<typename FMT=char, typename=std::enable_if_t<std::is_same<char, FMT>::value>>
    long putsn(const char *str, long len)  {
        len = putsn_(str, len);
        s[l] = 0;
        return l;
    }
    int sprintf(const char *fmt, ...) {
        size_t len;
        std::va_list ap;
        va_start(ap, fmt);
        len = vsnprintf(s + l, m - l, fmt, ap); // This line does not work with glibc 2.0. See `man snprintf'.
        if (len + 1 > m - l) {
            resize(len + 1);
            len = vsnprintf(s + l, m - l, fmt, ap);
        }
        va_end(ap);
        l += len;
        return len;
    }

    // Transfer ownership
    char  *release() {auto ret(s); l = m = 0; s = nullptr; return ret;}

    // STL imitation
    size_t size() const {return l;}
    auto  begin() const {return s;}
    auto    end() const {return s + l;}
    auto cbegin() const {return const_cast<const char *>(s);}
    auto   cend() const {return const_cast<const char *>(s + l);}
    char pop() {const char ret(s[--l]); s[l] = 0; return ret;}
    void pop(size_t n) {
        l = l > n ? l - n: 0;
        s[l] = 0;
    }

    void clear() {l = 0; s[0] = '\0';}

    const char     *data() const {return s;}
    char           *data()       {return s;}
    int resize(size_t size) {
        if (m < size) {
            char *tmp;
            m = size;
            roundup64(m);
            if ((tmp = static_cast<char*>(std::realloc(s, m * sizeof(char)))))
                s = tmp;
            else
                return -1;
        }
        return 0;
    }

    // std::string imitation: Appending
    // Append char
    auto &operator+=(const char c) {putc(c);  return *this;}

    // Append formatted int/unsigned/long
    auto &operator+=(int c)        {putw(c);  return *this;}
    auto &operator+=(unsigned c)   {putuw(c); return *this;}
    auto &operator+=(long c)       {putl(c);  return *this;}

    // Append string forms
#ifdef KSTRING_H
    template<typename FMT=char, typename=std::enable_if_t<std::is_same<char, FMT>::value>>
    auto &operator+=(const kstring_t *ks) {
        putsn(ks->s, ks->l);
        return *this;
    }
    template<typename FMT=char, typename=std::enable_if_t<std::is_same<char, FMT>::value>>
    auto &operator+=(const kstring_t &ks) {
        return operator+=(&ks);
    }
#endif
    template<typename FMT=char, typename=std::enable_if_t<std::is_same<char, FMT>::value>>
    auto &operator+=(const std::string &s) {
        putsn(s.data(), s.size());
        return *this;
    }
    template<typename FMT=char, typename=std::enable_if_t<std::is_same<char, FMT>::value>>
    auto &operator+=(const KString &other) {putsn(other.s, other.l); return *this;}
    template<typename FMT=char, typename=std::enable_if_t<std::is_same<char, FMT>::value>>
    auto &operator+=(const char *s)        {puts(s); return *this;}

    // Access
    const char &operator[](size_t index) const {return s[index];}
    char       &operator[](size_t index)       {return s[index];}

    int write(FILE *fp) const {return std::fwrite(s, sizeof(char), l, fp);}
    int write(int fd)   const {return     ::write(fd, s, l * sizeof(char));}
};

} // namespace ks

#endif // #ifndef _KS_WRAPPER_H__
