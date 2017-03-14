#ifndef _KS_WRAPPER_H__
#define _KS_WRAPPER_H__

#include <cstdint>
#include <cstdio>
#include "klib/kstring.h"

namespace ks {

class KString {
    kstring_t ks_;

public:

    explicit KString(std::size_t size): ks_({0, size, size ? (char *)std::malloc(size): nullptr}) {}
    explicit KString(std::size_t used, std::size_t max, char *str): ks_({used, max, str}) {}
    explicit KString(char *str): ks_({0, 0, str}) {}
    KString(): KString(nullptr) {}
    ~KString() {free(ks_.s);}

    const auto operator->() const {
        return const_cast<const kstring_t *>(&ks_);
    }
    kstring_t *operator->() {return &ks_;}

    const auto &operator*() const {
        return const_cast<const kstring_t &>(ks_);
    }
    kstring_t  &operator*() {return ks_;}

    // Copy
    KString(const KString &other): ks_{other->l, other->m, (char *)std::malloc(other->m)} {
        memcpy(ks_.s, other->s, other->m);
    }

    // Move
    KString(KString &&other) {
        memcpy(this, &other, sizeof(other));
        memset(&other, 0, sizeof(other));
    }
    int cmp(const char *s) {
        for(char *s2(ks_.s);*s && *s2;++s, ++s2) if(*s != *s2) return *s2 - *s;
        return 0;
    }
    int cmp(const KString &other) {return cmp(other->s);}

    bool operator==(const KString &other) {
        if(other->l != ks_.l) return 0;
        for(std::size_t i(0); i < ks_.l; ++i) if(ks_.s[i] != other->s[i]) return 0;
        return 1;
    }

    int putc(int c) {return kputc(c, &ks_);}
    int putw(int c) {return kputw(c, &ks_);}
    int sprintf(const char *fmt, ...) {
        va_list ap;
        va_start(ap, fmt);
        const int ret(kvsprintf(&ks_, fmt, ap));
        va_end(ap);
        return ret;
    }
    operator char *()   const {return ks_.s;}
    operator kstring_t *() {
        return &ks_;
    }
    operator const kstring_t *() const {
        return static_cast<const kstring_t *>(&ks_);
    }

    const char     *data() const {return ks_.s;}
    char           *data() {return ks_.s;}

    char  *release() {auto ret(ks_.s); ks_.l = ks_.m = 0; ks_.s = nullptr; return ret;}

    std::size_t  size() const {return ks_.l;}
    auto        begin() const {return ks_.s;}
    auto          end() const {return ks_.s + ks_.l;}
    const auto cbegin() const {return const_cast<const char *>(ks_.s);}
    const auto   cend() const {return const_cast<const char *>(ks_.s + ks_.l);}
    void pop() {ks_.s[--ks_.l] = 0;}
    void pop(std::size_t n) {
        ks_.l = ks_.l > n ? ks_.l - n: 0;
        ks_.s[ks_.l] = 0;
    }

    void clear() {
        ks_.l = 0;
        ks_.s = 0;
    }

    auto resize(std::size_t new_size) {
        return ks_resize(&ks_, new_size);
    }
};

} // namespace ks

#endif // #ifndef _KS_WRAPPER_H__
