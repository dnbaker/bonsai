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
    KString(): KString(0u, 0u, nullptr) {}

    ~KString() {free (ks_.s);}

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

    const auto operator->() const {
        return const_cast<const kstring_t *>(&ks_);
    }
    kstring_t *operator->() {return &ks_;}

    void clear() {
        ks_.l = 0;
    }

    auto resize(std::size_t new_size) {
        return ks_resize(&ks_, new_size);
    }
};

} // namespace ks

#endif // #ifndef _KS_WRAPPER_H__
