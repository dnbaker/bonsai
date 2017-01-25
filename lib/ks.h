#ifndef _KS_WRAPPER_H__
#define _KS_WRAPPER_H__

#include <cstdint>
#include <cstdio>

#include "klib/kstring.h"

namespace ks {

class KString {
    kstring_t ks;

public:
    explicit KString(std::size_t size): ks({0, size, size ? (char *)malloc(size): nullptr}) {
    }
    KString(): KString(0) {}
    ~KString() {free (ks.s);}

    int putc(int c) {return kputc(c, &ks);}
    int putw(int c) {return kputw(c, &ks);}
    int sprintf(const char *fmt, ...) {
        va_list ap;
        va_start(ap, fmt);
        const int ret(kvsprintf(&ks, fmt, ap));
        va_end(ap);
        return ret;
    }
    operator char *()      const {return ks.s;}
    operator kstring_t *() const {return &ks;}
};

} // namespace ks

#endif // #ifndef _KS_WRAPPER_H__
