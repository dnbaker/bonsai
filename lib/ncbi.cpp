#include "ncbi.h"
#include <cstdio>
#include <cstdlib>
#include <set>

namespace kpg {
khash_t(p) *build_parent_map(const char *fn) {
    size_t nlines(count_lines(fn));
    FILE *fp(fopen(fn, "r"));
    khash_t(p) *ret(kh_init(p));
    kh_resize(p, ret, nlines);
    size_t bufsz = 4096;
    char *buf((char *)malloc(bufsz));
    ssize_t len;
    khint_t ki;
    int khr;
    while((len = getline(&buf, &bufsz, fp)) >= 0) {
        switch(*buf) case '\n': case '0': case '#': continue;
        ki = kh_put(p, ret, atoi(buf), &khr);
        kh_val(ret, ki) = atoi(strchr(buf, '|') + 2);
    }
    ki = kh_put(p, ret, 1, &khr);
    kh_val(ret, ki) = 0; // Root of the tree.
    fclose(fp);
    return ret;
}

}
