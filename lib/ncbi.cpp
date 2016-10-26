#include "ncbi.h"
#include <cstdio>
#include <cstdlib>
#include <set>

namespace kpg {

size_t count_lines(const char *fn) {
    FILE *fp(fopen(fn, "r"));
    size_t bufsz = 4096;
    char *buf((char *)malloc(bufsz));
    ssize_t len;
    size_t n(0);
    while((len = getline(&buf, &bufsz, fp)) >= 0) ++n;
    free(buf);
    fclose(fp);
    return n;
}

khash_t(p) *build_parent_map(const char *fn) {
    size_t nlines(count_lines(fn));
    FILE *fp(fopen(fn, "r"));
    khash_t(p) *ret(kh_init(p));
    kh_resize(p, ret, nlines);
    size_t bufsz = 4096;
    char *buf((char *)malloc(bufsz)), *tmp(nullptr);
    ssize_t len;
    int *offsets(nullptr), noffsets(0);
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

// Rewritten from Kraken's source code.
// Consider rewriting to use kbtree instead of std::map.
// lh3's benchmarks indicate that his is only as fast as std::map,
// though it uses less than half the memory. The thing is, taxonomic trees
// aren't that deep.
// We don't gain much from optimizing that.
uint32_t lca(khash_t(p) *map, uint32_t a, uint32_t b) {
    // Use std::set to use RB trees for small set rather than hash table.
    std::set<uint32_t> nodes;
    int khr;
    khint_t ki;
    while(a) {
        nodes.insert(a);
        if((ki = kh_get(p, map, a)) == kh_end(map)) {
            fprintf(stderr, "Missing taxid %i. Abort!\n", a);
            exit(1);
        }
        a = kh_val(map, ki);
    }
    while(b) {
        if(nodes.find(b) != nodes.end()) return b;
        if((ki = kh_get(p, map, b)) == kh_end(map)) {
            fprintf(stderr, "Missing taxid %i. Abort!\n", b);
            exit(1);
        }
        b = kh_val(map, ki);
    }
    return 1;
}

}
