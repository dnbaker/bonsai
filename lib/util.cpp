#include "util.h"
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

// Rewritten from Kraken's source code.
// Consider rewriting to use kbtree instead of std::map.
// lh3's benchmarks indicate that his is only as fast as std::map,
// though it uses less than half the memory. The thing is, taxonomic trees
// aren't that deep.
// We don't gain much from optimizing that.
uint32_t lca(khash_t(p) *map, uint32_t a, uint32_t b) {
    // Use std::set to use RB trees for small set rather than hash table.
    std::set<uint32_t> nodes;
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

unsigned node_depth(khash_t(p) *map, uint32_t a) {
    unsigned ret(0);
    khint_t ki;
    while(a) {
        if((ki = kh_get(p, map, a)) == kh_end(map)) {
            fprintf(stderr, "Tax ID %u missing. Abort!\n", a);
            exit(1);
        }
        a = kh_val(map, ki);
        ++ret;
    }
    return ret;
}

khash_t(name) *build_name_hash(const char *fn) {
    size_t bufsz(2048), namelen;
    char *buf((char *)malloc(bufsz));
    ssize_t len;
    FILE *fp(fopen(fn, "r"));
    char *p;
    int khr;
    khint_t ki;
    khash_t(name) *ret(kh_init(name));
    kh_resize(name, ret, count_lines(fn));
    while((len = getline(&buf, &bufsz, fp)) >= 0) {
        switch(*buf) case '\n': case '#': continue;
        p = strchr(buf, '\t');
        *p = 0;
        ki = kh_put(name, ret, buf, &khr);
        namelen = p - buf;
        kh_key(ret, ki) = (char *)malloc(namelen + 1);
        memcpy((void*)kh_key(ret, ki), buf, namelen);
       ((char *)kh_key(ret, ki))[namelen] = 0;
        kh_val(ret, ki) = atoi(++p);
    }
    fclose(fp);
    return ret;
}

void destroy_name_hash(khash_t(name) *hash) {
    for(khint_t ki(kh_begin(hash)); ki != kh_end(hash); ++ki)
        if(kh_exist(hash, ki))
            free((void *)kh_key(hash, ki));
    kh_destroy(name, hash);
}

}
