#include <cctype>
#include "feature_min.h"

namespace kpg {


void lca2depth(khash_t(c) *lca_map, khash_t(p) *tax_map) {
    for(khiter_t ki(kh_begin(lca_map)); ki != kh_end(lca_map); ++ki)
        if(kh_exist(lca_map, ki))
            kh_val(lca_map, ki) = node_depth(tax_map, kh_val(lca_map, ki));
}

khash_t(c) *make_depth_hash(khash_t(c) *lca_map, khash_t(p) *tax_map) {
    khash_t(c) *ret(kh_init(c));
    kh_resize(c, ret, kh_size(lca_map));
    khiter_t ki1;
    int khr;
    for(khiter_t ki2(kh_begin(lca_map)); ki2 != kh_end(lca_map); ++ki2) {
        if(kh_exist(lca_map, ki2)) {
            ki1 = kh_put(c, ret, kh_key(lca_map, ki2), &khr);
            kh_val(ret, ki1) = node_depth(tax_map, kh_val(lca_map, ki2));
        }
    }
    return ret;
}

khash_t(64) *make_taxdepth_hash(khash_t(c) *kc, khash_t(p) *tax) {
    khash_t(64) *ret(kh_init(64));
    int khr;
    khiter_t kir;
    kh_resize(64, ret, kc->n_buckets);
    for(khiter_t ki(0); ki != kh_end(kc); ++ki) {
        if(!kh_exist(kc, ki)) continue;
        kir = kh_put(64, ret, kh_key(kc, ki), &khr);
        kh_val(ret, kir) = ((uint64_t)node_depth(tax, kh_val(kc, ki) << 32)) | kh_val(kc, ki);
    }
    return ret;
}

void update_lca_map(khash_t(c) *kc, khash_t(all) *set, khash_t(p) *tax, uint32_t taxid) {
    int khr;
    khint_t k2;
    for(khiter_t ki(kh_begin(set)); ki != kh_end(set); ++ki) {
        if(kh_exist(set, ki)) {
            if((k2 = kh_get(c, kc, kh_key(set, ki))) == kh_end(kc)) {
                k2 = kh_put(c, kc, kh_key(set, ki), &khr);
                kh_val(kc, k2) = taxid;
            } else if(kh_val(kc, k2) != taxid) kh_val(kc, k2) = lca(tax, taxid, kh_val(kc, k2));
        }
    }
}

uint32_t get_taxid(const char *fn, khash_t(name) *name_hash) {
    gzFile fp(gzopen(fn, "rb"));
    static const size_t bufsz(2048);
    khint_t ki;
    char buf[bufsz];
    char *line(gzgets(fp, buf, bufsz));
    char *p(++line);
    while(!isspace(*p)) ++p;
    *p = 0;
    if((ki = kh_get(name, name_hash, line)) == kh_end(name_hash)) {
        fprintf(stderr, "Missing taxid for %s.\n", line);
        exit(1);
    }
    const uint32_t ret(kh_val(name_hash, ki));
    gzclose(fp);
    return ret;
}


void add_to_feature_counter(khash_t(c) *kc, khash_t(all) *set) {
    int khr;
    khint_t k2;
    for(khiter_t ki(kh_begin(set)); ki != kh_end(set); ++ki) {
        if(kh_exist(set, ki)) {
            if((k2 = kh_get(c, kc, kh_key(set, ki))) == kh_end(kc)) {
                k2 = kh_put(c, kc, kh_key(set, ki), &khr);
                kh_val(kc, k2) = 1;
            } else ++kh_val(kc, k2);
        }
    }
}

} //namespace kpg
