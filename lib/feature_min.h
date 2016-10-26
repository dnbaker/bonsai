#ifndef _FEATURE_MIN_
#define _FEATURE_MIN_

#include "encoder.h"
#include "spacer.h"
#include "htslib/khash.h"

namespace kpg {


KHASH_SET_INIT_INT64(all)
KHASH_MAP_INIT_INT64(c, uint32_t)

size_t fill_set_genome(const char *path, const Spacer &sp, khash_t(all) *ret, size_t index);


// Return value: whether or not additional sequences were present and added.
int fill_set_seq(kseq_t *ks, const Spacer &sp, khash_t(all) *ret)
{
    assert(ret);
    Encoder<is_lt> enc(0, 0, sp, nullptr);
    int khr; // khash return value. Unused, really.
    if(kseq_read(ks) >= 0) {
        enc.assign(ks);
        while(enc.has_next_kmer()) kh_put(all, ret, enc.next_kmer(), &khr);
        return 1;
    } else return 0;
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

khash_t(c) *feature_count_map(std::vector<std::string> fns, const Spacer &sp, int num_threads=8);

} // namespace kpg
#endif // #ifdef _FEATURE_MIN_
