#include "feature_min.h"

namespace kpg {

size_t fill_set_genome(const char *path, const Spacer &sp, khash_t(all) *ret, size_t index) {
    assert(ret);
    gzFile ifp(gzopen(path, "rb"));
    Encoder<is_lt> enc(0, 0, sp, nullptr);
    kseq_t *ks(kseq_init(ifp));
    int khr; // khash return value. Unused, really.
    while(kseq_read(ks) >= 0) {
        enc.assign(ks);
        while(enc.has_next_kmer()) kh_put(all, ret, enc.next_kmer(), &khr);
    }
    kseq_destroy(ks);
    gzclose(ifp);
    return index;
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

khash_t(c) *feature_count_map(std::vector<std::string> fns, const Spacer &sp, int num_threads) {

    size_t submitted(0), completed(0), todo(fns.size());
    khash_t(all) **counters((khash_t(all) **)malloc(sizeof(khash_t(all) *) * todo));
    khash_t(c) *ret(kh_init(c));
    for(size_t i(0), end(fns.size()); i != end; ++i) counters[i] = kh_init(all);
    std::vector<std::future<size_t>> futures;

    // Submit the first set of jobs
    for(size_t i(0); i < num_threads && i < todo; ++i) {
        futures.emplace_back(std::async(
          std::launch::async, fill_set_genome, fns[i].data(), sp, counters[i], i));
        ++submitted;
    }

    // Daemon -- check the status of currently running jobs, submit new ones when available.
    while(submitted < todo) {
        for(auto &f: futures) {
            if(is_ready(f)) {
                const size_t index(f.get());
                f = std::async(
                  std::launch::async, fill_set_genome, fns[submitted+1].data(),
                  sp, counters[submitted+1], submitted + 1);
                ++submitted, ++completed;
                add_to_feature_counter(ret, counters[index]);
                kh_destroy(all, counters[index]); // Destroy set once we're done with it.
            }
        }
    }

    // Join
    for(auto &f: futures) if(f.valid()) {
        const size_t index(f.get());
        add_to_feature_counter(ret, counters[index]);
        kh_destroy(all, counters[index]);
        ++completed;
    }

    // Clean up
    free(counters);
    return ret;
}

} //namespace kpg
