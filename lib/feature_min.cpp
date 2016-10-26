#include "feature_min.h"

namespace kpg {

size_t fill_set_genome(const char *path, const Spacer &sp, khash_t(all) *ret, size_t index)
{
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
            }
        }
    }

    // Join
    for(auto &f: futures) if(f.valid()) add_to_feature_counter(ret, counters[f.get()]), ++completed;

    // Clean up
    for(size_t i(0); i < todo; ++i) kh_destroy(all, counters[i]);
    free(counters);
    return ret;
}

} //namespace kpg
