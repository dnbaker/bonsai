#ifndef _FEATURE_MIN_
#define _FEATURE_MIN_

#include "encoder.h"
#include "spacer.h"
#include "htslib/khash.h"
#include "util.h"

namespace kpg {


void add_to_feature_counter(khash_t(c) *kc, khash_t(all) *set);
khash_t(c) *feature_count_map(std::vector<std::string> fns, const Spacer &sp, int num_threads=8);
uint32_t get_taxid(const char *fn, khash_t(name) *name_hash);

template<uint64_t (*score)(uint64_t, void *)>
khash_t(c) *lca_map(std::vector<std::string> fns, const char *tax_map_path,
                    const char *seq2tax_path,
                    const Spacer &sp, int num_threads=8);
khash_t(c) *make_depth_hash(khash_t(c) *lca_map, khash_t(p) *tax_map);
void lca2depth(khash_t(c) *lca_map, khash_t(p) *tax_map);
template<uint64_t (*score)(uint64_t, void *)>
int fill_set_seq(kseq_t *ks, const Spacer &sp, khash_t(all) *ret);
template<uint64_t (*score)(uint64_t, void *)>
size_t fill_set_genome(const char *path, const Spacer &sp, khash_t(all) *ret, size_t index);
void update_lca_map(khash_t(c) *kc, khash_t(all) *set, khash_t(p) *tax, uint32_t taxid);


// Return value: whether or not additional sequences were present and added.
template<uint64_t (*score)(uint64_t, void *)>
int fill_set_seq(kseq_t *ks, const Spacer &sp, khash_t(all) *ret) {
    assert(ret);
    Encoder<score> enc(0, 0, sp, nullptr);
    int khr; // khash return value. Unused, really.
    if(kseq_read(ks) >= 0) {
        enc.assign(ks);
        while(enc.has_next_kmer()) kh_put(all, ret, enc.next_kmer(), &khr);
        return 1;
    } else return 0;
}

template<uint64_t (*score)(uint64_t, void *)>
size_t fill_set_genome(const char *path, const Spacer &sp, khash_t(all) *ret, size_t index) {
    assert(ret);
    gzFile ifp(gzopen(path, "rb"));
    Encoder<score> enc(0, 0, sp, nullptr);
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

template<uint64_t (*score)(uint64_t, void *)>
khash_t(c) *lca_map(std::vector<std::string> fns, const char *tax_map_path,
                    const char *seq2tax_path,
                    const Spacer &sp, int num_threads) {
    size_t submitted(0), completed(0), todo(fns.size());
    khash_t(all) **counters((khash_t(all) **)malloc(sizeof(khash_t(all) *) * todo));
    khash_t(c) *ret(kh_init(c));
    khash_t(name) *name_hash(build_name_hash(seq2tax_path));
    khash_t(p) *tax_map(load_khash_map<khash_t(p)>(tax_map_path));
    uint32_t taxid;
    for(size_t i(0), end(fns.size()); i != end; ++i) counters[i] = kh_init(all);
    std::vector<std::future<size_t>> futures;

    // Submit the first set of jobs
    for(size_t i(0); i < num_threads && i < todo; ++i) {
        futures.emplace_back(std::async(
          std::launch::async, fill_set_genome<score>, fns[i].data(), sp, counters[i], i));
        ++submitted;
    }

    // Daemon -- check the status of currently running jobs, submit new ones when available.
    while(submitted < todo) {
        for(auto &f: futures) {
            if(is_ready(f)) {
                const size_t index(f.get());
                f = std::async(
                  std::launch::async, fill_set_genome<score>, fns[submitted+1].data(),
                  sp, counters[submitted+1], submitted + 1);
                ++submitted, ++completed;
                update_lca_map(ret, counters[index], tax_map, get_taxid(fns[index].data(), name_hash));
                kh_destroy(all, counters[index]); // Destroy set once we're done with it.
            }
        }
    }

    // Join
    for(auto &f: futures) if(f.valid()) {
        const size_t index(f.get());
        update_lca_map(ret, counters[index], tax_map, get_taxid(fns[index].data(), name_hash));
        kh_destroy(all, counters[index]);
        ++completed;
    }

    // Clean up
    free(counters);
    kh_destroy(p, tax_map);
    destroy_name_hash(name_hash);
    return ret;
}
template <uint64_t (*score)(uint64_t, void *)>
khash_t(c) *feature_count_map(std::vector<std::string> fns, const Spacer &sp, int num_threads) {
    size_t submitted(0), completed(0), todo(fns.size());
    khash_t(all) **counters((khash_t(all) **)malloc(sizeof(khash_t(all) *) * todo));
    khash_t(c) *ret(kh_init(c));
    for(size_t i(0), end(fns.size()); i != end; ++i) counters[i] = kh_init(all);
    std::vector<std::future<size_t>> futures;

    // Submit the first set of jobs
    for(size_t i(0); i < num_threads && i < todo; ++i) {
        futures.emplace_back(std::async(
          std::launch::async, fill_set_genome<score>, fns[i].data(), sp, counters[i], i));
        ++submitted;
    }

    // Daemon -- check the status of currently running jobs, submit new ones when available.
    while(submitted < todo) {
        for(auto &f: futures) {
            if(is_ready(f)) {
                const size_t index(f.get());
                f = std::async(
                  std::launch::async, fill_set_genome<score>, fns[submitted+1].data(),
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



} // namespace kpg
#endif // #ifdef _FEATURE_MIN_
