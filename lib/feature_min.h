#ifndef _FEATURE_MIN_
#define _FEATURE_MIN_

#include "encoder.h"
#include "spacer.h"
#include "htslib/khash.h"
#include "util.h"

#include <set>

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
khash_t(64) *make_taxdepth_hash(khash_t(c) *kc, khash_t(p) *tax);

// Decode 64-bit hash (contains both tax id and taxonomy depth for id)
#define TDtax(key) ((uint32_t)key)
#define TDdepth(key) (key >> 32)

// Decode 64-bit hash for feature counting.
// TODO: add building of FeatureMin hash
#define FMtax(key) ((uint32_t)key)
#define FMcount(key) (key >> 32)

// Return value: whether or not additional sequences were present and added.
template<uint64_t (*score)(uint64_t, void *)>
int fill_set_seq(kseq_t *ks, const Spacer &sp, khash_t(all) *ret) {
    assert(ret);
    Encoder<score> enc(0, 0, sp, nullptr);
    int khr; // khash return value. Unused, really.
    uint64_t kmer;
    if(kseq_read(ks) >= 0) {
        enc.assign(ks);
        while(enc.has_next_kmer())
            if((kmer = enc.next_minimizer()) != BF)
                kh_put(all, ret, kmer, &khr);
        return 1;
    } else return 0;
}

template<uint64_t (*score)(uint64_t, void *)>
size_t fill_set_genome(const char *path, const Spacer &sp, khash_t(all) *ret, size_t index) {
    assert(ret);
    gzFile ifp(gzopen(path, "rb"));
    if(!ifp) {
        fprintf(stderr, "Could not open file %s. Abort!\n", path);
        exit(1);
    }
    Encoder<score> enc(0, 0, sp, nullptr);
    kseq_t *ks(kseq_init(ifp));
    if(!ks) {
    }
    int khr; // khash return value. Unused, really.
    uint64_t kmer;
    while(kseq_read(ks) >= 0) {
        enc.assign(ks);
        while(enc.has_next_kmer())
            if((kmer = enc.next_minimizer()) != BF)
                kh_put(all, ret, kmer, &khr);
    }
    kseq_destroy(ks);
    gzclose(ifp);
    return index;
}

template<uint64_t (*score)(uint64_t, void *)>
khash_t(64) *taxdepth_map(std::vector<std::string> fns, khash_t(p) *tax_map,
                    const char *seq2tax_path,
                    const Spacer &sp, int num_threads, size_t start_size=1<<10) {
    khash_t(c) *lca(lca_map<score>(fns, tax_map, seq2tax_path, sp, num_threads, start_size));
    khash_t(64) *ret(make_taxdepth_hash(lca, tax_map));
    kh_destroy(c, lca);
    return ret;
}

template<uint64_t (*score)(uint64_t, void *)>
khash_t(c) *lca_map(std::vector<std::string> fns, khash_t(p) *tax_map,
                    const char *seq2tax_path,
                    const Spacer &sp, int num_threads, size_t start_size=1<<10) {
    size_t submitted(0), completed(0), todo(fns.size());
    khash_t(all) **counters((khash_t(all) **)malloc(sizeof(khash_t(all) *) * todo));
    khash_t(c) *ret(kh_init(c));
    kh_resize(c, ret, start_size);
    khash_t(name) *name_hash(build_name_hash(seq2tax_path));
    fprintf(stderr, "%p: name_hash.\n", (void *)name_hash);
    fprintf(stderr, "%p: counters.\n", (void *)counters);
    fprintf(stderr, "%p: ret.\n", (void *)ret);
    /*
    for(khiter_t ki(kh_begin(name_hash)); ki != kh_end(name_hash); ++ki) {
        if(kh_exist(name_hash, ki))
            fprintf(stderr, "Name: %s. Value: %u.\n", kh_key(name_hash, ki), kh_val(name_hash, ki));
    }
    */
    for(size_t i(0), end(fns.size()); i != end; ++i) counters[i] = kh_init(all);
    std::vector<std::future<size_t>> futures;

    // Submit the first set of jobs
    std::set<size_t> subbed, used;
    for(int i(0); i < num_threads && i < (ssize_t)todo; ++i) {
        futures.emplace_back(std::async(
          std::launch::async, fill_set_genome<score>, fns[i].data(), sp, counters[i], i));
        LOG_DEBUG("Submitted for %zu.\n", submitted);
        subbed.insert(submitted);
        ++submitted;
    }

    // Daemon -- check the status of currently running jobs, submit new ones when available.
    while(submitted < todo) {
        LOG_DEBUG("Submitted %zu, todo %zu\n", submitted, todo);
        for(auto &f: futures) {
            if(is_ready(f)) {
                if(submitted == todo) break;
                const size_t index(f.get());
                if(used.find(index) != used.end()) continue;
                used.insert(index);
                if(subbed.find(submitted) != subbed.end()) throw "a party!";
                f = std::async(
                  std::launch::async, fill_set_genome<score>, fns[submitted].data(),
                  sp, counters[submitted], submitted);
                subbed.insert(submitted);
                LOG_DEBUG("Submitted for %zu. Updating map for %zu\n", submitted, index);
                ++submitted, ++completed;
                update_lca_map(ret, counters[index], tax_map, get_taxid(fns[index].data(), name_hash));
                kh_destroy(all, counters[index]); // Destroy set once we're done with it.
            }
        }
    }

    // Join
    for(auto &f: futures) if(f.valid()) {
        const size_t index(f.get());
        if(used.find(index) != used.end()) continue;
        used.insert(index);
        update_lca_map(ret, counters[index], tax_map, get_taxid(fns[index].data(), name_hash));
        kh_destroy(all, counters[index]);
        ++completed;
    }
    LOG_DEBUG("Finished LCA map building! Subbed %zu, completed %zu, size of futures %zu.\n", submitted, completed, used.size());
#if !NDEBUG    
    for(size_t i(0); i < todo; ++i) assert(used.find(i) != used.end());
#endif

    // Clean up
    free(counters);
    destroy_name_hash(name_hash);
    LOG_DEBUG("Cleaned up after LCA map building!\n")
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
    std::set<size_t> used;
    for(size_t i(0); i < num_threads && i < todo; ++i) {
        futures.emplace_back(std::async(
          std::launch::async, fill_set_genome<score>, fns[i].data(), sp, counters[i], i));
        LOG_DEBUG("Submitted for %zu.\n", submitted);
        ++submitted;
    }

    // Daemon -- check the status of currently running jobs, submit new ones when available.
    while(submitted < todo) {
        for(auto &f: futures) {
            if(is_ready(f)) {
                const size_t index(f.get());
                if(submitted == todo) break;
                if(used.find(index) != used.end()) continue;
                used.insert(index);
                LOG_DEBUG("Submitted for %zu.\n", submitted);
                f = std::async(
                  std::launch::async, fill_set_genome<score>, fns[submitted].data(),
                  sp, counters[submitted], submitted);
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
