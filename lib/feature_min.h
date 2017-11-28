#ifndef _FEATURE_MIN__
#define _FEATURE_MIN__

#include "encoder.h"
#include "spacer.h"
#include "khash64.h"
#include "util.h"
#include "klib/kthread.h"
#include <set>

// Decode 64-bit hash (contains both tax id and taxonomy depth for id)
#define TDtax(key) ((tax_t)key)
#define TDdepth(key) ((tax_t)~0 - (key >> 32))
#define TDencode(depth, taxid) (((u64)((tax_t)~0 - depth) << 32) | taxid)

// Decode 64-bit hash for feature counting.
// TODO: add building of FeatureMin hash
#define FMtax(key) ((tax_t)key)
#define FMcount(key) (key >> 32)

#define FMencode(count, taxid) (((u64)count << 32) | taxid)


namespace emp {


template<u64 (*score)(u64, void *)>
khash_t(64) *feature_count_map(std::vector<std::string> fns, const Spacer &sp, int num_threads=8);

khash_t(c) *make_depth_hash(khash_t(c) *lca_map, khash_t(p) *tax_map);
void lca2depth(khash_t(c) *lca_map, khash_t(p) *tax_map);
template<u64 (*score)(u64, void *)>
int fill_set_seq(kseq_t *ks, const Spacer &sp, khash_t(all) *ret);

void update_lca_map(khash_t(c) *kc, khash_t(all) *set, khash_t(p) *tax, tax_t taxid, std::shared_mutex &m);
void update_td_map(khash_t(64) *kc, khash_t(all) *set, khash_t(p) *tax, tax_t taxid);
khash_t(64) *make_taxdepth_hash(khash_t(c) *kc, khash_t(p) *tax);
void update_feature_counter(khash_t(64) *kc, khash_t(all) *set, khash_t(p) *tax, const tax_t taxid);
void update_minimized_map(khash_t(all) *set, khash_t(64) *full_map, khash_t(c) *ret);

#define next_minimizer next_kmer

// Return value: whether or not additional sequences were present and added.
template<u64 (*score)(u64, void *)>
int fill_set_seq(kseq_t *ks, const Spacer &sp, khash_t(all) *ret) {
    assert(ret);
    Encoder<score> enc(0, 0, sp, nullptr);
    int khr; // khash return value. Unused, really.
    u64 kmer;
    if(kseq_read(ks) >= 0) {
        enc.assign(ks);
        while(enc.has_next_kmer())
            if((kmer = enc.next_minimizer()) != BF)
                kh_put(all, ret, kmer, &khr);
        return 1;
    } else return 0;
}

template<u64 (*score)(u64, void *)>
size_t fill_set_genome(const char *path, const Spacer &sp, khash_t(all) *ret, size_t index, void *data) {
    LOG_ASSERT(ret);
    LOG_INFO("Filling from genome at path %s\n", path);
    gzFile ifp(gzopen(path, "rb"));
    if(!ifp) {
        fprintf(stderr, "Could not open file %s for index %zu. Abort!\n", path, index);
        std::exit(EXIT_FAILURE);
    }
    Encoder<score> enc(0, 0, sp, data);
    kseq_t *ks(kseq_init(ifp));
    int khr; // khash return value. Unused, really.
    u64 kmer;
    if(sp.w_ > sp.k_) {
        while(kseq_read(ks) >= 0) {
            enc.assign(ks);
            while(likely(enc.has_next_kmer())) {
                if((kmer = enc.next_minimizer()) != BF)
                    kh_put(all, ret, kmer, &khr);
            }
        }
    } else {
        while(kseq_read(ks) >= 0) {
            enc.assign(ks);
            while(likely(enc.has_next_kmer()))
                if((kmer = enc.next_kmer()) != BF)
                    kh_put(all, ret, kmer, &khr);
        }
    }
    kseq_destroy(ks);
    gzclose(ifp);
    LOG_INFO("Set of size %lu filled from genome at path %s\n", kh_size(ret), path);
    return index;
}

#undef next_minimizer

template<u64 (*score)(u64, void *), class Container>
size_t fill_set_genome_container(Container &container, const Spacer &sp, khash_t(all) *ret, void *data) {
    size_t sz(0);
    for(std::string &str: container)
        sz += fill_set_genome<score>(str.data(), sp, ret, 0, data);
    return sz;
}

template<u64 (*score)(u64, void *)>
khash_t(64) *ftct_map(std::vector<std::string> &fns, khash_t(p) *tax_map,
                      const char *seq2tax_path,
                      const Spacer &sp, int num_threads, size_t start_size) {
    return feature_count_map<score>(fns, tax_map, seq2tax_path, sp, num_threads, start_size);
}

#if USE_PAR_HELPERS
template<typename T>
struct helper {
    const Spacer                   &sp_;
    const std::vector<std::string> &fns_;
    T                              *ret_;
    khash_t(name)                  *name_hash_;
    khash_t(p)                     *tax_map_;
    std::vector<khash_t(all) *>     hashes_;
    void                           *data_;
    std::shared_mutex              *m_;
};
using lca_helper = helper<khash_t(c)>;
using td_helper  = helper<khash_t(64)>;
template<u64 (*score)(u64, void *)>
void fct_for_helper(void *data_, long index, int tid) {
    LOG_DEBUG("Helper with index %ld and tid %i starting\n", index, tid);
    td_helper &tdh(*(td_helper *)data_);
    khash_t(all) *h(tdh.hashes_[tid]);
    khash_t(64) *ret(tdh.ret_);
    if(kh_size(h)) kh_clear(all, h);
    fill_set_genome<score>(tdh.fns_[index].data(), tdh.sp_, h, index, nullptr);
    tax_t taxid(get_taxid(tdh.fns_[index].data(), tdh.name_hash_));
    if(taxid == UINT32_C(-1)) {
        LOG_WARNING("Taxid for %s not listed in summary.txt. Not including.", tdh.fns_[index].data());
    } else update_feature_counter(ret, h, tdh.tax_map_, taxid);
}

template<u64 (*score)(u64, void *)>
void td_for_helper(void *data_, long index, int tid) {
    LOG_DEBUG("Helper with index %ld and tid %i starting\n", index, tid);
    td_helper &tdh(*(td_helper *)data_);
    khash_t(all) *h(tdh.hashes_[tid]);
    khash_t(64) *ret(tdh.ret_);
    if(kh_size(h)) kh_clear(all, h);
    fill_set_genome<score>(tdh.fns_[index].data(), tdh.sp_, h, index, nullptr);
    tax_t taxid(get_taxid(tdh.fns_[index].data(), tdh.name_hash_));
    if(taxid == UINT32_C(-1)) {
        LOG_WARNING("Taxid for %s not listed in summary.txt. Not including.", tdh.fns_[index].data());
    } else update_td_map(ret, h, tdh.tax_map_, taxid);
}

template<u64 (*score)(u64, void *)>
void lca_for_helper(void *data_, long index, int tid) {
    LOG_DEBUG("Helper with index %ld and tid %i starting\n", index, tid);
    lca_helper &lh(*(lca_helper *)data_);
    khash_t(all) *h(lh.hashes_[tid]);
    khash_t(c) *ret(lh.ret_);
    if(kh_size(h)) kh_clear(all, h);
    fill_set_genome<score>(lh.fns_[index].data(), lh.sp_, h, index, nullptr);
    tax_t taxid(get_taxid(lh.fns_[index].data(), lh.name_hash_));
    if(taxid == UINT32_C(-1)) {
        LOG_WARNING("Taxid for %s not listed in summary.txt. Not including.", lh.fns_[index].data());
    } else update_lca_map(ret, h, lh.tax_map_, taxid, *lh.m_);
}

template<u64 (*score)(u64, void *)>
void min_for_helper(void *data_, long index, int tid) {
    LOG_DEBUG("Helper with index %ld and tid %i starting\n", index, tid);
    lca_helper &lh(*(lca_helper *)data_);
    khash_t(all) *h(lh.hashes_[tid]);
    khash_t(c) *ret(lh.ret_);
    if(kh_size(h)) kh_clear(all, h);
    fill_set_genome<score>(lh.fns_[index].data(), lh.sp_, h, index, (khash_t(64) *)lh.data_);
    update_minimized_map(h, (khash_t(64) *)lh.data_, ret);
}

template<u64 (*score)(u64, void *)>
khash_t(c) *minimized_map(std::vector<std::string> fns,
                          khash_t(64) *full_map,
                          const Spacer &sp, int num_threads, size_t start_size) {
    khash_t(c) *ret(kh_init(c));
    kh_resize(c, ret, start_size);
    std::vector<khash_t(all) *> counters;
    counters.reserve(num_threads);
    for(size_t i(0), end(fns.size()); i != end; ++i) counters.emplace_back(kh_init(all));
    lca_helper lh{sp, fns, ret, nullptr, nullptr, counters, (void *)full_map};
    kt_for(num_threads, &min_for_helper<score>, (void *)&lh, fns.size());
    for(auto c: counters) khash_destroy(c);
    LOG_DEBUG("Cleaned up after minimized map building! Size: %zu\n", kh_size(ret));
    return ret;
}


template<u64 (*score)(u64, void *)>
khash_t(64) *taxdepth_map(std::vector<std::string> &fns, khash_t(p) *tax_map,
                          const char *seq2tax_path, const Spacer &sp,
                          int num_threads, size_t start_size) {
    khash_t(64) *ret(kh_init(64));
    kh_resize(64, ret, start_size);
    khash_t(name) *name_hash(build_name_hash(seq2tax_path));
    std::vector<khash_t(all) *> counters;
    counters.reserve(num_threads);
    for(size_t i(0), end(fns.size()); i != end; ++i) counters.emplace_back(kh_init(all));
    td_helper th{sp, fns, ret, name_hash, tax_map, counters, nullptr};
    kt_for(num_threads, &td_for_helper<score>, (void *)&th, fns.size());
    destroy_name_hash(name_hash);
    for(auto c: counters) khash_destroy(c);
    LOG_DEBUG("Cleaned up after taxdepth map building!\n");
    return ret;
}


template<u64 (*score)(u64, void *)>
khash_t(c) *lca_map(const std::vector<std::string> &fns, khash_t(p) *tax_map,
                    const char *seq2tax_path,
                    const Spacer &sp, int num_threads, size_t start_size) {
    khash_t(c) *ret(kh_init(c));
    kh_resize(c, ret, start_size);
    khash_t(name) *name_hash(build_name_hash(seq2tax_path));
    std::vector<khash_t(all) *> counters;
    std::shared_mutex m;
    counters.reserve(num_threads);
    for(size_t i(0), end(fns.size()); i != end; ++i) counters.emplace_back(kh_init(all));
    lca_helper lh{sp, fns, ret, name_hash, tax_map, counters, nullptr, &m};
    kt_for(num_threads, &lca_for_helper<score>, (void *)&lh, fns.size());
    destroy_name_hash(name_hash);
    for(auto c: counters) khash_destroy(c);
    LOG_DEBUG("Cleaned up after LCA map building! Size: %zu\n", kh_size(ret));
    return ret;
}


template <u64 (*score)(u64, void *)>
khash_t(64) *feature_count_map(std::vector<std::string> fns, khash_t(p) *tax_map,
                               const char *seq2tax_path,
                               const Spacer &sp, int num_threads, size_t start_size) {
    // TODO: make scale better.
    khash_t(64) *ret(kh_init(64));
    kh_resize(64, ret, start_size);
    khash_t(name) *name_hash(build_name_hash(seq2tax_path));
    std::vector<khash_t(all) *> counters;
    counters.reserve(num_threads);
    for(size_t i(0), end(fns.size()); i != end; ++i) counters.emplace_back(kh_init(all));
    td_helper th{sp, fns, ret, name_hash, tax_map, counters, nullptr};
    kt_for(num_threads, &fct_for_helper<score>, (void *)&th, fns.size());
    for(auto c: counters) khash_destroy(c);
    destroy_name_hash(name_hash);
    return ret;
}
#else
template<u64 (*score)(u64, void *)>
khash_t(c) *minimized_map(std::vector<std::string> fns,
                          khash_t(64) *full_map,
                          const Spacer &sp, int num_threads, size_t start_size) {
    size_t submitted(0), completed(0), todo(fns.size());
    std::vector<khash_t(all) *> counters(todo, nullptr);
    khash_t(c) *ret(kh_init(c));
    kh_resize(c, ret, start_size);
    std::vector<std::future<size_t>> futures;
    //for(auto &i: fns) fprintf(stderr, "Filename: %s\n", i.data());

    if(num_threads < 0) num_threads = 16;

    LOG_DEBUG("Number of items to do: %zu\n", todo);

    for(size_t i(0); i < todo; ++i) counters[i] = kh_init(all), kh_resize(all, counters[i], start_size);

    // Submit the first set of jobs
    for(int i(0), e(std::min(num_threads, (int)todo)); i < e; ++i) {
        futures.emplace_back(std::async(
          std::launch::async, fill_set_genome<score>, fns[i].data(), sp, counters[i], i, (void *)full_map));
        ++submitted;
    }

    // Daemon -- check the status of currently running jobs, submit new ones when available.
    while(submitted < todo) {
        //LOG_DEBUG("Submitted %zu, todo %zu\n", submitted, todo);
        for(auto f(futures.begin()), fend(futures.end()); f != fend; ++f) {
            if(is_ready(*f)) {
                const size_t index(f->get());
                futures.erase(f);
                futures.emplace_back(std::async(
                     std::launch::async, fill_set_genome<score>, fns[submitted].data(),
                     sp, counters[submitted], submitted, (void *)full_map));
                LOG_INFO("Submitted for %zu. Updating map for %zu. Total completed/all: %zu/%zu. Current size: %zu\n",
                         submitted, index, completed, todo, kh_size(ret));
                ++submitted, ++completed;
                update_minimized_map(counters[index], full_map, ret);
                kh_destroy(all, counters[index]); // Destroy set once we're done with it.
                break;
            }
        }
    }

    // Join
    for(auto &f: futures) if(f.valid()) {
        const size_t index(f.get());
        update_minimized_map(counters[index], full_map, ret);
        kh_destroy(all, counters[index]);
        ++completed;
        LOG_DEBUG("Number left to do: %zu\n", todo - completed);
    }
    LOG_DEBUG("Finished minimized map building! Subbed %zu, completed %zu.\n", submitted, completed);

    // Clean up
    LOG_DEBUG("Cleaned up after LCA map building!\n");
    return ret;
}

template<u64 (*score)(u64, void *)>
khash_t(64) *taxdepth_map(std::vector<std::string> &fns, khash_t(p) *tax_map,
                          const char *seq2tax_path, const Spacer &sp,
                          int num_threads, size_t start_size=1<<10) {
    size_t submitted(0), completed(0), todo(fns.size());
    std::vector<khash_t(all) *> counters(todo, nullptr);
    khash_t(64) *ret(kh_init(64));
    kh_resize(64, ret, start_size);
    khash_t(name) *name_hash(build_name_hash(seq2tax_path));
    std::vector<std::future<size_t>> futures;

    for(size_t i(0), end(fns.size()); i != end; ++i) counters[i] = kh_init(all);

    // Submit the first set of jobs
    std::set<size_t> subbed, used;
    for(int i(0), e(std::min(num_threads, (int)todo)); i < e; ++i) {
        LOG_DEBUG("Launching thread to read from file %s.\n", fns[i].data());
        futures.emplace_back(std::async(
          std::launch::async, fill_set_genome<score>, fns[i].data(), sp, counters[i], i, nullptr));
        //LOG_DEBUG("Submitted for %zu.\n", submitted);
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
                if(subbed.find(submitted) != subbed.end()) throw std::runtime_error("Could not find what I was looking for....");
                LOG_DEBUG("Launching thread to read from file %s.\n", fns[submitted].data());
                f = std::async(
                  std::launch::async, fill_set_genome<score>, fns[submitted].data(),
                  sp, counters[submitted], submitted, nullptr);
                subbed.insert(submitted);
                LOG_INFO("Submitted for %zu. Updating map for %zu. Total completed/all: %zu/%zu. Total size: %zu.\n",
                         submitted, index, completed, todo, kh_size(ret));
                ++submitted, ++completed;
                const tax_t taxid(get_taxid(fns[index].data(), name_hash));
                LOG_DEBUG("Just fetched taxid from file %s %u.\n", fns[index].data(), taxid);
                update_td_map(ret, counters[index], tax_map, taxid);
                kh_destroy(all, counters[index]); // Destroy set once we're done with it.
            }
        }
    }

    // Join
    for(auto &f: futures) if(f.valid()) {
        const size_t index(f.get());
        if(used.find(index) != used.end()) continue;
        used.insert(index);
        update_td_map(ret, counters[index], tax_map, get_taxid(fns[index].data(), name_hash));
        kh_destroy(all, counters[index]);
        ++completed;
        LOG_DEBUG("Number left to do: %zu\n", todo - completed);
    }
    LOG_DEBUG("Finished LCA map building! Subbed %zu, completed %zu, size of futures %zu.\n", submitted, completed, used.size());
#if !NDEBUG
    for(size_t i(0); i < todo; ++i) assert(used.find(i) != used.end());
#endif

    // Clean up
    destroy_name_hash(name_hash);
    LOG_DEBUG("Cleaned up after LCA map building!\n");
    return ret;
}

template<u64 (*score)(u64, void *)>
khash_t(c) *lca_map(const std::vector<std::string> &fns, khash_t(p) *tax_map,
                    const char *seq2tax_path,
                    const Spacer &sp, int num_threads, size_t start_size=1<<10) {
    size_t submitted(0), completed(0), todo(fns.size());
    std::vector<khash_t(all) *> counters;
    khash_t(c) *ret(kh_init(c));
    kh_resize(c, ret, start_size);
    counters.reserve(todo);
    khash_t(name) *name_hash(build_name_hash(seq2tax_path));
    std::vector<std::future<size_t>> futures;
    std::shared_mutex m;

    for(size_t i(0), end(fns.size()); i != end; ++i) counters.emplace_back(kh_init(all));


    // Submit the first set of jobs
    std::set<size_t> subbed, used;
    for(int i(0), e(std::min(num_threads, (int)todo)); i < e; ++i) {
        LOG_DEBUG("Launching thread to read from file %s.\n", fns[i].data());
        futures.emplace_back(std::async(
          std::launch::async, fill_set_genome<score>, fns[i].data(), sp, counters[i], i, nullptr));
        subbed.insert(submitted);
        ++submitted;
    }

    // Daemon -- check the status of currently running jobs, submit new ones when available.
    while(submitted < todo) {
        LOG_DEBUG("Submitted %zu, todo %zu\n", submitted, todo);
        tax_t taxid;
        for(auto &f: futures) {
            if(is_ready(f)) {
                if(submitted == todo) break;
                const size_t index(f.get());
                if(used.find(index) != used.end()) continue;
                used.insert(index);
                if(subbed.find(submitted) != subbed.end()) throw std::runtime_error("Could not find what I was looking for....");
                LOG_DEBUG("Launching thread to read from file %s.\n", fns[submitted].data());
                f = std::async(
                  std::launch::async, fill_set_genome<score>, fns[submitted].data(),
                  sp, counters[submitted], submitted, nullptr);
                subbed.insert(submitted);
                LOG_INFO("Submitted for %zu. Updating map for %zu. Total completed/all: %zu/%zu. Total size: %zu\n",
                         submitted, index, completed, todo, kh_size(ret));
                ++submitted, ++completed;
                if((taxid = get_taxid(fns[index].data(), name_hash)) == UINT32_C(-1)) {
                    LOG_WARNING("Taxid for %s not listed in summary.txt. Not including.\n", fns[index].data());
                } else update_lca_map(ret, counters[index], tax_map, taxid, m);
                kh_destroy(all, counters[index]); // Destroy set once we're done with it.
            }
        }
    }

    // Join
    for(auto &f: futures) if(f.valid()) {
        const size_t index(f.get());
        if(used.find(index) != used.end()) continue;
        used.insert(index);
        tax_t taxid(get_taxid(fns[index].data(), name_hash));
        if(taxid == UINT32_C(-1)) {
            LOG_WARNING("Taxid for %s not listed in summary.txt. Not including.", fns[index].data());
        } else update_lca_map(ret, counters[index], tax_map, taxid, m);
        kh_destroy(all, counters[index]);
        ++completed;
        LOG_DEBUG("Number left to do: %zu\n", todo - completed);
    }
    LOG_DEBUG("Finished LCA map building! Subbed %zu, completed %zu, size of futures %zu.\n", submitted, completed, used.size());

    // Clean up
    destroy_name_hash(name_hash);
    LOG_DEBUG("Cleaned up after LCA map building!\n");
    return ret;
}

template <u64 (*score)(u64, void *)>
khash_t(64) *feature_count_map(std::vector<std::string> fns, khash_t(p) *tax_map, const char *seq2tax_path, const Spacer &sp, int num_threads, size_t start_size) {
    // Update this to include tax ids in the hash map.
    size_t submitted(0), completed(0), todo(fns.size());
    std::vector<khash_t(all) *> counters(todo, nullptr);
    khash_t(64) *ret(kh_init(64));
    kh_resize(64, ret, start_size);
    khash_t(name) *name_hash(build_name_hash(seq2tax_path));
    for(size_t i(0), end(fns.size()); i != end; ++i) counters[i] = kh_init(all);
    std::vector<std::future<size_t>> futures;
    fprintf(stderr, "Will use tax_map (%p) and seq2tax_map (%s) to assign "
                    "feature-minimized values to all kmers.\n", (void *)tax_map, seq2tax_path);

    // Submit the first set of jobs
    std::set<size_t> used;
    for(size_t i(0); i < (unsigned)num_threads && i < todo; ++i) {
        futures.emplace_back(std::async(
          std::launch::async, fill_set_genome<score>, fns[i].data(), sp, counters[i], i, nullptr));
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
                  sp, counters[submitted], submitted, nullptr);
                ++submitted, ++completed;
                const tax_t taxid(get_taxid(fns[index].data(), name_hash));
                update_feature_counter(ret, counters[index], tax_map, taxid);
                kh_destroy(all, counters[index]); // Destroy set once we're done with it.
            }
        }
    }

    // Join
    for(auto &f: futures) if(f.valid()) {
        const size_t index(f.get());
        const tax_t taxid(get_taxid(fns[index].data(), name_hash));
        update_feature_counter(ret, counters[index], tax_map, taxid);
        kh_destroy(all, counters[index]);
        ++completed;
    }

    // Clean up
    kh_destroy(name, name_hash);
    return ret;
}
#endif


} // namespace emp
#endif // #ifdef _FEATURE_MIN__

