#pragma once

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


namespace bns {


inline khash_t(c) *make_depth_hash(khash_t(c) *lca_map, khash_t(p) *tax_map);
inline void lca2depth(khash_t(c) *lca_map, khash_t(p) *tax_map);

inline khash_t(64) *make_taxdepth_hash(khash_t(c) *kc, const khash_t(p) *tax);


inline void update_lca_map(khash_t(c) *kc, const khash_t(all) *set, const khash_t(p) *tax, tax_t taxid);
inline void update_td_map(khash_t(64) *kc, const khash_t(all) *set, const khash_t(p) *tax, tax_t taxid);
inline void update_feature_counter(khash_t(64) *kc, const khash_t(all) *set, const khash_t(p) *tax, tax_t taxid);
inline void update_minimized_map(const khash_t(all) *set, const khash_t(64) *full_map, khash_t(c) *ret);

// Wrap these in structs so that downstream code can be managed as a set, not updated one-by-one.
struct LcaMap {
    using ReturnType = khash_t(c) *;
    static constexpr size_t ValSize = sizeof(*(ReturnType{0})->vals);
    static void update(const khash_t(p) *tax, const khash_t(all) *set, const khash_t(64) *d64, khash_t(c) *r32, khash_t(64) *r64, tax_t taxid) {
        update_lca_map(r32, set, tax, taxid);
    }
};
struct TdMap {
    using ReturnType = khash_t(64) *;
    static constexpr size_t ValSize = sizeof(*(ReturnType{0})->vals);
    static void update(const khash_t(p) *tax, const khash_t(all) *set, const khash_t(64) *d64, khash_t(c) *r32, khash_t(64) *r64, tax_t taxid) {
        update_td_map(r64, set, tax, taxid);
    }
};
struct FcMap {
    using ReturnType = khash_t(64) *;
    static constexpr size_t ValSize = sizeof(*(ReturnType{0})->vals);
    static void update(const khash_t(p) *tax, const khash_t(all) *set, const khash_t(64) *d64, khash_t(c) *r32, khash_t(64) *r64, tax_t taxid) {
        update_feature_counter(r64, set, tax, taxid);
    }
};
struct MinMap {
    using ReturnType = khash_t(c) *;
    static constexpr size_t ValSize = sizeof(*(ReturnType{0})->vals);
    static void update(const khash_t(p) *tax, const khash_t(all) *set, const khash_t(64) *d64, khash_t(c) *r32, khash_t(64) *r64, tax_t taxid) {
        update_minimized_map(set, d64, r32);
    }
};

template<typename ScoreType>
size_t fill_set_genome(const char *path, const Spacer &sp, khash_t(all) *ret, size_t index, void *data, bool canon, kseq_t *ks=nullptr) {
    LOG_ASSERT(ret);
    LOG_DEBUG("Filling from genome at path %s. kseq is pre-allocated ? %s. %p\n", path, ks ? "true": "false", (void *)ks);

    Encoder<ScoreType> enc(0, 0, sp, data, canon);
    enc.add(ret, path, ks);
    LOG_DEBUG("Set of size %lu filled from genome at path %s\n", kh_size(ret), path);
    return index;
}

template<typename Container, typename ScoreType>
size_t fill_set_genome_container(Container &container, const Spacer &sp, khash_t(all) *ret, void *data, bool canon, kseq_t *ks=nullptr) {
    bool destroy;
    size_t sz(0);
    for(std::string &str: container)
        sz += fill_set_genome<ScoreType>(str.data(), sp, ret, 0, data, canon);
    return sz;
}


template<typename ScoreType, typename MapUpdater>
typename MapUpdater::ReturnType
make_map(const std::vector<std::string> fns, const khash_t(p) *tax_map, const char *seq2tax_path, const Spacer &sp, int num_threads, bool canon, size_t start_size, const khash_t(64) *data) {
    MapUpdater mu;

    khash_t(c) *r32 = nullptr;
    khash_t(64) *r64 = nullptr;
    // Update this to include tax ids in the hash map.
    size_t submitted(0), completed(0), todo(fns.size());
    std::vector<khash_t(all)> counters(num_threads);
    
#if 0
		khint_t n_buckets, size, n_occupied, upper_bound; \
		khint32_t *flags; \
		khkey_t *keys; \
		khval_t *vals; \
		mutable std::shared_mutex m;
#endif
    std::memset(counters.data(), 0, sizeof(khash_t(all)) * counters.size());
    LOG_DEBUG("Started things\n");
    if constexpr(MapUpdater::ValSize == 8) {
        r64 = static_cast<typename MapUpdater::ReturnType>(std::calloc(sizeof(khash_t(64)), 1));
        kh_resize(64, r64, start_size);
    } else {
        r32 = static_cast<typename MapUpdater::ReturnType>(std::calloc(sizeof(khash_t(c)), 1));
        kh_resize(c, r32, start_size);
    }
    LOG_DEBUG("Allocated memory\n");
    khash_t(name) *name_hash(build_name_hash(seq2tax_path));
    std::vector<std::future<size_t>> futures;
    // Mkae the future return the kseq pointer and then use it for resubmission.
    // TODO: Also use a fixed st of kh_all sets to reduce memory allocations.
    KSeqBufferHolder kseqs(num_threads);
    std::vector<uint32_t> counter_map;

    // Submit the first set of jobs
    std::set<size_t> used;
    for(size_t i(0); i < (unsigned)num_threads && i < todo; ++i) {
        futures.emplace_back(std::async(
          std::launch::async, fill_set_genome<ScoreType>, fns[i].data(), sp, counters.data() + i, i, (void *)data, canon, kseqs.data() + submitted));
        counter_map.emplace_back(submitted);
        LOG_DEBUG("Submitted for %zu.\n", submitted);
        ++submitted;
    }

    // Daemon -- check the status of currently running jobs, submit new ones when available.
    while(submitted < todo) {
        if(auto it = std::find_if(std::begin(futures), std::end(futures), [](auto &f) {return is_ready(f);}); it != futures.end()) {
            auto &f(*it);
            const size_t index(f.get());
            if(submitted == todo) break;
            if(used.find(index) != used.end()) continue;
            used.insert(index);
            const auto coffset = counter_map.at(index);
            khash_t(all) *counter = counters.data() + coffset; // Pointer to the counter to use
            kseq_t *ks_to_submit = kseqs.data() + coffset;
            f = std::async(
              std::launch::async, fill_set_genome<ScoreType>, fns[submitted].data(),
              sp, counter, submitted, (void *)data, canon, ks_to_submit);
            counter_map.emplace_back(coffset);
            ++submitted, ++completed;
            LOG_DEBUG("Have now submitted %zu element\n", submitted);
            const tax_t taxid(get_taxid(fns[index].data(), name_hash));
            mu.update(tax_map, counter, data, r32, r64, taxid);
        }
    }

    // Join
    for(auto &f: futures) if(f.valid()) {
        const size_t index(f.get());
        const tax_t taxid(get_taxid(fns[index].data(), name_hash));
        mu.update(tax_map, counters.data() + index, data, r32, r64, taxid);
        ++completed;
    }

    // Clean up
    for(auto &counter: counters) {
        std::free(counter.flags);
        std::free(counter.keys);
    }
    kh_destroy(name, name_hash);
    LOG_DEBUG("Finished making map!\n");
    if constexpr(MapUpdater::ValSize == 8)
        return r64;
    else
        return r32;
}
template<typename ScoreType>
auto feature_count_map(const std::vector<std::string> fns, const khash_t(p) *tax_map, const char *seq2tax_path, const Spacer &sp, int num_threads, bool canon, size_t start_size) {
    return make_map<ScoreType, FcMap>(fns, tax_map, seq2tax_path, sp, num_threads, canon, start_size, nullptr);
}

template<typename ScoreType>
khash_t(c) *lca_map(const std::vector<std::string> &fns, const khash_t(p) *tax_map,
                    const char *seq2tax_path,
                    const Spacer &sp, int num_threads, bool canon, size_t start_size) {
    return make_map<ScoreType, LcaMap>(fns, tax_map, seq2tax_path, sp, num_threads, canon, start_size, nullptr);
}

template<typename ScoreType>
khash_t(c) *minimized_map(std::vector<std::string> fns,
                          const khash_t(64) *full_map, const char *seq2tax_path, const khash_t(p) *tax_map,
                          const Spacer &sp, int num_threads, size_t start_size, bool canon) {
    return make_map<ScoreType, MinMap>(fns, tax_map, seq2tax_path, sp, num_threads, canon, start_size, full_map);
}

template<typename ScoreType>
khash_t(64) *ftct_map(const std::vector<std::string> &fns, const khash_t(p) *tax_map,
                      const char *seq2tax_path,
                      const Spacer &sp, int num_threads, bool canon, size_t start_size) {
    return feature_count_map<ScoreType>(fns, tax_map, seq2tax_path, sp, num_threads, canon, start_size);
}
template<typename ScoreType>
khash_t(64) *taxdepth_map(const std::vector<std::string> &fns, const khash_t(p) *tax_map,
                          const char *seq2tax_path, const Spacer &sp,
                          int num_threads, bool canon, size_t start_size=1<<10) {
    return make_map<ScoreType, TdMap>(fns, tax_map, seq2tax_path, sp, num_threads, canon, start_size, nullptr);
}

inline void update_lca_map(khash_t(c) *kc, const khash_t(all) *set, const khash_t(p) *tax, tax_t taxid) {
    int khr;
    khint_t k2;
    LOG_DEBUG("Adding set of size %zu to total set of current size %zu.\n", kh_size(set), kh_size(kc));
    for(khiter_t ki(kh_begin(set)); ki != kh_end(set); ++ki) {
        if(kh_exist(set, ki)) {
            if((k2 = kh_get(c, kc, kh_key(set, ki))) == kh_end(kc)) {
                k2 = kh_put(c, kc, kh_key(set, ki), &khr);
                kh_val(kc, k2) = taxid;
#if !NDEBUG
                if(unlikely(kh_size(kc) % 1000000 == 0)) LOG_DEBUG("Final hash size %zu\n", kh_size(kc));
#endif
            } else if(kh_val(kc, k2) != taxid) {
                kh_val(kc, k2) = lca(tax, taxid, kh_val(kc, k2));
                if(kh_val(kc, k2) == UINT32_C(-1)) kh_val(kc, k2) = 1, LOG_WARNING("Missing taxid %u. Setting lca to 1\n", taxid);
            }
        }
    }
    LOG_DEBUG("After updating with set of size %zu, total set current size is %zu.\n", kh_size(set), kh_size(kc));
}

inline void update_td_map(khash_t(64) *kc, const khash_t(all) *set, const khash_t(p) *tax, tax_t taxid) {
    int khr;
    khint_t k2;
    tax_t val;
    LOG_DEBUG("Adding set of size %zu to total set of current size %zu.\n", kh_size(set), kh_size(kc));
    for(khiter_t ki(kh_begin(set)); ki != kh_end(set); ++ki) {
        if(kh_exist(set, ki)) {
            if((k2 = kh_get(64, kc, kh_key(set, ki))) == kh_end(kc)) {
                k2 = kh_put(64, kc, kh_key(set, ki), &khr);
                kh_val(kc, k2) = TDencode(node_depth(tax, kh_val(kc, ki)), kh_val(kc, ki));
                if(unlikely(kh_size(kc) % 1000000 == 0)) LOG_INFO("Final hash size %zu\n", kh_size(kc));
            } else if(kh_val(kc, k2) != taxid) {
                do val = lca(tax, taxid, kh_val(kc, k2));
                while(!kh_try_set(64, kc, k2, val == (tax_t)-1 ? 1: TDencode(node_depth(tax, val), val)));
            }
        }
    }
    LOG_DEBUG("After updating with set of size %zu, total set current size is %zu.\n", kh_size(set), kh_size(kc));
}
inline void update_feature_counter(khash_t(64) *kc, const khash_t(all) *set, const khash_t(p) *tax, const tax_t taxid) {
    // TODO: make this threadsafe.
    int khr;
    khint_t k2;
    for(khiter_t ki(kh_begin(set)); ki != kh_end(set); ++ki) {
        if(kh_exist(set, ki)) {
           if((k2 = kh_get(64, kc, kh_key(set, ki))) == kh_end(kc)) {
                k2 = kh_put(64, kc, kh_key(set, ki), &khr);
                kh_val(kc, k2) = FMencode(1, node_depth(tax, taxid));
            } else while(!kh_try_set(64, kc, k2, FMencode(FMcount(kh_val(kc, k2)), lca(tax, taxid, kh_val(kc, k2)))));
        }
    }
}

inline void update_minimized_map(const khash_t(all) *set, const khash_t(64) *full_map, khash_t(c) *ret) {
    khiter_t kif;
    LOG_DEBUG("Size of set: %zu\n", kh_size(set));
    for(khiter_t ki(0); ki != kh_end(set); ++ki) {
        if(!kh_exist(set, ki) || kh_get(c, ret, kh_key(set, ki)) != kh_end(ret))
            continue;
            // If the key is already in the main map, what's the problem?
        if(unlikely((kif = kh_get(64, full_map, kh_key(set, ki))) == kh_end(full_map)))
            LOG_EXIT("Missing kmer from database... Check for matching spacer and kmer size.\n");
        kh_set(c, ret, kh_key(full_map, kif), kh_val(full_map, kif));
        if(unlikely(kh_size(ret) % 1000000 == 0)) LOG_INFO("Final hash size %zu\n", kh_size(ret));
    }
    return;
}

inline void lca2depth(khash_t(c) *lca_map, const khash_t(p) *tax_map) {
    for(khiter_t ki(kh_begin(lca_map)); ki != kh_end(lca_map); ++ki)
        if(kh_exist(lca_map, ki))
            kh_val(lca_map, ki) = node_depth(tax_map, kh_val(lca_map, ki));
}

inline khash_t(c) *make_depth_hash(khash_t(c) *lca_map, const khash_t(p) *tax_map) {
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


inline khash_t(64) *make_taxdepth_hash(khash_t(c) *kc, const khash_t(p) *tax) {
    khash_t(64) *ret(kh_init(64));
    int khr;
    khiter_t kir;
    kh_resize(64, ret, kc->n_buckets);
    for(khiter_t ki(0); ki != kh_end(kc); ++ki) {
        if(kh_exist(kc, ki)) {
            kir = kh_put(64, ret, kh_key(kc, ki), &khr);
            kh_val(ret, kir) = TDencode(node_depth(tax, kh_val(kc, ki)), kh_val(kc, ki));
        }
    }
    return ret;
}

} // namespace bns
