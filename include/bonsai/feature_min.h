#pragma once

#include "encoder.h"
#include "spacer.h"
#include "khash64.h"
#include "util.h"
#include "klib/kthread.h"
#include "kmermap.h"
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
inline void update_frequency_counter(khash_t(64) *kc, const khash_t(64) *counts, const khash_t(p) *tax, tax_t taxid);
inline void update_minimized_map(const khash_t(all) *set, const khash_t(64) *full_map, khash_t(c) *ret);

template<typename ScoreType>
// Presize a freshly-allocated khash from the input's byte size so the fill
// loop rehashes once or twice instead of ~log2(n) times. Compressed inputs
// underestimate the sequence length, which only leaves a couple of growth
// rounds. No-op for non-empty or already-large tables (reused tables).
#define BNS_KHASH_PRESIZE(name, h, path) do { \
    if(kh_size(h) == 0 && (h)->n_buckets < 256) { \
        struct stat sb; \
        if(::stat(path, &sb) == 0 && sb.st_size > 0) { \
            const u64 want(u64(sb.st_size) * 2 + 64); \
            kh_resize(name, h, want > (u64(1) << 26) ? khint_t(1u << 26): khint_t(want)); \
        } \
    } \
} while(0)

size_t fill_set_genome(const char *path, const Spacer &sp, khash_t(all) *ret, size_t index, void *data, bool canon, kseq_t *ks=nullptr);
template<typename ScoreType>
size_t fill_count_genome(const char *path, const Spacer &sp, khash_t(64) *ret, size_t index, void *data, bool canon, kseq_t *ks=nullptr);

// Wrap these in structs so that downstream code can be managed as a set, not updated one-by-one.
struct LcaMap {
    using ReturnType = khash_t(c) *;
    using CounterType = khash_t(all);
    static constexpr size_t ValSize = sizeof(*(ReturnType{0})->vals);
    template<typename ScoreType>
    static size_t fill(const char *path, const Spacer &sp, khash_t(all) *ret, size_t index, void *data, bool canon, kseq_t *ks) {
        return fill_set_genome<ScoreType>(path, sp, ret, index, data, canon, ks);
    }
    static void update(const khash_t(p) *tax, const khash_t(all) *set, const khash_t(64) *, khash_t(c) *r32, khash_t(64) *, tax_t taxid) {
        update_lca_map(r32, set, tax, taxid);
    }
};
// Builds the LCA database straight into a CompactKmerMap (see
// CompactLcaBuilder, kmermap.h) through make_map's scheduling; no khash_t(c)
// is allocated. A "custom sink" updater: make_map hands it each genome's set
// and returns release().
struct CompactLcaMap {
    using ReturnType = std::unique_ptr<CompactKmerMap>;
    using CounterType = khash_t(all);
    static constexpr size_t ValSize = 4;
    static constexpr bool kCustomSink = true;
    std::unique_ptr<CompactLcaBuilder> b;
    template<typename ScoreType>
    static size_t fill(const char *path, const Spacer &sp, khash_t(all) *ret, size_t index, void *data, bool canon, kseq_t *ks) {
        return fill_set_genome<ScoreType>(path, sp, ret, index, data, canon, ks);
    }
    void update(const khash_t(p) *tax, const khash_t(all) *set, const khash_t(64) *, khash_t(c) *, khash_t(64) *, tax_t taxid) {
        b->update(set, tax, taxid);
    }
    ReturnType release() {return b->finish();}
};
template<typename T, typename = void> struct has_custom_sink: std::false_type {};
template<typename T> struct has_custom_sink<T, std::void_t<decltype(T::kCustomSink)>>: std::integral_constant<bool, T::kCustomSink> {};

struct TdMap {
    using ReturnType = khash_t(64) *;
    using CounterType = khash_t(all);
    static constexpr size_t ValSize = sizeof(*(ReturnType{0})->vals);
    template<typename ScoreType>
    static size_t fill(const char *path, const Spacer &sp, khash_t(all) *ret, size_t index, void *data, bool canon, kseq_t *ks) {
        return fill_set_genome<ScoreType>(path, sp, ret, index, data, canon, ks);
    }
    static void update(const khash_t(p) *tax, const khash_t(all) *set, const khash_t(64) *, khash_t(c) *, khash_t(64) *r64, tax_t taxid) {
        update_td_map(r64, set, tax, taxid);
    }
};
struct FcMap {
    using ReturnType = khash_t(64) *;
    using CounterType = khash_t(all);
    static constexpr size_t ValSize = sizeof(*(ReturnType{0})->vals);
    template<typename ScoreType>
    static size_t fill(const char *path, const Spacer &sp, khash_t(all) *ret, size_t index, void *data, bool canon, kseq_t *ks) {
        return fill_set_genome<ScoreType>(path, sp, ret, index, data, canon, ks);
    }
    static void update(const khash_t(p) *tax, const khash_t(all) *set, const khash_t(64) *, khash_t(c) *, khash_t(64) *r64, tax_t taxid) {
        update_feature_counter(r64, set, tax, taxid);
    }
};
struct MinMap {
    using ReturnType = khash_t(c) *;
    using CounterType = khash_t(all);
    static constexpr size_t ValSize = sizeof(*(ReturnType{0})->vals);
    template<typename ScoreType>
    static size_t fill(const char *path, const Spacer &sp, khash_t(all) *ret, size_t index, void *data, bool canon, kseq_t *ks) {
        return fill_set_genome<ScoreType>(path, sp, ret, index, data, canon, ks);
    }
    static void update(const khash_t(p) *, const khash_t(all) *set, const khash_t(64) *d64, khash_t(c) *r32, khash_t(64) *, tax_t) {
        update_minimized_map(set, d64, r32);
    }
};
struct FqMap {
    using ReturnType = khash_t(64) *;
    using CounterType = khash_t(64);
    static constexpr size_t ValSize = sizeof(*(ReturnType{0})->vals);
    template<typename ScoreType>
    static size_t fill(const char *path, const Spacer &sp, khash_t(64) *ret, size_t index, void *data, bool canon, kseq_t *ks) {
        return fill_count_genome<ScoreType>(path, sp, ret, index, data, canon, ks);
    }
    static void update(const khash_t(p) *tax, const khash_t(64) *counts, const khash_t(64) *, khash_t(c) *, khash_t(64) *r64, tax_t taxid) {
        update_frequency_counter(r64, counts, tax, taxid);
    }
};

static inline void free_counter(khash_t(all) &c) {
    std::free(c.flags);
    std::free(c.keys);
}
static inline void free_counter(khash_t(64) &c) {
    std::free(c.flags);
    std::free(c.keys);
    std::free(c.vals);
}
static inline void clear_counter(khash_t(all) *c) { kh_clear(all, c); }
static inline void clear_counter(khash_t(64) *c) { kh_clear(64, c); }

template<typename ScoreType>
size_t fill_set_genome(const char *path, const Spacer &sp, khash_t(all) *ret, size_t index, void *data, bool canon, kseq_t *ks) {
    LOG_ASSERT(ret);
    LOG_DEBUG("Filling from genome at path %s. kseq is pre-allocated ? %s. %p\n", path, ks ? "true": "false", (void *)ks);
    if(!ks) BNS_KHASH_PRESIZE(all, ret, path);

    KmerKeys<ScoreType> enc(sp, canon, data);   // exact k-mers for k <= 32, rolling-hash keys beyond
    enc.for_each([&](u64 x) {
        auto it = kh_get(all, ret, x);
        if(it == kh_end(ret)) {
            int khr;
            kh_put(all, ret, x, &khr);
        }
    }, path, ks);
    LOG_DEBUG("Set of size %lu filled from genome at path %s\n", kh_size(ret), path);
    return index;
}

template<typename Container, typename ScoreType>
size_t fill_set_genome_container(Container &container, const Spacer &sp, khash_t(all) *ret, void *data, bool canon, kseq_t *ks=nullptr) {
    size_t sz(0);
    for(std::string &str: container)
        sz += fill_set_genome<ScoreType>(str.data(), sp, ret, 0, data, canon);
    return sz;
}

// Like fill_set_genome, but counts absolute k-mer frequencies within the genome
// (set-style maps collapse repeated occurrences; this does not).
template<typename ScoreType>
size_t fill_count_genome(const char *path, const Spacer &sp, khash_t(64) *ret, size_t index, void *data, bool canon, kseq_t *ks) {
    LOG_ASSERT(ret);
    LOG_DEBUG("Filling k-mer counts from genome at path %s\n", path);
    if(!ks) BNS_KHASH_PRESIZE(64, ret, path);
    Encoder<ScoreType> enc(0, 0, sp, data, canon);
    enc.for_each([&](auto x) {
        auto it = kh_get(64, ret, x);
        if(it == kh_end(ret)) {
            int khr;
            it = kh_put(64, ret, x, &khr);
            if(unlikely(khr < 0))
                RUNTIME_ERROR(ks::sprintf("Could not insert key %" PRIu64 " to table of size %zu.", x, kh_size(ret)).data());
            kh_val(ret, it) = 1;
        } else ++kh_val(ret, it);
    }, path, ks);
    LOG_DEBUG("Counted %zu distinct k-mers from genome at path %s\n", kh_size(ret), path);
    return index;
}


template<typename ScoreType, typename MapUpdater, typename CounterT=typename MapUpdater::CounterType>
typename MapUpdater::ReturnType
make_map(const std::vector<std::string> fns, const khash_t(p) *tax_map, const char *seq2tax_path, const Spacer &sp, int num_threads, bool canon, size_t start_size, const khash_t(64) *data,
         MapUpdater mu = MapUpdater{}) {
    constexpr bool custom_sink = has_custom_sink<MapUpdater>::value;

    khash_t(c) *r32 = nullptr;
    khash_t(64) *r64 = nullptr;
    // Update this to include tax ids in the hash map.
    size_t submitted(0), completed(0), todo(fns.size());
    if(num_threads < 0) num_threads = std::thread::hardware_concurrency();
    std::vector<CounterT> counters(num_threads);
    std::memset(counters.data(), 0, sizeof(CounterT) * counters.size());
    LOG_DEBUG("Started things\n");
    if constexpr(!custom_sink) {
        if(MapUpdater::ValSize == 8) {
            r64 = static_cast<khash_t(64) *>(std::calloc(1, sizeof(khash_t(64))));
            kh_resize(64, r64, start_size);
        } else {
            r32 = static_cast<khash_t(c) *>(std::calloc(1, sizeof(khash_t(c))));
            kh_resize(c, r32, start_size);
        }
    }
    khash_t(name) *name_hash(build_name_hash(seq2tax_path));
    std::vector<std::future<size_t>> futures;
    // Mkae the future return the kseq pointer and then use it for resubmission.
    // TODO: Also use a fixed st of kh_all sets to reduce memory allocations.
    KSeqBufferHolder kseqs(num_threads);

    // Submit the first set of jobs
    std::set<size_t> used;
    std::vector<uint32_t> genome_counter(todo);
    for(size_t i(0); i < (unsigned)num_threads && i < todo; ++i) {
        futures.emplace_back(std::async(
          std::launch::async, &MapUpdater::template fill<ScoreType>, fns[i].data(), sp, counters.data() + i, i, (void *)data, canon, kseqs.data() + submitted));
        genome_counter[i] = i;
        LOG_DEBUG("Submitted for %zu.\n", submitted);
        ++submitted;
    }

    // Daemon -- check the status of currently running jobs, submit new ones when available.
    while(submitted < todo) {
        auto it = std::find_if(std::begin(futures), std::end(futures), [](auto &f) {return is_ready(f);});
        if(it != futures.end()) {
            auto &f(*it);
            const size_t index(f.get());
            if(submitted == todo) break;
            if(used.find(index) != used.end()) continue;
            used.insert(index);
            const auto coffset = genome_counter[index];
            CounterT *counter = counters.data() + coffset; // Pointer to the counter to use
            const tax_t taxid(get_taxid(fns[index].data(), name_hash));
            mu.update(tax_map, counter, data, r32, r64, taxid);
            clear_counter(counter); // Reset for the next genome assigned to this slot.
            kseq_t *ks_to_submit = kseqs.data() + coffset;
            f = std::async(
              std::launch::async, &MapUpdater::template fill<ScoreType>, fns[submitted].data(),
              sp, counter, submitted, (void *)data, canon, ks_to_submit);
            genome_counter[submitted] = coffset;
            ++submitted, ++completed;
            LOG_DEBUG("Have now submitted %zu element\n", submitted);
        }
    }

    // Join
    for(auto &f: futures) if(f.valid()) {
        const size_t index(f.get());
        const tax_t taxid(get_taxid(fns[index].data(), name_hash));
        mu.update(tax_map, counters.data() + genome_counter[index], data, r32, r64, taxid);
        ++completed;
    }

    // Clean up
    for(auto &counter: counters) free_counter(counter);
    kh_destroy(name, name_hash);
    LOG_DEBUG("Finished making map!\n");
    if constexpr(custom_sink) {
        return mu.release();
    } else {
        if (MapUpdater::ValSize == 8)
            r32 = reinterpret_cast<khash_t(c) *>(r64);
        return reinterpret_cast<typename MapUpdater::ReturnType>(r32);
    }
}

template<typename ScoreType>
auto feature_count_map(const std::vector<std::string> fns, const khash_t(p) *tax_map, const char *seq2tax_path, const Spacer &sp, int num_threads, bool canon, size_t start_size) {
    return make_map<ScoreType, FcMap>(fns, tax_map, seq2tax_path, sp, num_threads, canon, start_size, nullptr);
}

template<typename ScoreType>
khash_t(64) *frequency_map(const std::vector<std::string> fns, const khash_t(p) *tax_map,
                           const char *seq2tax_path, const Spacer &sp, int num_threads, bool canon, size_t start_size) {
    return make_map<ScoreType, FqMap>(fns, tax_map, seq2tax_path, sp, num_threads, canon, start_size, nullptr);
}

template<typename ScoreType>
khash_t(c) *lca_map(const std::vector<std::string> &fns, const khash_t(p) *tax_map,
                    const char *seq2tax_path,
                    const Spacer &sp, int num_threads, bool canon, size_t start_size) {
    return make_map<ScoreType, LcaMap>(fns, tax_map, seq2tax_path, sp, num_threads, canon, start_size, nullptr);
}

// lca_map straight into a compact database (no khash_t(c)): estimated_kmers
// from the cardinality estimate, key_bits from KmerKeys::key_bits_for(k).
template<typename ScoreType>
std::unique_ptr<CompactKmerMap> compact_lca_map(const std::vector<std::string> &fns, const khash_t(p) *tax_map,
                                                const char *seq2tax_path, const Spacer &sp, int num_threads, bool canon,
                                                size_t estimated_kmers, unsigned key_bits) {
    CompactLcaMap mu;
    mu.b = std::make_unique<CompactLcaBuilder>(estimated_kmers, kh_size(tax_map), key_bits);
    return make_map<ScoreType, CompactLcaMap>(fns, tax_map, seq2tax_path, sp, num_threads, canon, 0, nullptr, std::move(mu));
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
    static int warn_missing = 1;
    LOG_DEBUG("Adding set of size %zu to total set of current size %zu.\n", kh_size(set), kh_size(kc));
    for(khiter_t ki(kh_begin(set)); ki < kh_end(set); ++ki) {
        if(kh_exist(set, ki)) {
            if((k2 = kh_get(c, kc, kh_key(set, ki))) == kh_end(kc)) {
                k2 = kh_put(c, kc, kh_key(set, ki), &khr);
                if(unlikely(khr < 0))
                    RUNTIME_ERROR(ks::sprintf("Could not insert key %" PRIu64 " to table of size %zu.", kh_key(set, ki), kh_size(kc)).data());
                kh_val(kc, k2) = taxid;
#if !NDEBUG
                if(unlikely(kh_size(kc) % 1000000 == 0)) LOG_DEBUG("Final hash size %zu\n", kh_size(kc));
#endif
            } else if(kh_val(kc, k2) != taxid) {
                kh_val(kc, k2) = lca(tax, taxid, kh_val(kc, k2));
                if(kh_val(kc, k2) == UINT32_C(1))
                    if(warn_missing) LOG_WARNING("ancestor of %u missing from taxonomy. This is not unexpected considering the issues the NCBI taxonomy has.\n", taxid), warn_missing = 0;
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
    for(khiter_t ki(kh_begin(set)); ki < kh_end(set); ++ki) {
        if(kh_exist(set, ki)) {
            if((k2 = kh_get(64, kc, kh_key(set, ki))) == kh_end(kc)) {
                k2 = kh_put(64, kc, kh_key(set, ki), &khr);
                if(unlikely(khr < 0))
                    RUNTIME_ERROR(ks::sprintf("Could not insert key %" PRIu64 " to table of size %zu.", kh_key(set, ki), kh_size(kc)).data());
                kh_val(kc, k2) = TDencode(node_depth(tax, taxid), taxid);
                if(unlikely(kh_size(kc) % 1000000 == 0)) LOG_INFO("Final hash size %zu\n", kh_size(kc));
            } else if(TDtax(kh_val(kc, k2)) != taxid) {
                do val = lca(tax, taxid, TDtax(kh_val(kc, k2)));
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
    for(khiter_t ki(kh_begin(set)); ki < kh_end(set); ++ki) {
        if(kh_exist(set, ki)) {
           if((k2 = kh_get(64, kc, kh_key(set, ki))) == kh_end(kc)) {
                k2 = kh_put(64, kc, kh_key(set, ki), &khr);
                if(unlikely(khr < 0))
                    RUNTIME_ERROR(ks::sprintf("Could not insert key %" PRIu64 " to table of size %zu.", kh_key(set, ki), kh_size(kc)).data());
                kh_val(kc, k2) = FMencode(1, node_depth(tax, taxid));
            } else while(!kh_try_set(64, kc, k2, FMencode(FMcount(kh_val(kc, k2)), lca(tax, taxid, kh_val(kc, k2)))));
        }
    }
}
// Frequency-based minimization (TODO.md item 5): value tracks absolute k-mer
// frequency summed over genomes (FMcount) alongside the LCA of the genomes
// containing it (FMtax), rather than genome count or taxonomic depth.
inline void update_frequency_counter(khash_t(64) *kc, const khash_t(64) *counts, const khash_t(p) *tax, const tax_t taxid) {
    int khr;
    khint_t k2;
    for(khiter_t ki(kh_begin(counts)); ki < kh_end(counts); ++ki) {
        if(kh_exist(counts, ki)) {
            const u64 count(kh_val(counts, ki));
            if((k2 = kh_get(64, kc, kh_key(counts, ki))) == kh_end(kc)) {
                k2 = kh_put(64, kc, kh_key(counts, ki), &khr);
                if(unlikely(khr < 0))
                    RUNTIME_ERROR(ks::sprintf("Could not insert key %" PRIu64 " to table of size %zu.", kh_key(counts, ki), kh_size(kc)).data());
                kh_val(kc, k2) = FMencode(count, taxid);
            } else {
                const u64 newcount(FMcount(kh_val(kc, k2)) + count);
                kh_val(kc, k2) = FMencode(newcount,
                                          lca(tax, taxid, FMtax(kh_val(kc, k2))));
            }
        }
    }
}

inline void update_minimized_map(const khash_t(all) *set, const khash_t(64) *full_map, khash_t(c) *ret) {
    khiter_t kif;
    LOG_DEBUG("Size of set: %zu\n", kh_size(set));
    for(khiter_t ki(0); ki < kh_end(set); ++ki) {
        if(!kh_exist(set, ki) || kh_get(c, ret, kh_key(set, ki)) != kh_end(ret))
            continue;
            // If the key is already in the main map, what's the problem?
        if(unlikely((kif = kh_get(64, full_map, kh_key(set, ki))) == kh_end(full_map)))
            LOG_EXIT("Missing kmer from database... Check for matching spacer and kmer size.\n");
        if(unlikely(kh_set(c, ret, kh_key(full_map, kif), kh_val(full_map, kif)) < 0))
            RUNTIME_ERROR(ks::sprintf("Failed to update minimized map with kh_set to table of size %zu.", kh_size(ret)).data());
        if(unlikely((kh_size(ret) & 0xFFFFF) == 0)) LOG_INFO("Final hash size %zu\n", kh_size(ret));
    }
    return;
}

inline void lca2depth(khash_t(c) *lca_map, const khash_t(p) *tax_map) {
    for(khiter_t ki(kh_begin(lca_map)); ki < kh_end(lca_map); ++ki)
        if(kh_exist(lca_map, ki))
            kh_val(lca_map, ki) = node_depth(tax_map, kh_val(lca_map, ki));
}

inline khash_t(c) *make_depth_hash(khash_t(c) *lca_map, const khash_t(p) *tax_map) {
    khash_t(c) *ret(kh_init(c));
    kh_resize(c, ret, kh_size(lca_map));
    khiter_t ki1;
    int khr;
    for(khiter_t ki2(kh_begin(lca_map)); ki2 < kh_end(lca_map); ++ki2) {
        if(kh_exist(lca_map, ki2)) {
            ki1 = kh_put(c, ret, kh_key(lca_map, ki2), &khr);
            if(unlikely(khr < 0)) throw std::runtime_error(ks::sprintf("[%s:%d] Failed to insert into to table of size %zu.", __PRETTY_FUNCTION__, __LINE__, kh_size(ret)).data());
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
    for(khiter_t ki(0); ki < kh_end(kc); ++ki) {
        if(kh_exist(kc, ki)) {
            kir = kh_put(64, ret, kh_key(kc, ki), &khr);
            kh_val(ret, kir) = TDencode(node_depth(tax, kh_val(kc, ki)), kh_val(kc, ki));
        }
    }
    return ret;
}

} // namespace bns
