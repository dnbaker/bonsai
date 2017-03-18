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
#define TDencode(depth, taxid) (((std::uint64_t)((tax_t)~0 - depth) << 32) | taxid)

// Decode 64-bit hash for feature counting.
// TODO: add building of FeatureMin hash
#define FMtax(key) ((tax_t)key)
#define FMcount(key) (key >> 32)

#define FMencode(count, taxid) (((std::uint64_t)count << 32) | taxid)


namespace emp {


template<std::uint64_t (*score)(std::uint64_t, void *)>
khash_t(64) *feature_count_map(std::vector<std::string> fns, const Spacer &sp, int num_threads=8);

khash_t(c) *make_depth_hash(khash_t(c) *lca_map, khash_t(p) *tax_map);
void lca2depth(khash_t(c) *lca_map, khash_t(p) *tax_map);
template<std::uint64_t (*score)(std::uint64_t, void *)>
int fill_set_seq(kseq_t *ks, const Spacer &sp, khash_t(all) *ret);

void update_lca_map(khash_t(c) *kc, khash_t(all) *set, khash_t(p) *tax, tax_t taxid, std::shared_mutex &m);
void update_td_map(khash_t(64) *kc, khash_t(all) *set, khash_t(p) *tax, tax_t taxid);
khash_t(64) *make_taxdepth_hash(khash_t(c) *kc, khash_t(p) *tax);
void update_feature_counter(khash_t(64) *kc, khash_t(all) *set, khash_t(p) *tax, const tax_t taxid);


// Return value: whether or not additional sequences were present and added.
template<std::uint64_t (*score)(std::uint64_t, void *)>
int fill_set_seq(kseq_t *ks, const Spacer &sp, khash_t(all) *ret) {
    assert(ret);
    Encoder<score> enc(0, 0, sp, nullptr);
    int khr; // khash return value. Unused, really.
    std::uint64_t kmer;
    if(kseq_read(ks) >= 0) {
        enc.assign(ks);
        while(enc.has_next_kmer())
            if((kmer = enc.next_minimizer()) != BF)
                kh_put(all, ret, kmer, &khr);
        return 1;
    } else return 0;
}

template<std::uint64_t (*score)(std::uint64_t, void *)>
std::size_t fill_set_genome(const char *path, const Spacer &sp, khash_t(all) *ret, std::size_t index, void *data) {
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
    std::uint64_t kmer;
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

template<std::uint64_t (*score)(std::uint64_t, void *), class Container>
std::size_t fill_set_genome_container(Container container, const Spacer &sp, khash_t(all) *ret, void *data) {
    std::size_t sz(0);
    for(std::string &str: container)
        sz += fill_set_genome<score>(str.data(), sp, ret, 0, data);
    return sz;
}

template<std::uint64_t (*score)(std::uint64_t, void *)>
khash_t(64) *ftct_map(std::vector<std::string> &fns, khash_t(p) *tax_map,
                      const char *seq2tax_path,
                      const Spacer &sp, int num_threads, std::size_t start_size) {
    return feature_count_map<score>(fns, tax_map, seq2tax_path, sp, num_threads, start_size);
}

void update_minimized_map(khash_t(all) *set, khash_t(64) *full_map, khash_t(c) *ret);


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

template<std::uint64_t (*score)(std::uint64_t, void *)>
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

template<std::uint64_t (*score)(std::uint64_t, void *)>
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

template<std::uint64_t (*score)(std::uint64_t, void *)>
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

template<std::uint64_t (*score)(std::uint64_t, void *)>
void min_for_helper(void *data_, long index, int tid) {
    LOG_DEBUG("Helper with index %ld and tid %i starting\n", index, tid);
    lca_helper &lh(*(lca_helper *)data_);
    khash_t(all) *h(lh.hashes_[tid]);
    khash_t(c) *ret(lh.ret_);
    if(kh_size(h)) kh_clear(all, h);
    fill_set_genome<score>(lh.fns_[index].data(), lh.sp_, h, index, (khash_t(64) *)lh.data_);
    update_minimized_map(h, (khash_t(64) *)lh.data_, ret);
}

template<std::uint64_t (*score)(std::uint64_t, void *)>
khash_t(c) *minimized_map(std::vector<std::string> fns,
                          khash_t(64) *full_map,
                          const Spacer &sp, int num_threads, std::size_t start_size) {
    khash_t(c) *ret(kh_init(c));
    kh_resize(c, ret, start_size);
    std::vector<khash_t(all) *> counters;
    counters.reserve(num_threads);
    for(std::size_t i(0), end(fns.size()); i != end; ++i) counters.emplace_back(kh_init(all));
    lca_helper lh{sp, fns, ret, nullptr, nullptr, counters, (void *)full_map};
    kt_for(num_threads, &min_for_helper<score>, (void *)&lh, fns.size());
    for(auto c: counters) khash_destroy(c);
    LOG_DEBUG("Cleaned up after minimized map building! Size: %zu\n", kh_size(ret));
    return ret;
}


template<std::uint64_t (*score)(std::uint64_t, void *)>
khash_t(64) *taxdepth_map(std::vector<std::string> &fns, khash_t(p) *tax_map,
                          const char *seq2tax_path, const Spacer &sp,
                          int num_threads, std::size_t start_size) {
    khash_t(64) *ret(kh_init(64));
    kh_resize(64, ret, start_size);
    khash_t(name) *name_hash(build_name_hash(seq2tax_path));
    std::vector<khash_t(all) *> counters;
    counters.reserve(num_threads);
    for(std::size_t i(0), end(fns.size()); i != end; ++i) counters.emplace_back(kh_init(all));
    td_helper th{sp, fns, ret, name_hash, tax_map, counters, nullptr};
    kt_for(num_threads, &td_for_helper<score>, (void *)&th, fns.size());
    destroy_name_hash(name_hash);
    for(auto c: counters) khash_destroy(c);
    LOG_DEBUG("Cleaned up after taxdepth map building!\n");
    return ret;
}


template<std::uint64_t (*score)(std::uint64_t, void *)>
khash_t(c) *lca_map(const std::vector<std::string> &fns, khash_t(p) *tax_map,
                    const char *seq2tax_path,
                    const Spacer &sp, int num_threads, std::size_t start_size) {
    khash_t(c) *ret(kh_init(c));
    kh_resize(c, ret, start_size);
    khash_t(name) *name_hash(build_name_hash(seq2tax_path));
    std::vector<khash_t(all) *> counters;
    std::shared_mutex m;
    counters.reserve(num_threads);
    for(std::size_t i(0), end(fns.size()); i != end; ++i) counters.emplace_back(kh_init(all));
    lca_helper lh{sp, fns, ret, name_hash, tax_map, counters, nullptr, &m};
    kt_for(num_threads, &lca_for_helper<score>, (void *)&lh, fns.size());
    destroy_name_hash(name_hash);
    for(auto c: counters) khash_destroy(c);
    LOG_DEBUG("Cleaned up after LCA map building! Size: %zu\n", kh_size(ret));
    return ret;
}


template <std::uint64_t (*score)(std::uint64_t, void *)>
khash_t(64) *feature_count_map(std::vector<std::string> fns, khash_t(p) *tax_map,
                               const char *seq2tax_path,
                               const Spacer &sp, int num_threads, std::size_t start_size) {
    // TODO: make scale better.
    khash_t(64) *ret(kh_init(64));
    kh_resize(64, ret, start_size);
    khash_t(name) *name_hash(build_name_hash(seq2tax_path));
    std::vector<khash_t(all) *> counters;
    counters.reserve(num_threads);
    for(std::size_t i(0), end(fns.size()); i != end; ++i) counters.emplace_back(kh_init(all));
    td_helper th{sp, fns, ret, name_hash, tax_map, counters, nullptr};
    kt_for(num_threads, &fct_for_helper<score>, (void *)&th, fns.size());
    for(auto c: counters) khash_destroy(c);
    destroy_name_hash(name_hash);
    return ret;
}

} // namespace emp
#endif // #ifdef _FEATURE_MIN__
