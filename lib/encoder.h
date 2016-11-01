#ifndef _ENCODER_H_
#define _ENCODER_H_
#include <thread>
#include <future>

#include <zlib.h>
#include <unistd.h>

#include "htslib/kstring.h"

#include "hash.h"
#include "hll.h"
#include "kseq_declare.h"
#include "qmap.h"
#include "spacer.h"
#include "util.h"


namespace kpg {

enum score_scheme {
    LEX = 0,
    TAX_DEPTH = 1,
    FEATURE_COUNT = 2
};

struct lca_tax_t {
    khash_t(c) *lca;
    khash_t(p) *tax;
};

static INLINE int is_lt(uint64_t i, uint64_t j, UNUSED(void *data)) {
    return i < j;
}

static INLINE int tax_is_lt(uint64_t i, uint64_t j, void *data) {
    lca_tax_t *hashes((lca_tax_t *)data);
    khint_t k1, k2;
    if(unlikely((k1 = kh_get(c, hashes->lca, i)) == kh_end(hashes->lca))) goto fail;
    if(unlikely((k2 = kh_get(c, hashes->lca, j)) == kh_end(hashes->lca))) goto fail;
    return node_depth(hashes->tax, kh_val(hashes->lca, k1)) < node_depth(hashes->tax, kh_val(hashes->lca, k2));
    fail:
        fprintf(stderr, "i: %" PRIu64 ". j: %" PRIu64 ". Failed? %s.\n", i, j,
                i == kh_end(hashes->lca) ? j == kh_end(hashes->lca) ? "both": "i" : "j");
        exit(EXIT_FAILURE);
}

static INLINE int hashval_is_lt(uint64_t i, uint64_t j, void *data) {
    khash_t(c) *h((khash_t(c) *)data);
    khint_t k1, k2;
    if(unlikely((k1 = kh_get(c, h, i)) == kh_end(h))) goto fail;
    if(unlikely((k2 = kh_get(c, h, j)) == kh_end(h))) goto fail;
    return kh_val(h, k1) < kh_val(h, k2);
    fail:
        fprintf(stderr, "i: %" PRIu64 ". j: %" PRIu64 ". Failed? %s.\n", i, j,
                i == kh_end(h) ? j == kh_end(h) ? "both": "i" : "j");
        exit(EXIT_FAILURE);
}

static INLINE uint64_t lex_score(uint64_t i, UNUSED(void *data)) {
    return i;
}

static INLINE uint64_t tax_score(uint64_t i, void *data) {
    lca_tax_t *hashes((lca_tax_t *)data);
    khint_t k1;
    if(unlikely((k1 = kh_get(c, hashes->lca, i)) == kh_end(hashes->lca))) goto fail;
    return node_depth(hashes->tax, kh_val(hashes->lca, k1));
    fail:
        fprintf(stderr, "i: %" PRIu64 "\n", i);
        exit(EXIT_FAILURE);
        return 0uL;
}

static INLINE uint64_t hash_score(uint64_t i, void *data) {
    khash_t(64) *hash((khash_t(64) *)data);
    khint_t k1;
    if(unlikely((k1 = kh_get(64, hash, i)) == kh_end(hash))) goto fail;
    return kh_val(hash, k1);
    fail:
        fprintf(stderr, "i: %" PRIu64 "\n", i);
        exit(EXIT_FAILURE);
        return 0uL;
}


template<uint64_t (*score)(uint64_t, void *)=lex_score>
class Encoder {
    const char *s_;
    size_t l_;
    const Spacer sp_;
    unsigned pos_;
    void *data_;
    qmap_t qmap_;
public:
#if !NDEBUG
    elscore_t max_in_queue() {
        return qmap_.map_.begin()->first;
    }
    elscore_t max_in_queue_manual() {
        elscore_t t1{BF, BF};
        for(auto &i: qmap_.map_) {
            if(i.first < t1) t1 = i.first;
        }
        return t1;
    }
#endif
    Encoder(char *s, size_t l, const Spacer &sp, void *data=nullptr):
      s_(s),
      l_(l),
      sp_(sp),
      pos_(0),
      data_(data),
      qmap_(sp_.w_ - sp_.k_ + 1)
    {
    }
    INLINE void assign(char *s, size_t l) {
        s_ = s; l_ = l; pos_ = 0;
        qmap_.reset();
    }
    INLINE void assign(kstring_t *ks) {assign(ks->s, ks->l);}
    INLINE void assign(kseq_t *ks) {assign(ks->seq.s, ks->seq.l);}
    INLINE uint64_t kmer(unsigned start) {
        assert(start <= l_ - sp_.c_ + 1);
        uint64_t new_kmer(cstr_lut[s_[start]]);
        for(const auto s: sp_.s_) {
            new_kmer <<= 2;
            start += s;
            new_kmer |= cstr_lut[s_[start]];
            if(start > l_) fprintf(stderr, "Gone too far...\n"), exit(1);
        }
        new_kmer = canonical_representation(new_kmer, sp_.k_) ^ XOR_MASK;
        return new_kmer;
    }
    // Algorithmic inefficiencies
    // 1. Not skipping over previousy discovered ambiguous bases
    // 2. (Solved)
    //    Recalculating kmers for positions shared between neighboring windows.
    //    [Addressed 2 using a linked list of elements with scores and
    //     a red-black tree with scores as keys and counts as values.]
    INLINE uint64_t decode(uint64_t kmer) {return kmer ^ XOR_MASK;}
    INLINE int has_next_kmer() {
        return pos_ < l_ - sp_.c_ + 1;
    }
    INLINE uint64_t next_kmer() {
        assert(has_next_kmer());
        return kmer(pos_++);
    }
    INLINE uint64_t next_minimizer() {
        assert(has_next_kmer());
        //fprintf(stderr, "passed assert\n");
        const uint64_t k(kmer(pos_++)), kscore(score(k, data_));
        //fprintf(stderr, "Got k! %s.\n", sp_.to_string(k).data());
        return qmap_.next_value(k, kscore);
    }
};


template<uint64_t (*score)(uint64_t, void *)>
khash_t(all) *hashcount_lmers(const std::string &path, const Spacer &space,
                      void *data=nullptr) {

    Encoder<score> enc(nullptr, 0, space, data);
    gzFile fp(gzopen(path.data(), "rb"));
    kseq_t *ks(kseq_init(fp));
    khash_t(all) *ret(kh_init(all));
    int khr;
    uint64_t min;
    //fprintf(stderr, "About to start loop. gzfp %p, ksfp %p.\n", (void *)fp, (void *)ks);
    while(kseq_read(ks) >= 0) {
        enc.assign(ks);
        //fprintf(stderr, "Assignd to kseq. name: %s.\n", ks->name.s);
        while(enc.has_next_kmer()) {
        //fprintf(stderr, "Has!: %s.\n", ks->name.s);
            if((min = enc.next_minimizer()) != BF) {
                //fprintf(stderr, "kmer encoded: %s.\n", space.to_string(min).data());
                kh_put(all, ret, min, &khr);
            }
        }
    }
    //fprintf(stderr, "Cleaning up!\n");
    kseq_destroy(ks);
    gzclose(fp);
    //fprintf(stderr, "Exiting!\n");
    return ret;
}

template<uint64_t (*score)(uint64_t, void *), size_t np=22>
hll_t hllcount_lmers(const std::string &path, const Spacer &space,
                      void *data=nullptr) {

    Encoder<score> enc(nullptr, 0, space, data);
    gzFile fp(gzopen(path.data(), "rb"));
    kseq_t *ks(kseq_init(fp));
    hll_t ret(np);
    uint64_t min;
    //fprintf(stderr, "About to start loop. gzfp %p, ksfp %p.\n", (void *)fp, (void *)ks);
    while(kseq_read(ks) >= 0) {
        enc.assign(ks);
        //fprintf(stderr, "Assignd to kseq. name: %s.\n", ks->name.s);
        while(enc.has_next_kmer()) {
        //fprintf(stderr, "Has!: %s.\n", ks->name.s);
            if((min = enc.next_minimizer()) != BF) {
                //fprintf(stderr, "kmer encoded: %s.\n", space.to_string(min).data());
                ret.add(u64hash(min));
            }
        }
    }
    //fprintf(stderr, "Cleaning up!\n");
    kseq_destroy(ks);
    gzclose(fp);
    //fprintf(stderr, "Estimated: %lf with %zu.\n", ret.report(), np);
    //fprintf(stderr, "Exiting!\n");
    return ret;
}

template<uint64_t (*score)(uint64_t, void *)=lex_score>
size_t count_cardinality(const std::vector<std::string> paths,
                         unsigned k, uint16_t w, spvec_t spaces,
                         void *data=nullptr, int num_threads=-1) {
    // Default to using all available threads.
    if(num_threads < 0) num_threads = sysconf(_SC_NPROCESSORS_ONLN);
    const Spacer space(k, w, spaces);
    size_t submitted(0), completed(0), todo(paths.size());
    std::vector<std::future<khash_t(all) *>> futures;
    std::vector<khash_t(all) *> hashes;
    // Submit the first set of jobs
    for(size_t i(0); i < (unsigned)num_threads && i < todo; ++i) {
        futures.emplace_back(std::async(
          std::launch::async, hashcount_lmers<score>, paths[i], space, data));
        ++submitted;
    }
    // Daemon -- check the status of currently running jobs, submit new ones when available.
    while(submitted < todo) {
        for(auto &f: futures) {
            if(is_ready(f)) {
                hashes.push_back(f.get());
                f = std::async(
                  std::launch::async, hashcount_lmers<score>, paths[submitted++],
                  space, data);
                ++completed;
            }
        }
    }
    // Get values from the rest of these threads.
    for(auto &f: futures) if(f.valid()) hashes.push_back(f.get());
    // Combine them all for a final count
    for(auto i(hashes.begin() + 1), end = hashes.end(); i != end; ++i) {
        kset_union(hashes[0], *i);
    }
    size_t ret(hashes[0]->n_occupied);
    for(auto i: hashes) kh_destroy(all, i);
    return ret;
}

template<uint64_t (*score)(uint64_t, void *)=lex_score, size_t np=20>
size_t estimate_cardinality(const std::vector<std::string> paths,
                            unsigned k, uint16_t w, spvec_t spaces,
                            void *data=nullptr, int num_threads=-1) {
    // Default to using all available threads.
    if(num_threads < 0) num_threads = sysconf(_SC_NPROCESSORS_ONLN);
    const Spacer space(k, w, spaces);
    size_t submitted(0), completed(0), todo(paths.size());
    std::vector<std::future<hll_t>> futures;
    std::vector<hll_t> hlls;
    // Submit the first set of jobs
    for(size_t i(0); i < (unsigned)num_threads && i < todo; ++i) {
        futures.emplace_back(std::async(
          std::launch::async, hllcount_lmers<score, np>, paths[i], space, data));
        ++submitted;
    }
    LOG_DEBUG("About to start daemon.\n");
    // Daemon -- check the status of currently running jobs, submit new ones when available.
    while(submitted < todo) {
        for(auto f(futures.begin()); f != futures.end(); ++f) {
            if(submitted == todo) break;
            if(is_ready(*f)) {
                hlls.push_back(f->get());
                futures.erase(f);
                futures.emplace_back(std::async(
                  std::launch::async, hllcount_lmers<score, np>, paths[submitted++],
                  space, data));
                ++completed;
                break; // Iterators are invalid. Start again.
            }
        }
    }
    // Get values from the rest of these threads.
    for(auto &f: futures) if(f.valid()) hlls.push_back(f.get());
    // Combine them all for a final count
    // Note: This could be parallelized by dividing into subsets, summing those subsets,
    // and then summing those sums. In practice, this is already obscenely fast.

    for(auto i(hlls.begin() + 1), end = hlls.end(); i != end; ++i) hlls[0] += *i;
    if(hlls.size() != todo) throw "a party!";
    hlls[0].sum();
    return (size_t)hlls[0].report();
}

} //namespace kpg
#endif // _ENCODER_H_
