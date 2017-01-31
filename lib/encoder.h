#ifndef _EMP_ENCODER_H__
#define _EMP_ENCODER_H__
#include <thread>
#include <future>
#include <limits>

#include <zlib.h>
#include <unistd.h>

#include "klib/kstring.h"

#include "hash.h"
#include "hll/hll.h"
#include "kseq_declare.h"
#include "qmap.h"
#include "spacer.h"
#include "util.h"


namespace emp {

enum score_scheme {
    LEX = 0,
    TAX_DEPTH = 1,
    FEATURE_COUNT = 2
};

template<typename T>
static INLINE int is_lt(T i, T j, UNUSED(void *data)) {
    return i < j;
}

static INLINE std::uint64_t lex_score(std::uint64_t i, UNUSED(void *data)) {
    return i;
}

static INLINE std::uint64_t hash_score(std::uint64_t i, void *data) {
    khash_t(64) *hash((khash_t(64) *)data);
    khint_t k1;
    if(unlikely((k1 = kh_get(64, hash, i)) == kh_end(hash))) goto fail;
    return kh_val(hash, k1);
    fail:
        for(k1 = 0; k1 != kh_end(hash); ++k1) {
            LOG_DEBUG("Did not find key. Scanning.\n");
            if(kh_key(hash, k1) == i) __ac_set_isdel_false(hash->flags, k1);
            return kh_val(hash, k1);
        }
        fprintf(stderr, "i: %" PRIu64 "\n", i);
        exit(EXIT_FAILURE);
        return 0uL;
}

/*
 *Encoder:
 * Uses a Spacer to control spacing.
 * It keeps a sliding window of best-scoring kmers and their scores.
 * To switch between sequences, use the assign functions.
 *
 * BF signals overflow.
 */
template<std::uint64_t (*score)(std::uint64_t, void *)=lex_score>
class Encoder {
    const char *s_; // String from which we are encoding our kmers.
    std::int64_t l_; // Length of the string
    const Spacer sp_; // Defines window size, spacing, and kmer size.
    int pos_; // Current position within the string s_ we're working with.
    void *data_; // A void pointer for using with scoring. Needed for hash_score.
    qmap_t qmap_; // queue of max scores and std::map which keeps kmers, scores, and counts so that we can select the top kmer for a window.

public:
    // These functions check that the queue is really keeping track of the best kmer for a window.
    elscore_t max_in_queue() {
        return qmap_.begin()->first;
    }
    elscore_t max_in_queue_manual() {
        elscore_t t1{BF, BF};
        for(auto &i: qmap_) {
            if(i.first < t1) t1 = i.first;
        }
        return t1;
    }
    Encoder(char *s, std::size_t l, const Spacer &sp, void *data=nullptr):
      s_(s),
      l_(l),
      sp_(sp),
      pos_(0),
      data_(data),
      qmap_(sp_.w_ - sp_.c_ + 1)
    {
    }
    Encoder(const Spacer &sp, void *data): Encoder(nullptr, 0, sp, data) {
    }
    Encoder(const Spacer &sp): Encoder(sp, nullptr) {
    }
    Encoder(Encoder<score> &other): Encoder(other.sp_, other.data_) {
    }

    // Assign functions: These tell the encoder to fetch kmers from this string.
    // kstring and kseq are overloads which call assign(char *s, std::size_t l) on
    // the correct portions of the structs.
    INLINE void assign(char *s, std::size_t l) {
        s_ = s; l_ = l; pos_ = 0;
        qmap_.reset();
        assert(l_ >= sp_.c_ || !has_next_kmer());
    }
    INLINE void assign(kstring_t *ks) {assign(ks->s, ks->l);}
    INLINE void assign(kseq_t    *ks) {assign(ks->seq.s, ks->seq.l);}

    // Encodes a kmer starting at `start` within string `s_`.
    INLINE std::uint64_t kmer(unsigned start) {
        assert(start <= l_ - sp_.c_ + 1);
        if(l_ < sp_.c_) return BF;
        std::uint64_t new_kmer(cstr_lut[s_[start]]);
        for(const auto s: sp_.s_) {
            new_kmer <<= 2;
            start += s;
            new_kmer |= cstr_lut[s_[start]];
        }
        new_kmer = canonical_representation(new_kmer, sp_.k_) ^ XOR_MASK;
        //LOG_DEBUG("Kmer about to be classified: %s.\n", sp_.to_string(new_kmer).data());
        return new_kmer;
    }
    // When we encode kmers, we XOR it with this XOR_MASK for good key dispersion
    // in spite of potential redundacies in the sequence.
    INLINE std::uint64_t decode(std::uint64_t kmer) {return kmer ^ XOR_MASK;}
    // Whether or not an additional kmer is present in the sequence being encoded.
    INLINE int has_next_kmer() {
        return pos_ < l_ - sp_.c_ + 1;
        static_assert(std::is_same<decltype((std::int64_t)l_ - sp_.c_ + 1), std::int64_t>::value, "is not same");
    }
    // This fetches our next kmer for our window. It is immediately placed in the qmap_t,
    // which is a tree map containing kmers and scores so we can keep track of the best-scoring
    // kmer in the window.
    INLINE std::uint64_t next_kmer() {
        assert(has_next_kmer());
        return kmer(pos_++);
    }
    // This is the actual point of entry for fetching our minimizers.
    // It wraps encoding and scoring a kmer, updates qmap, and returns the minimizer
    // for the next window.
    INLINE std::uint64_t next_minimizer() {
        assert(has_next_kmer());
        const std::uint64_t k(kmer(pos_++)), kscore(score(k, data_));
        return qmap_.next_value(k, kscore);
    }
};


template<std::uint64_t (*score)(std::uint64_t, void *)>
khash_t(all) *hashcount_lmers(const std::string &path, const Spacer &space,
                              void *data=nullptr) {

    Encoder<score> enc(nullptr, 0, space, data);
    gzFile fp(gzopen(path.data(), "rb"));
    kseq_t *ks(kseq_init(fp));
    khash_t(all) *ret(kh_init(all));
    int khr;
    std::uint64_t min;
    while(kseq_read(ks) >= 0) {
        enc.assign(ks);
        while(enc.has_next_kmer()) {
            if((min = enc.next_minimizer()) != BF) {
                kh_put(all, ret, min, &khr);
            }
        }
    }
    kseq_destroy(ks);
    gzclose(fp);
    return ret;
}

template<std::uint64_t (*score)(std::uint64_t, void *)>
hll::hll_t hllcount_lmers(const std::string &path, const Spacer &space,
                          std::size_t np=22, void *data=nullptr) {

    Encoder<score> enc(nullptr, 0, space, data);
    gzFile fp(gzopen(path.data(), "rb"));
    kseq_t *ks(kseq_init(fp));
    hll::hll_t ret(np);
    std::uint64_t min;
    //fprintf(stderr, "About to start loop. gzfp %p, ksfp %p.\n", (void *)fp, (void *)ks);
    while(kseq_read(ks) >= 0) {
        enc.assign(ks);
        //fprintf(stderr, "Assignd to kseq. name: %s.\n", ks->name.s);
        while(enc.has_next_kmer()) {
        //fprintf(stderr, "Has!: %s.\n", ks->name.s);
            if((min = enc.next_minimizer()) != BF) {
                //fprintf(stderr, "kmer encoded: %s.\n", space.to_string(min).data());
                ret.add(wang_hash(min));
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

template<std::uint64_t (*score)(std::uint64_t, void *)=lex_score>
std::size_t count_cardinality(const std::vector<std::string> paths,
                         unsigned k, uint16_t w, spvec_t spaces,
                         void *data=nullptr, int num_threads=-1) {
    // Default to using all available threads.
    if(num_threads < 0) num_threads = sysconf(_SC_NPROCESSORS_ONLN);
    const Spacer space(k, w, spaces);
    std::size_t submitted(0), completed(0), todo(paths.size());
    std::vector<std::future<khash_t(all) *>> futures;
    std::vector<khash_t(all) *> hashes;
    // Submit the first set of jobs
    for(std::size_t i(0); i < (unsigned)num_threads && i < todo; ++i) {
        futures.emplace_back(std::async(
          std::launch::async, hashcount_lmers<score>, paths[i], space, data));
        ++submitted;
    }
    // Daemon -- check the status of currently running jobs, submit new ones when available.
    while(submitted < todo) {
        static const int max_retries = 10;
        for(auto &f: futures) {
            if(is_ready(f)) {
                hashes.push_back(f.get());
                int success(0), tries(0);
                while(!success) {
                    try {
                        f = std::async(
                          std::launch::async, hashcount_lmers<score>, paths[submitted],
                          space, data);
                        ++submitted;
                        ++completed;
                        success = 1;
                    } catch (std::system_error &se) {
                          LOG_DEBUG("System error: resource temporarily available. Retry #%i\n", ++tries);
                          sleep(5);
                          if(tries >= max_retries) {LOG_EXIT("Exceeded maximum retries\n"); throw;}
                    }
                }
            }
        }
    }
    // Get values from the rest of these threads.
    for(auto &f: futures) if(f.valid()) hashes.push_back(f.get());
    // Combine them all for a final count
    for(auto i(hashes.begin() + 1), end = hashes.end(); i != end; ++i) {
        kset_union(hashes[0], *i);
    }
    std::size_t ret(hashes[0]->n_occupied);
    for(auto i: hashes) kh_destroy(all, i);
    return ret;
}

template<std::uint64_t (*score)(std::uint64_t, void *)=lex_score>
std::size_t estimate_cardinality(const std::vector<std::string> &paths,
                                 unsigned k, uint16_t w, spvec_t spaces,
                                 void *data=nullptr, int num_threads=-1, std::size_t np=23, const int high_memory=0) {
    // Default to using all available threads.
    if(num_threads < 0) {
        num_threads = sysconf(_SC_NPROCESSORS_ONLN);
    }
    const Spacer space(k, w, spaces);
    std::size_t submitted(0), completed(0), todo(paths.size());
    std::vector<std::future<hll::hll_t>> futures;
    std::vector<hll::hll_t> hlls;
    // Submit the first set of jobs
    for(std::size_t i(0); i < (unsigned)num_threads && i < todo; ++i) {
        futures.emplace_back(std::async(
          std::launch::async, hllcount_lmers<score>, paths[i], space, np, data));
        ++submitted;
    }
    hll::hll_t master(np * (!high_memory));
    LOG_DEBUG("About to start daemon.\n");
    // Daemon -- check the status of currently running jobs, submit new ones when available.
    while(submitted < todo) {
        static const int max_retries = 10;
        for(auto f(futures.begin()); f != futures.end(); ++f) {
            if(submitted == todo) break;
            if(is_ready(*f)) {
                if(high_memory) hlls.push_back(f->get());
                else            master += f->get();
                futures.erase(f);
#if !NDEBUG
                master.sum();
                LOG_DEBUG("Latest estimate: %lf\n", master.report());
#endif
                int success(0), tries(0);
                while(!success) {
                    try {
                    futures.emplace_back(std::async(
                      std::launch::async, hllcount_lmers<score>, paths[submitted],
                      space, np, data));
                    success = 1;
                    ++submitted;
                    ++completed;
                    } catch(std::system_error &se) {
                      LOG_DEBUG("System error: resource temporarily available. Retry #%i\n", ++tries);
                      sleep(5);
                      if(tries >= max_retries) {LOG_EXIT("Exceeded maximum retries\n"); throw;}
                    }
                }
                break; // Iterators are invalid. Start again.
            }
        }
    }
    // Get values from the rest of these threads.
    // Combine them all for a final count
    // Note: This could be parallelized by dividing into subsets, summing those subsets,
    // and then summing those sums. In practice, this is already obscenely fast.
    if(high_memory) {
        for(auto &f: futures) if(f.valid()) hlls.push_back(f.get());
        for(auto i(hlls.begin() + 1), end = hlls.end(); i != end; ++i) hlls[0] += *i;
        if(hlls.size() != todo) throw "a party!";
        hlls[0].sum();
        return (std::size_t)hlls[0].report();
    }
    for(auto f(std::begin(futures)); f != std::end(futures); ++f) {
        if(f->valid()) {
            hll::hll_t tmp(std::move(f->get()));
            master += tmp;
            master.sum();
            LOG_DEBUG("Latest estimate: %lf\n", master.report());
        }
    }
    master.sum();
    return (std::size_t)master.report();
}

} //namespace emp
#endif // _EMP_ENCODER_H__
