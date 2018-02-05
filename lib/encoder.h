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
#include "klib/kthread.h"
#include <mutex>


namespace emp {

enum score_scheme {
    LEX = 0,
    ENTROPY,
    TAX_DEPTH,
    FEATURE_COUNT
};

template<typename T>
static INLINE int is_lt(T i, T j, UNUSED(void *data)) {
    return i < j;
}

using ScoringFunction = u64 (*)(u64, void*);

static INLINE u64 lex_score(u64 i, UNUSED(void *data)) {return i;}
static INLINE u64 ent_score(u64 i, void *data) {
    // For this, the highest-entropy kmers will be selected as "minimizers".
    LOG_WARNING("Note: this will likely need a tie breaker.");
    return UINT64_C(-1) - static_cast<u64>(UINT64_C(7958933093282078720) * kmer_entropy(i, *(unsigned *)data));
}
static INLINE u64 hash_score(u64 i, void *data) {
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
        std::fprintf(stderr, "i: %" PRIu64 "\n", i);
        std::exit(EXIT_FAILURE);
        return 0uL;
}

namespace score {
#define DECHASH(name, fn) struct name {u64 operator()(u64 i, void *data) const { return fn(i, data);}}
DECHASH(Lex, lex_score);
DECHASH(Entropy, ent_score);
DECHASH(Hash, hash_score);
#undef DECHASH
} // namespace score



/*
 *Encoder:
 * Uses a Spacer to control spacing.
 * It keeps a sliding window of best-scoring kmers and their scores.
 * To switch between sequences, use the assign functions.
 *
 * BF signals overflow.
 */
template<typename ScoreType=score::Lex>
class Encoder {
    const char   *s_; // String from which we are encoding our kmers.
    u64           l_; // Length of the string
public:
    const Spacer sp_; // Defines window size, spacing, and kmer size.
private:
    u64         pos_; // Current position within the string s_ we're working with.
    void      *data_; // A void pointer for using with scoring. Needed for hash_score.
    qmap_t     qmap_; // queue of max scores and std::map which keeps kmers, scores, and counts so that we can select the top kmer for a window.
    const ScoreType scorer_; // scoring struct

public:
    Encoder(char *s, size_t l, const Spacer &sp, void *data=nullptr):
      s_(s),
      l_(l),
      sp_(sp),
      pos_(0),
      data_(data),
      qmap_(sp_.w_ - sp_.c_ + 1),
      scorer_{} {}
    Encoder(const Spacer &sp, void *data): Encoder(nullptr, 0, sp, data) {}
    Encoder(const Spacer &sp): Encoder(sp, nullptr) {}
    Encoder(const Encoder &other): Encoder(other.sp_, other.data_) {}
    Encoder(unsigned k): Encoder(nullptr, 0, Spacer(k), nullptr) {}

    // Assign functions: These tell the encoder to fetch kmers from this string.
    // kstring and kseq are overloads which call assign(char *s, size_t l) on
    // the correct portions of the structs.
    INLINE void assign(char *s, size_t l) {
        s_ = s; l_ = l; pos_ = 0;
        qmap_.reset();
        assert(l_ >= sp_.c_ || !has_next_kmer());
    }
    INLINE void assign(kstring_t *ks) {assign(ks->s, ks->l);}
    INLINE void assign(kseq_t    *ks) {assign(&ks->seq);}

    // Encodes a kmer starting at `start` within string `s_`.
    INLINE u64 kmer(unsigned start) {
        assert(start <= l_ - sp_.c_ + 1);
        if(l_ < sp_.c_) return BF;
        u64 new_kmer(cstr_lut[s_[start]]);
        for(const auto s: sp_.s_) {
            new_kmer <<= 2;
            start += s;
            new_kmer |= cstr_lut[s_[start]];
        }
        new_kmer = canonical_representation(new_kmer, sp_.k_) ^ XOR_MASK;
        return new_kmer;
    }
    // When we encode kmers, we XOR it with this XOR_MASK for good key dispersion
    // in spite of potential redundacies in the sequence.
    INLINE u64 decode(u64 kmer) const {return kmer ^ XOR_MASK;}
    // Whether or not an additional kmer is present in the sequence being encoded.
    INLINE int has_next_kmer() const {
        return pos_ < l_ - sp_.c_ + 1;
        static_assert(std::is_same_v<decltype((std::int64_t)l_ - sp_.c_ + 1), std::int64_t>, "is not same");
    }
    // This fetches our next kmer for our window. It is immediately placed in the qmap_t,
    // which is a tree map containing kmers and scores so we can keep track of the best-scoring
    // kmer in the window.
    INLINE u64 next_kmer() {
        assert(has_next_kmer());
        return kmer(pos_++);
    }
    INLINE u64 next_unspaced_kmer(u64 last_kmer) {
        assert(has_next_kmer() && sp_.unspaced());
        if(unlikely(last_kmer == BF)) {
            if(likely((pos_ += sp_.k_) < l_ - sp_.c_)) {
                last_kmer = UINT64_C(0);
                for(unsigned i(0); i < sp_.k_; ++i) {
                    last_kmer |= cstr_lut[s_[pos_++]];
                    last_kmer <<= 2;
                    ++i;
                }
            } else {
                pos_ = l_ - sp_.c_ + 1;
                last_kmer = BF; // To signal to not use this kmer.
            }
        } else last_kmer <<= 2, last_kmer |= cstr_lut[s_[pos_++ + sp_.k_]];
        return last_kmer;
    }
    // This is the actual point of entry for fetching our minimizers.
    // It wraps encoding and scoring a kmer, updates qmap, and returns the minimizer
    // for the next window.
    INLINE u64 next_minimizer() {
        assert(has_next_kmer());
        const u64 k(kmer(pos_++)), kscore(scorer_(k, data_));
        return qmap_.next_value(k, kscore);
    }
    elscore_t max_in_queue() const {
        return qmap_.begin()->first;
    }
};

template<typename ScoreType, typename KhashType>
void add_to_khash(KhashType *kh, Encoder<ScoreType> &enc, kseq_t *ks) {
    u64 min(BF);
    int khr;
    if(enc.sp_.unwindowed()) {
        if(enc.sp_.unspaced()) {
            while(kseq_read(ks) >= 0) {
                enc.assign(ks);
                while(enc.has_next_kmer())
                    if((min = enc.next_unspaced_kmer(min)) != BF)
                        khash_put(kh, min, &khr);
            }
        } else {
            while(kseq_read(ks) >= 0) {
                enc.assign(ks);
                while(enc.has_next_kmer())
                    if((min = enc.next_kmer()) != BF)
                        khash_put(kh, min, &khr);
            }
        }
    } else {
        while(kseq_read(ks) >= 0) {
            enc.assign(ks);
            while(enc.has_next_kmer())
                if((min = enc.next_minimizer()) != BF)
                    khash_put(kh, min, &khr);
        }
    }
}


template<typename ScoreType>
khash_t(all) *hashcount_lmers(const std::string &path, const Spacer &space,
                              void *data=nullptr) {

    Encoder<ScoreType> enc(nullptr, 0, space, data);
    gzFile fp(gzopen(path.data(), "rb"));
    kseq_t *ks(kseq_init(fp));
    khash_t(all) *ret(kh_init(all));
    add_to_khash<ScoreType, khash_t(all)>(ret, enc, ks);
    kseq_destroy(ks);
    gzclose(fp);
    return ret;
}

template<typename ScoreType>
void add_to_hll(hll::hll_t &hll, kseq_t *ks, Encoder<ScoreType> &enc) {
    u64 min(BF);
    if(enc.sp_.unwindowed()) {
        if(enc.sp_.unspaced()) {
            while(kseq_read(ks) >= 0) {
                enc.assign(ks);
                while(enc.has_next_kmer())
                    if((min = enc.next_unspaced_kmer(min)) != BF)
                        hll.addh(min);
            }
        } else {
            while(kseq_read(ks) >= 0) {
                enc.assign(ks);
                while(enc.has_next_kmer())
                    if((min = enc.next_kmer()) != BF)
                        hll.addh(min);
            }
        }
    } else {
        while(kseq_read(ks) >= 0) {
            enc.assign(ks);
            while(enc.has_next_kmer())
                if((min = enc.next_minimizer()) != BF)
                    hll.addh(min);
        }
    }
}

template<typename ScoreType>
void hll_fill_lmers(hll::hll_t &hll, const std::string &path, const Spacer &space,
                    void *data=nullptr) {
    hll.not_ready();
    Encoder<ScoreType> enc(nullptr, 0, space, data);
    gzFile fp(gzopen(path.data(), "rb"));
    if(fp == nullptr) LOG_EXIT("Could not open file at %s\n", path.data());
    kseq_t *ks(kseq_init(fp));
    LOG_DEBUG("I have initialized ks: %p, gzfp: %p, and am about to fill lmers.\n", (void *)ks, (void *)fp);
    add_to_hll<ScoreType>(hll, ks, enc);
    kseq_destroy(ks);
    gzclose(fp);
}

#define SUB_CALL \
    std::async(std::launch::async,\
               hashcount_lmers<ScoreType>, paths[submitted], space, data)

template<typename ScoreType>
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
    while(futures.size() < (unsigned)num_threads && futures.size() < todo)
        futures.emplace_back(SUB_CALL), ++submitted;
    // Daemon -- check the status of currently running jobs, submit new ones when available.
    while(submitted < todo) {
        static const int max_retries = 10;
        for(auto &f: futures) {
            if(is_ready(f)) {
                hashes.push_back(f.get());
                int success(0), tries(0);
                while(!success) {
                    try {
                        f = SUB_CALL;
                        ++submitted;
                        ++completed;
                        success = 1;
                    } catch (std::system_error &se) {
                          LOG_DEBUG("System error: resource temporarily available. Retry #%i\n", tries + 1);
                          if(++tries >= max_retries) {LOG_EXIT("Exceeded maximum retries\n"); throw;}
                          sleep(1);
                    }
                }
            }
        }
    }
    // Get values from the rest of these threads.
    for(auto &f: futures) if(f.valid()) hashes.push_back(f.get());
    // Combine them all for a final count
    for(auto i(hashes.begin() + 1), end = hashes.end(); i != end; ++i) kset_union(hashes[0], *i);
    size_t ret(hashes[0]->n_occupied);
    for(auto i: hashes) khash_destroy(i);
    return ret;
}

struct est_helper {
    const Spacer                   &sp_;
    const std::vector<std::string> &paths_;
    std::mutex                     &m_;
    const size_t               np_;
    void                           *data_;
    hll::hll_t                     &master_;
};

template<typename ScoreType=score::Lex>
void est_helper_fn(void *data_, long index, int tid) {
    est_helper &h(*(est_helper *)(data_));
    hll_fill_lmers<ScoreType>(h.master_, h.paths_[index], h.sp_, h.data_);
}

template<typename ScoreType=score::Lex>
void fill_hll(hll::hll_t &ret, const std::vector<std::string> &paths,
              unsigned k, uint16_t w, const spvec_t &spaces,
              void *data=nullptr, int num_threads=1, size_t np=23) {
    // Default to using all available threads if num_threads is negative.
#if 0
    LOG_DEBUG("Filling hll of %zu/%zu size, %zu paths, k%u, w%u, data %p, nt %u, sketch size %zu",
              ret.size(), ret.p(), paths.size(), k, w, data, num_threads, np);
    LOG_DEBUG("Spacer: %s\n", ::emp::str(spaces).data());
#endif
    if(num_threads < 0) {
        num_threads = sysconf(_SC_NPROCESSORS_ONLN);
        LOG_WARNING("Number of threads was negative and has been adjusted to all available threads (%i).\n", num_threads);
    }
    const Spacer space(k, w, spaces);
    if(num_threads <= 1) {
        LOG_DEBUG("Starting serial\n");
        for(size_t i(0); i < paths.size(); hll_fill_lmers<ScoreType>(ret, paths[i++], space, data));
    } else {
        LOG_DEBUG("Starting parallel\n");
        std::mutex m;
        est_helper helper{space, paths, m, np, data, ret};
        kt_for(num_threads, &est_helper_fn<ScoreType>, &helper, paths.size());
    }
    
}

template<typename T>
void hll_from_khash(hll::hll_t &ret, const T *kh, bool clear=true) {
    if(clear) {
        LOG_DEBUG("Clearing hll::hll_t ret with address %p\n", (void *)&ret);
        ret.clear();
    }
    for(khiter_t i(0); i < kh_size(kh); ++i)
        if(kh_exist(kh, i))
            ret.addh(kh->keys[i]);
}

template<typename ScoreType=score::Lex>
hll::hll_t make_hll(const std::vector<std::string> &paths,
                unsigned k, uint16_t w, spvec_t spaces,
                void *data=nullptr, int num_threads=1, size_t np=23) {
    hll::hll_t master(np);
    fill_hll(master, paths, k, w, spaces, data, num_threads, np);
    return master;
}

template<typename ScoreType=score::Lex>
size_t estimate_cardinality(const std::vector<std::string> &paths,
                            unsigned k, uint16_t w, spvec_t spaces,
                            void *data=nullptr, int num_threads=-1, size_t np=23) {
    auto tmp(make_hll<ScoreType>(paths, k, w, spaces, data, num_threads, np));
    return tmp.report();
}

} //namespace emp
#endif // _EMP_ENCODER_H__
