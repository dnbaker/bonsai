#ifndef _EMP_ENCODER_H__
#define _EMP_ENCODER_H__
#include <thread>
#include <future>
#include <limits>

#include <unistd.h>
#include "klib/kstring.h"
#include "hash.h"
#include "hll/hll.h"
#include "hll/filterhll.h"
#include "entropy.h"
#include "kseq_declare.h"
#include "qmap.h"
#include "spacer.h"
#include "util.h"
#include "klib/kthread.h"
#include <mutex>


namespace bns {
using namespace sketch;


#if USE_HASH_FILLER
namespace detail {
template<typename SketchType, size_t UNROLL_COUNT=4>
struct HashFiller {
    using UType = typename vec::SIMDTypes<u64>::VType;
    static constexpr size_t COUNT = vec::SIMDTypes<u64>::COUNT;
    SketchType &ref_;
    std::array<UType, UNROLL_COUNT> vec_;
    unsigned int count_;
    INLINE HashFiller(SketchType &ref): ref_(ref), count_(0) {}
    INLINE void add(u64 val) {
        vec_[count_ / COUNT].arr_[count_ % COUNT] = val;
        if(++count_ == COUNT * UNROLL_COUNT) {
            for(const auto el: vec_)
                ref_.addh(el.simd_);
            count_ = 0;
        }
    }
    INLINE ~HashFiller() {
        for(;count_;ref_.addh(vec_[count_ / COUNT].arr_[(count_ - 1) % COUNT]), --count_);
    }
};
} // namespace detail
#endif

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

static INLINE u64 lex_score(u64 i, UNUSED(void *data)) {return i ^ XOR_MASK;}
static INLINE u64 ent_score(u64 i, void *data) {
    // For this, the highest-entropy kmers will be selected as "minimizers".
    return UINT64_C(-1) - static_cast<u64>(UINT64_C(7958933093282078720) * kmer_entropy(i, *(unsigned *)data));
}
static INLINE u64 hash_score(u64 i, void *data) {
    khint_t k1;
    khash_t(64) *hash((khash_t(64) *)data);
    if(likely((k1 = kh_get(64, hash, i)) == kh_end(hash))) return kh_val(hash, k1);
    for(k1 = 0; k1 != kh_end(hash); ++k1) {
        LOG_DEBUG("Did not find key. Scanning.\n");
        if(kh_key(hash, k1) == i) __ac_set_isdel_false(hash->flags, k1);
        return kh_val(hash, k1);
    }
    std::fprintf(stderr, "i: %" PRIu64 "\n", i);
    std::exit(EXIT_FAILURE);
    __builtin_unreachable();
    return 0uL;
}

namespace score {
#define DECHASH(name, fn) struct name {u64 operator()(u64 i, void *data) const {return fn(i, data);}}
DECHASH(Lex, lex_score);
DECHASH(Entropy, ent_score);
DECHASH(Hash, hash_score);
#undef DECHASH
} // namespace score



static std::array<u64, 256> make_nthash_lut(u64 seedseed) {
    aes::AesCtr<uint64_t> gen(seedseed);
    uint64_t a = gen(), c = gen(), g = gen(), t = gen();
    std::array<u64, 256> ret;
    std::fill(ret.begin(), ret.end(), 0);
    ret[4] = ret['a'] = ret['A'] = a;
    ret[7] = ret['c'] = ret['C'] = c;
    ret[3] = ret['g'] = ret['G'] = g;
    ret[1] = ret['t'] = ret['T'] = t;
    return ret;
}

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
    std::array<u64, 256> *rolling_seeds_;
    qmap_t     qmap_; // queue of max scores and std::map which keeps kmers, scores, and counts so that we can select the top kmer for a window.
    const ScoreType  scorer_; // scoring struct
    bool canonicalize_;
#if 0
template<unsigned b1, unsigned b2, unsigned range=1>
INLINE uint64_t swapbits(uint64_t x) {
template<size_t lbits>
INLINE uint64_t rol(uint64_t x, unsigned s) {
template<size_t n>
INLINE uint64_t lrot(uint64_t x) {
    return (x << n) | (x >> (64 - n));
}
template<size_t n>
INLINE uint64_t rrot(uint64_t x) {
    return (x >> n) | (x << (64 - n));
}
#endif
static constexpr uint64_t DEFAULT_ROLLING_SEED = UINT64_C(0xB0BAFE77C001D00D);

public:
    Encoder(char *s, u64 l, const Spacer &sp, void *data=nullptr,
            bool canonicalize=true, uint64_t rolling_seed=0):
      s_(s),
      l_(l),
      sp_(sp),
      pos_(0),
      data_(data),
      rolling_seeds_(nullptr),
      qmap_(sp_.w_ - sp_.c_ + 1),
      scorer_{},
      canonicalize_(canonicalize)
    {
        LOG_DEBUG("Canonicalizing: %s\n", canonicalize_ ? "True": "False");
        if(std::is_same_v<ScoreType, score::Entropy> && sp_.unspaced() && !sp_.unwindowed()) {
            if(data_) RUNTIME_ERROR("No data pointer must be provided for lex::Entropy minimization.");
            data_ = static_cast<void *>(new CircusEnt(sp_.k_));
        }
        if(rolling_seed) {
            if((rolling_seeds_ = static_cast<std::array<u64, 256> *>(std::malloc(sizeof(*rolling_seeds_)))) == nullptr) throw std::bad_alloc();
            *rolling_seeds_ = make_nthash_lut(rolling_seed);
        }
    }
    Encoder(const Spacer &sp, void *data, bool canonicalize=true, uint64_t rolling_seed=0): Encoder(nullptr, 0, sp, data, canonicalize, rolling_seed) {}
    Encoder(const Spacer &sp, bool canonicalize=true, uint64_t rolling_seed=0): Encoder(sp, nullptr, canonicalize, rolling_seed) {}
    Encoder(const Encoder &other): Encoder(other.sp_, other.data_) {
        canonicalize_ = other.canonicalize_;
        if(other.rolling_seeds_) {
            if((rolling_seeds_ = static_cast<std::array<u64, 256> *>(std::malloc(sizeof(*rolling_seeds_)))) == nullptr) throw std::bad_alloc();
            std::memcpy(rolling_seeds_, other.rolling_seeds_, sizeof(std::array<u64, 256>));
        }
    }
    Encoder(unsigned k, bool canonicalize=true, uint64_t rolling_seed=0): Encoder(nullptr, 0, Spacer(k), nullptr, canonicalize, rolling_seed) {}

    // Assign functions: These tell the encoder to fetch kmers from this string.
    // kstring and kseq are overloads which call assign(char *s, u64 l) on
    // the correct portions of the structs.
    INLINE void assign(const char *s, u64 l) {
        s_ = s; l_ = l; pos_ = 0;
        if(!sp_.unwindowed())
            qmap_.reset();
        assert((l_ >= sp_.c_ || (!has_next_kmer())) || std::fprintf(stderr, "l: %zu. c: %zu. pos: %zu\n", size_t(l), size_t(sp_.c_), size_t(pos_)) == 0);
    }
    INLINE void assign(kstring_t *ks) {assign(ks->s, ks->l);}
    INLINE void assign(kseq_t    *ks) {assign(&ks->seq);}
#if THEY_SEE_ME_ROLLIN
    template<typename Functor>
    INLINE void for_each_rolling_hash(const Functor &func, u64 seed=0xB0BAF377C001D00D) {
        // Uses nthash, a rolling hash for nucleotides.
        if(unlikely(!sp_.unspaced())) throw std::runtime_error("Can't perform "s + __PRETTY_FUNCTION__ + " for spaced seeds. Use for_each and hash it yourself.");
        if((!sp_.unwindowed())) throw std::runtime_error("Can't perform "s + __PRETTY_FUNCTION__ + " for windowed seeds. Use for_each and hash it yourself.");
        const u64 mask((UINT64_C(-1)) >> (64 - (sp_.k_ << 1)));
        u64 min;
        unsigned filled;
        loop_start:
        min = filled = 0;
        // TODO: this
#if 0
        while(likely(pos_ < l_)) {
            while(filled < sp_.k_ && likely(pos_ < l_)) {
                min <<= 2;
                //std::fprintf(stderr, "Encoding character %c with value %u at last position.\n", s_[pos_], (unsigned)cstr_lut[s_[pos_]]);
                if(unlikely((min |= cstr_lut[s_[pos_++]]) == BF) && (sp_.k_ < 31 || cstr_lut[s_[pos_ - 1]] != 'T')) {
                    goto loop_start;
                }
                ++filled;
            }
            if(likely(filled == sp_.k_)) {
                min &= mask;
                func(min);
                --filled;
            }
        }
#endif
    }
    template<typename Functor>
    INLINE void for_each_rolling_hash_canon(const Functor &func) {
        // TODO: this
    }
    template<typename Functor, size_t N>
    INLINE void for_each_rolling_hash_seed_multiply(const Functor &func) {
        // TODO: this (uses a stack of values with seeds to perform the hashing)
    }
    template<typename Functor>
    INLINE void for_each_rolling_hash_seed_multiply_n(const Functor &func, size_t n) {

        // TODO: this (uses a heap-allocated array of seeds to perform the hashing)
    }
#endif /* THEY_SEE_ME_ROLLIN */

    template<typename Functor>
    INLINE void for_each_canon_windowed(const Functor &func) {
        u64 min;
        while(likely(has_next_kmer()))
            if((min = next_canonicalized_minimizer()) != BF)
                func(min);
    }
    template<typename Functor>
    INLINE void for_each_canon_unwindowed(const Functor &func) {
        if(sp_.unspaced())
            for_each_uncanon_unspaced_unwindowed([&](u64 min) {return func(canonical_representation(min, sp_.k_));});
        else {
            u64 min;
            while(likely(has_next_kmer()))
                if((min = next_kmer()) != BF)
                    func(canonical_representation(min, sp_.k_));
        }
    }
    template<typename Functor>
    INLINE void for_each_uncanon_spaced(const Functor &func) {
        u64 min;
        while(likely(has_next_kmer()))
            if((min = next_minimizer()) != BF)
                func(min);
    }
    template<typename Functor>
    INLINE void for_each_uncanon_unspaced_unwindowed(const Functor &func) {
        const u64 mask((UINT64_C(-1)) >> (64 - (sp_.k_ << 1)));
        u64 min;
        unsigned filled;
        loop_start:
        min = filled = 0;
        while(likely(pos_ < l_)) {
            while(filled < sp_.k_ && likely(pos_ < l_)) {
                min <<= 2;
                if(unlikely((min |= cstr_lut[s_[pos_++]]) == BF) && (sp_.k_ < 32 || cstr_lut[s_[pos_ - 1]] != 'T')) {
                    goto loop_start;
                }
                ++filled;
            }
            if(likely(filled == sp_.k_)) {
                min &= mask;
                func(min);
                --filled;
            }
        }
    }
    template<typename Functor>
    INLINE void for_each_uncanon_unspaced_windowed(const Functor &func) {
        const u64 mask((UINT64_C(-1)) >> (64 - (sp_.k_ << 1)));
        u64 min, kmer;
        unsigned filled;
        windowed_loop_start:
        min = filled = 0;
        while(likely(pos_ < l_)) {
            while(filled < sp_.k_ && likely(pos_ < l_)) {
                min <<= 2;
                if(unlikely((min |= cstr_lut[s_[pos_++]]) == BF) && likely(sp_.k_ < 32 || cstr_lut[s_[pos_ - 1]] != 'T')) {
                    goto windowed_loop_start;
                }
                ++filled;
            }
            if(likely(filled == sp_.k_)) {
                min &= mask;
                if((kmer = qmap_.next_value(min, scorer_(min, data_))) != BF) func(kmer);
                --filled;
            }
        }
    }
    template<typename Functor>
    INLINE void for_each_uncanon_unspaced_windowed_entropy_(const Functor &func) {
        // NEVER CALL THIS DIRECTLY.
        // This contains instructions for generating uncanonicalized but windowed entropy-minimized kmers.
        const u64 mask((UINT64_C(-1)) >> (64 - (sp_.k_ << 1)));
        u64 min, kmer;
        unsigned filled;
        CircusEnt &ent = *(static_cast<CircusEnt *>(data_));
        windowed_loop_start:
        filled = min = 0;
        while(likely(pos_ < l_)) {
            while(filled < sp_.k_ && likely(pos_ < l_)) {
                min <<= 2;
                if(unlikely((min |= cstr_lut[s_[pos_]]) == BF) && likely(sp_.k_ < 31 || cstr_lut[s_[pos_]] != 'T')) {
                    ++pos_;
                    goto windowed_loop_start;
                }
                ent.push(s_[pos_]);
                ++pos_;
                ++filled;
            }
            if(likely(filled == sp_.k_)) {
                min &= mask;
                if((kmer = qmap_.next_value(min, ent.score())) != BF) func(kmer);
                --filled;
            }
        }
    }
    template<typename Functor>
    INLINE void for_each_canon_unspaced_windowed_entropy_(const Functor &func) {
        this->for_each_uncanon_unspaced_windowed_entropy_([&](u64 &min) {return func(canonical_representation(min, sp_.k_));});
    }
    // Utility 'for-each'-like functions.
    template<typename Functor>
    INLINE void for_each(const Functor &func, const char *str, u64 l) {
        this->assign(str, l);
        if(!has_next_kmer()) return;
        if(canonicalize_) {
            if(sp_.unwindowed()) {
                 for_each_canon_unwindowed(func);
            } else {
                if constexpr(std::is_same_v<ScoreType, score::Entropy>) {
                    if(sp_.unspaced()) for_each_canon_unspaced_windowed_entropy_(func);
                    else               for_each_canon_windowed(func);
                } else for_each_canon_windowed(func);
            }
        } else {
            // Note that an entropy-based score calculation can be sped up for this case.
            // This will benefit from either a special function or an if constexpr
            if(sp_.unspaced()) {
                if(sp_.unwindowed()) for_each_uncanon_unspaced_unwindowed(func);
                else {
                    if constexpr(std::is_same_v<ScoreType, score::Entropy>)
                        for_each_uncanon_unspaced_windowed_entropy_(func);
                    else for_each_uncanon_unspaced_windowed(func);
                }
            } else for_each_uncanon_spaced(func);
        }
    }
    template<typename Functor>
    INLINE void for_each(const Functor &func, kseq_t *ks) {
        while(kseq_read(ks) >= 0) assign(ks), for_each<Functor>(func, ks->seq.s, ks->seq.l);
    }
    template<typename Functor>
    INLINE void for_each_canon(const Functor &func, kseq_t *ks) {
        if(sp_.unwindowed())
            while(kseq_read(ks) >= 0) assign(ks), for_each_canon_unwindowed<Functor>(func);
        else
            while(kseq_read(ks) >= 0) assign(ks), for_each_canon_windowed<Functor>(func);
    }
    template<typename Functor>
    INLINE void for_each_uncanon(const Functor &func, kseq_t *ks) {
        if(sp_.unspaced()) {
            if(sp_.unwindowed()) while(kseq_read(ks) >= 0) assign(ks), for_each_uncanon_unspaced_unwindowed(func);
            else                 while(kseq_read(ks) >= 0) assign(ks), for_each_uncanon_unspaced_windowed(func);
        } else while(kseq_read(ks) >= 0) assign(ks), for_each_uncanon_spaced(func);
    }
    template<typename Functor>
    void for_each_canon(const Functor &func, gzFile fp, kseq_t *ks=nullptr) {
        bool destroy;
        if(ks == nullptr) ks = kseq_init(fp), destroy = true;
        else            kseq_assign(ks, fp), destroy = false;
        //LOG_DEBUG("Destroy is %s. I have just assigned kseq's f->f field to the pointer %p\n", destroy ? "true": "false", fp);
        for_each_canon<Functor>(func, ks);
        if(destroy) kseq_destroy(ks);
    }
    template<typename Functor>
    void for_each_uncanon(const Functor &func, gzFile fp, kseq_t *ks=nullptr) {
        bool destroy;
        if(ks == nullptr) ks = kseq_init(fp), destroy = true;
        else            kseq_assign(ks, fp), destroy = false;
        for_each_uncanon<Functor>(func, ks);
        if(destroy) kseq_destroy(ks);
    }
    template<typename Functor>
    void for_each_canon(const Functor &func, const char *path, kseq_t *ks=nullptr) {
        gzFile fp(gzopen(path, "rb"));
        if(!fp) RUNTIME_ERROR(ks::sprintf("Could not open file at %s. Abort!\n", path).data());
        for_each_canon<Functor>(func, fp, ks);
        gzclose(fp);
    }
    template<typename Functor>
    void for_each_uncanon(const Functor &func, const char *path, kseq_t *ks=nullptr) {
        gzFile fp(gzopen(path, "rb"));
        if(!fp) RUNTIME_ERROR(ks::sprintf("Could not open file at %s. Abort!\n", path).data());
        for_each_uncanon<Functor>(func, fp, ks);
        gzclose(fp);
    }
    template<typename Functor>
    void for_each(const Functor &func, gzFile fp, kseq_t *ks=nullptr) {
        bool destroy;
        if(ks == nullptr) ks = kseq_init(fp), destroy = true;
        else            kseq_assign(ks, fp), destroy = false;
        if(canonicalize_) for_each_canon<Functor>(func, ks);
        else              for_each_uncanon<Functor>(func, ks);
        if(destroy) kseq_destroy(ks);
    }
    template<typename Functor>
    void for_each(const Functor &func, const char *path, kseq_t *ks=nullptr) {
        gzFile fp(gzopen(path, "rb"));
        if(!fp) RUNTIME_ERROR(ks::sprintf("Could not open file at %s. Abort!\n", path).data());
        if(canonicalize_) for_each_canon<Functor>(func, fp, ks);
        else              for_each_uncanon<Functor>(func, fp, ks);
        gzclose(fp);
    }
    template<typename Functor, typename ContainerType,
             typename=std::enable_if_t<std::is_same_v<typename ContainerType::value_type::value_type, char> ||
                                       std::is_same_v<std::decay_t<typename ContainerType::value_type>, char *>
                                      >
            >
    void for_each(const Functor &func, const ContainerType &strcon, kseq_t *ks=nullptr) {
        for(const auto &el: strcon) {
            LOG_DEBUG("Loading from file %s\n", get_cstr(el));
            for_each<Functor>(func, get_cstr(el), ks);
        }
    }

    template<typename Target>
    void add(hll::hll_t &hll, const Target &target, kseq_t *ks=nullptr) {
        this->for_each([&](u64 min) {hll.addh(min);}, target, ks);
    }

    template<typename Target, typename HashType>
    void add(fhll::fhllbase_t<HashType> &hll, const Target &target, kseq_t *ks=nullptr) {
        this->for_each([&](u64 min) {hll.addh(min);}, target, ks);
    }

    template<typename Target>
    void add(khash_t(all) *set, const Target &target, kseq_t *ks=nullptr) {
        int khr;
        this->for_each([&] (u64 min) {
            kh_put(all, set, min, &khr);
            if(unlikely(khr < 0)) throw std::runtime_error(ks::sprintf("Failed to insert key %" PRIu64 " into hash map. Size of map: %zu\n", min, kh_size(set)).data());
        }, target, ks);
    }

    template<typename ContainerType>
    void add(ContainerType &con, ContainerType &strcon, kseq_t *ks=nullptr) {
        for(const auto &el: strcon) add(get_cstr(con, get_cstr(el), ks));
    }

    // Encodes a kmer starting at `start` within string `s_`.
    INLINE u64 kmer(unsigned start) {
        assert(start <= l_ - sp_.c_ + 1);
        if(l_ < sp_.c_)    return BF;
        u64 new_kmer(cstr_lut[s_[start]]);
        if(new_kmer == BF) return BF;
        u64 len(sp_.s_.size());
        const u8 *spaces(sp_.s_.data());
#define ITER do {new_kmer <<= 2, start += *spaces++; if((new_kmer |= cstr_lut[s_[start]]) == BF) return new_kmer;} while(0);
        DO_DUFF(len, ITER);
#undef ITER
        return new_kmer;
    }
    // Whether or not an additional kmer is present in the sequence being encoded.
    INLINE int has_next_kmer() const {
        static_assert(std::is_same_v<decltype((std::int64_t)l_ - sp_.c_ + 1), std::int64_t>, "is not same");
        return (pos_ + sp_.c_ - 1) < l_;
    }
    // This fetches our next kmer for our window. It is immediately placed in the qmap_t,
    // which is a tree map containing kmers and scores so we can keep track of the best-scoring
    // kmer in the window.
    INLINE u64 next_kmer() {
        assert(has_next_kmer());
        return kmer(pos_++);
    }
    // This is the actual point of entry for fetching our minimizers.
    // It wraps encoding and scoring a kmer, updates qmap, and returns the minimizer
    // for the next window.
    INLINE u64 next_minimizer() {
        //if(unlikely(!has_next_kmer())) return BF;
        const u64 k(kmer(pos_++)), kscore(scorer_(k, data_));
        return qmap_.next_value(k, kscore);
    }
    INLINE u64 next_canonicalized_minimizer() {
        assert(has_next_kmer());
        const u64 k(canonical_representation(kmer(pos_++), sp_.k_)), kscore(scorer_(k, data_));
        return qmap_.next_value(k, kscore);
    }
    elscore_t max_in_queue() const {
        return qmap_.begin()->first;
    }
    bool canonicalize() const {return canonicalize_;}
    void set_canonicalize(bool value) {canonicalize_ = value;}
    auto pos()   const {return pos_;}
    uint32_t k() const {return sp_.k_;}
    ~Encoder() {
        if(std::is_same_v<ScoreType, score::Entropy> && sp_.unspaced() && !sp_.unwindowed()) {
            delete static_cast<CircusEnt *>(data_);
        }
        std::free(rolling_seeds_);
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
                              bool canonicalize, void *data=nullptr) {

    Encoder<ScoreType> enc(nullptr, 0, space, data, canonicalize);
    khash_t(all) *ret(kh_init(all));
    enc.add(ret, path.data());
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

template<typename ScoreType, typename SketchType>
void fill_lmers(SketchType &sketch, const std::string &path, const Spacer &space, bool canonicalize=true,
                void *data=nullptr, kseq_t *ks=nullptr) {
#if USE_HASH_FILLER
    detail::HashFiller<SketchType> hf(sketch);
#endif
    LOG_DEBUG("Canonicalizing: %s\n", canonicalize ? "true": "false");
    sketch.not_ready();
    Encoder<ScoreType> enc(nullptr, 0, space, data, canonicalize);
#if USE_HASH_FILLER
    enc.for_each([&](u64 min) {hf.add(min);}, path.data(), ks);
#else
    enc.for_each([&](u64 min) {sketch.addh(min);}, path.data(), ks);
#endif
}

template<typename ScoreType>
void hll_fill_lmers(hll::hll_t &hll, const std::string &path, const Spacer &space, bool canonicalize=true,
                    void *data=nullptr, kseq_t *ks=nullptr) {
    fill_lmers<ScoreType>(hll, path, space, canonicalize, data, nullptr);
}

#define SUB_CALL \
    std::async(std::launch::async,\
               hashcount_lmers<ScoreType>, paths[submitted], space, canonicalize, data)

template<typename ScoreType>
u64 count_cardinality(const std::vector<std::string> paths,
                      unsigned k, uint16_t w, spvec_t spaces,
                      bool canonicalize,
                      void *data=nullptr, int num_threads=-1) {
    // Default to using all available threads.
    if(num_threads < 0) num_threads = sysconf(_SC_NPROCESSORS_ONLN);
    const Spacer space(k, w, spaces);
    u64 submitted(0), completed(0), todo(paths.size());
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
    u64 ret(hashes[0]->n_occupied);
    for(auto i: hashes) khash_destroy(i);
    return ret;
}

template<typename SketchType>
inline hll::hll_t &get_hll(SketchType &s);
template<>
inline hll::hll_t &get_hll<hll::hll_t>(hll::hll_t &s) {return s;}
template<>
inline hll::hll_t &get_hll<fhll::fhll_t>(fhll::fhll_t &s) {return s.hll();}
template<>
inline hll::hll_t &get_hll<fhll::pcfhll_t>(fhll::pcfhll_t &s) {return s.hll();}

template<typename SketchType> inline const hll::hll_t &get_hll(const SketchType &s);
template<>
inline const hll::hll_t &get_hll<hll::hll_t>(const hll::hll_t &s) {return s;}
template<>
inline const hll::hll_t &get_hll<fhll::fhll_t>(const fhll::fhll_t &s) {return s.hll();}

template<typename SketchType>
struct est_helper {
    const Spacer                      &sp_;
    const std::vector<std::string> &paths_;
    std::mutex                         &m_;
    const u64                          np_;
    const bool                      canon_;
    void                            *data_;
    std::vector<SketchType>         &hlls_;
    kseq_t                            *ks_;
};

template<typename SketchType, typename ScoreType=score::Lex>
void est_helper_fn(void *data_, long index, int tid) {
    est_helper<SketchType> &h(*(est_helper<SketchType> *)(data_));
    fill_lmers<ScoreType, SketchType>(h.hlls_[tid], h.paths_[index], h.sp_, h.canon_, h.data_, h.ks_ + tid);
}

template<typename SketchType, typename ScoreType=score::Lex>
void fill_sketch(SketchType &ret, const std::vector<std::string> &paths,
              unsigned k, uint16_t w, const spvec_t &spaces, bool canon=true,
              void *data=nullptr, int num_threads=1, u64 np=23, kseq_t *ks=nullptr) {
    // Default to using all available threads if num_threads is negative.
    if(num_threads < 0) {
        num_threads = std::thread::hardware_concurrency();
        LOG_INFO("Number of threads was negative and has been adjusted to all available threads (%i).\n", num_threads);
    }
    const Spacer space(k, w, spaces);
    if(num_threads <= 1) {
        LOG_DEBUG("Starting serial\n");
        for(u64 i(0); i < paths.size(); fill_lmers<ScoreType, SketchType>(ret, paths[i++], space, canon, data, ks));
    } else {
        LOG_DEBUG("Starting parallel\n");
        std::mutex m;
        KSeqBufferHolder kseqs(num_threads);
        std::vector<SketchType> sketches;
        while(sketches.size() < (unsigned)num_threads) sketches.emplace_back(ret.clone());
        if(ret.size() != sketches.back().size()) throw "a party";
        est_helper<SketchType> helper{space, paths, m, np, canon, data, sketches, kseqs.data()};
        {
            ForPool pool(num_threads);
            pool.forpool(&est_helper_fn<SketchType, ScoreType>, &helper, paths.size());
        }
        auto &rhll = get_hll(ret);
        for(auto &sketch: sketches) rhll += get_hll(sketch);
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
                unsigned k, uint16_t w, spvec_t spaces, bool canon=true,
                void *data=nullptr, int num_threads=1, u64 np=23, kseq_t *ks=false, hll::EstimationMethod estim=hll::EstimationMethod::ERTL_MLE, uint16_t jestim=hll::JointEstimationMethod::ERTL_JOINT_MLE, bool clamp=true) {
    hll::hll_t master(np, estim, (hll::JointEstimationMethod)jestim, 1, clamp);
    fill_sketch<hll::hll_t, ScoreType>(master, paths, k, w, spaces, canon, data, num_threads, np, ks);
    return master;
}

template<typename ScoreType=score::Lex>
u64 estimate_cardinality(const std::vector<std::string> &paths,
                            unsigned k, uint16_t w, spvec_t spaces, bool canon,
                            void *data=nullptr, int num_threads=-1, u64 np=23, kseq_t *ks=nullptr, hll::EstimationMethod estim=hll::EstimationMethod::ERTL_MLE) {
    auto tmp(make_hll<ScoreType>(paths, k, w, spaces, canon, data, num_threads, np, ks, estim));
    return tmp.report();
}

} //namespace bns
#endif // _EMP_ENCODER_H__
