#ifndef _EMP_ENCODER_H__
#define _EMP_ENCODER_H__
#include <thread>
#include <future>
#include <limits>

#include <unistd.h>
#include "klib/kstring.h"
#include "hash.h"
#include "entropy.h"
#include "kseq_declare.h"
#include "qmap.h"
#include "spacer.h"
#include "klib/kthread.h"
#include <mutex>
#include "rollinghash/rabinkarphash.h"
#include "rollinghash/cyclichash.h"
#include "ntHash/nthash.hpp"
#include "alphabet.h"
#include "rhtraits.h"
#include "sketch/hash.h"
#include "sketch/div.h"
#include "sketch/exception.h"

namespace bns {
using namespace sketch;


// Aliases
using RollingHashType = InputType;
using RollingHashingType = InputType;
using RHT = InputType;

enum score_scheme {
    LEX = 0,
    ENTROPY,
    TAX_DEPTH,
    FEATURE_COUNT
};

template<typename T>
static INLINE int is_lt(T i, T j, void *) {
    return i < j;
}

using ScoringFunction = u64 (*)(u64, void*);
using FRev64 = sketch::hash::CEIFused<CEIXOR<0x533f8c2151b20f97>, CEIMul<0x9a98567ed20c127d>, RotL<31>, CEIXOR<0x691a9d706391077a>>;

static INLINE u128 lex_score(u128 i) {
    return sketch::hash::CEHasher()(i);
}
static INLINE u128 lex_score(u128 i, void *) {
    return lex_score(i);
}
template<typename T>
static INLINE T ent_score(T i, void *data) {
    return i / (reinterpret_cast<CircusEnt *>(data)->value() + 1e-4);
}
static INLINE u64 lex_score(u64 i) {return FRev64()(i);}
static INLINE u64 lex_score(u64 i, void *) {return FRev64()(i);}
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
struct Lex {
    u64 operator()(u64 i, void *data) const {return lex_score(i, data);}
    u128 operator()(u128 i, void *data) const {return lex_score(i, data);}
};
struct Entropy {
    u64 operator()(u64 i, void *data) const {return ent_score(i, data);}
    u128 operator()(u128 i, void *data) const {return ent_score(i, data);}
};
struct Hash {
    u64 operator()(u64 i, void *data) const {return hash_score(i, data);}
    u128 operator()(u128 i, void *data) const {return hash_score(i, data);}
};
} // namespace score



static std::array<u64, 256> make_nthash_lut(u64 seedseed) {
    RNGType gen(seedseed);
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
 * ENCODE_OVERFLOW signals overflow.
 */
template<typename ScoreType=score::Lex, typename KmerT=uint64_t>
class Encoder {
    const char   *s_; // String from which we are encoding our kmers.
    u64           l_; // Length of the string
public:
    Spacer sp_; // Defines window size, spacing, and kmer size.
    static constexpr KmerT ENCODE_OVERFLOW = static_cast<KmerT>(-1);
private:
    u64         pos_; // Current position within the string s_ we're working with.
    void      *data_ = nullptr; // A void pointer for using with scoring. Needed for hash_score.
    QueueMap<KmerT, KmerT> qmap_; // queue of max scores and std::map which keeps kmers, scores, and counts so that we can select the top kmer for a window.
    const ScoreType  scorer_; // scoring struct
    bool canonicalize_;
    InputType rht = InputType::DNA;
    const int8_t *lutptr = (const int8_t *)DNA4.data();
    size_t nremper = sizeof(KmerT) * 4;
    std::unique_ptr<CircusEnt> ent_tracker_;
    static_assert(std::is_unsigned<KmerT>::value || std::is_same<KmerT, u128>::value, "Must be unsigned integers");

public:

    Encoder(char *s, u64 l, const Spacer &sp, void *data=nullptr,
            bool canonicalize=true):
      s_(s),
      l_(l),
      sp_(sp),
      pos_(0),
      data_(data),
      qmap_(sp_.w_ - sp_.c_ + 1),
      scorer_{},
      canonicalize_(canonicalize)
    {
        if(std::is_same<ScoreType, score::Entropy>::value) {
            ent_tracker_.reset(new CircusEnt(sp_.k_));
        }
        if(!sp_.unspaced() && canonicalize_) {
            std::fprintf(stderr, "If a spaced seed is set, k-mers cannot be canonicalized\n");
            canonicalize_ = false;
        }
    }
    Encoder(const Spacer &sp, void *data, bool canonicalize=true): Encoder(nullptr, 0, sp, data, canonicalize) {}
    Encoder(const Spacer &sp, bool canonicalize=true): Encoder(sp, nullptr, canonicalize) {}
    Encoder(const Encoder &o): s_(o.s_), l_(o.l_), sp_(o.sp_), pos_(o.pos_), data_(o.data_), scorer_(o.scorer_), canonicalize_(o.canonicalize_), rht(o.rht) {
        if(sp_.w_ > sp_.c_)
            qmap_.resize(sp_.w_ - sp_.c_ + 1);
    }
    Encoder(Encoder<ScoreType, KmerT> &&o): s_(o.s_), l_(o.l_), sp_(o.sp_), pos_(o.pos_), data_(o.data_),
            qmap_(std::move(o.qmap_)), scorer_{}, canonicalize_(o.canonicalize_), rht(o.rht) {
        if(o.ent_tracker_) ent_tracker_.reset(new CircusEnt(*o.ent_tracker_));
    }
    Encoder &operator=(const Encoder<ScoreType, KmerT> &o) {
        s_ = o.s_; l_ = o.l_;
        sp_ = o.sp_;
        pos_ = o.pos_;
        data_ = o.data_;
        qmap_ = o.qmap_;
        canonicalize_ = o.canonicalize_;
        rht = o.rht;
        if(o.ent_tracker_) ent_tracker_.reset(new CircusEnt(std::move(*o.ent_tracker_)));
        return *this;
    }
    void hashtype(RollingHashType newrht) {
        rht = newrht; lutptr = rh2lp(rht);
        nremper = rh2n(rht, sizeof(KmerT));
    }
    RollingHashType hashtype() const {return rht;}
    size_t nremperres() const {return nremper;}
    size_t nremperres64() const {return rh2n(rht, 8);}
    size_t nremperres128() const {return rh2n(rht, 16);}
    Encoder(unsigned k, bool canonicalize=true): Encoder(nullptr, 0, Spacer(k), nullptr, canonicalize) {}
    Encoder<score::Entropy, u128> to_entmin128() const {
        Encoder<score::Entropy, u128> ret(sp_, data_, canonicalize_);
        ret.hashtype(this->rht);
        return ret;
    }
    Encoder<score::Entropy, u64> to_entmin64() const {
        Encoder<score::Entropy, u64> ret(sp_, data_, canonicalize_);
        ret.hashtype(this->rht);
        return ret;
    }
    Encoder<ScoreType, u128> to_u128() const {
        Encoder<ScoreType, u128> ret(sp_, data_, canonicalize_);
        ret.hashtype(this->rht);
        return ret;
    }

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


    template<typename Functor>
    INLINE void for_each_canon_windowed(const Functor &func) {
        KmerT min;
        while(likely(has_next_kmer()))
            if((min = next_canonicalized_minimizer()) != ENCODE_OVERFLOW)
                func(min);
    }
    template<typename Functor>
    INLINE void for_each_canon_unwindowed(const Functor &func) {
        auto ffunc = [k=sp_.k_,rht=this->rht,&func](KmerT min) {
            if(rht == DNA) min = canonical_representation(min, k);
            return func(min);
        };
        if(sp_.unspaced()) {
            for_each_uncanon_unspaced_unwindowed(ffunc);
        } else {
            KmerT min;
            while(likely(has_next_kmer()))
                if((min = next_kmer()) != ENCODE_OVERFLOW)
                    ffunc(min);
        }
    }
    template<typename Functor>
    INLINE void for_each_uncanon_spaced(const Functor &func) {
        KmerT min;
        while(likely(has_next_kmer()))
            if((min = next_minimizer()) != ENCODE_OVERFLOW)
                func(min);
    }
    template<typename Functor>
    INLINE void for_each_uncanon_unspaced_unwindowed(const Functor &func) {
        const KmerT mask(rhmask<KmerT>(rht, sp_.k_));
        schism::Schismatic<std::conditional_t<(sizeof(KmerT) <= 8), KmerT, uint64_t>> div(mask);
        KmerT min;
        unsigned filled;
        const size_t mul = rhmul();
        int8_t nv;
        loop_start:
        min = filled = 0;
        while(likely(pos_ < l_)) {
            while(filled < sp_.k_ && likely(pos_ < l_)) {
                nv = lutptr[s_[pos_++]];
                if(nv == int8_t(-1)) {min = ENCODE_OVERFLOW;goto loop_start;}
                min = (min * mul) | nv;
                ++filled;
            }
            if(likely(filled == sp_.k_)) {
                if(rht == DNA || rht == DNA2 || rht == PROTEIN_3BIT) min &= mask;
                else {
                    assert(div.mod(min) == min % mask);
                    CONST_IF(sizeof(KmerT) <= 8) {
                        min = div.mod(min);
                    } else {
                        min %= mask;
                    }
                }
                func(min);
                --filled;
            }
        }
    }
    template<typename Functor>
    INLINE void for_each_uncanon_unspaced_windowed(const Functor &func) {
        const KmerT mask(rhmask<KmerT>(rht, sp_.k_));
        schism::Schismatic<std::conditional_t<(sizeof(KmerT) <= 8), KmerT, uint64_t>> div(mask);
        KmerT min, kmer;
        unsigned filled;
        const size_t mul = rhmul();
        windowed_loop_start:
        min = filled = 0;
        while(likely(pos_ < l_)) {
            while(filled < sp_.k_ && likely(pos_ < l_)) {
                min *= mul;
                if(unlikely((min |= lutptr[s_[pos_++]]) == ENCODE_OVERFLOW) && likely(sp_.k_ < sizeof(KmerT) * 4 || lutptr[s_[pos_ - 1]] != 'T')) {
                    goto windowed_loop_start;
                }
                ++filled;
            }
            if(likely(filled == sp_.k_)) {
                if(rht == DNA || rht == DNA2 || rht == PROTEIN_3BIT) min &= mask;
                else {
                    assert(div.mod(min) == min % mask);
                    CONST_IF(sizeof(KmerT) <= 8) {
                        min = div.mod(min);
                    } else {
                        min %= mask;
                    }
                }
                if((kmer = qmap_.next_value(min, scorer_(min, getdata()))) != ENCODE_OVERFLOW) func(kmer);
                --filled;
            }
        }
        if(qmap_.partially_full())
            func(qmap_.max_in_queue().el_);
    }
    template<typename Functor>
    INLINE void for_each_uncanon_unspaced_windowed_entropy_(const Functor &func) {
        // NEVER CALL THIS DIRECTLY.
        // This contains instructions for generating uncanonicalized but windowed entropy-minimized kmers.
        const KmerT mask(rhmask<KmerT>(rht, sp_.k_));
        KmerT min, kmer;
        unsigned filled;
        schism::Schismatic<std::conditional_t<(sizeof(KmerT) <= 8), KmerT, uint64_t>> div(mask);
        if(!ent_tracker_)
            ent_tracker_.reset(new CircusEnt(this->k()));
        CircusEnt &ent = *ent_tracker_;
        windowed_loop_start:
        ent.clear();
        const size_t mul = rhmul();
        filled = min = 0;
        while(likely(pos_ < l_)) {
            while(filled < sp_.k_ && likely(pos_ < l_)) {
                const auto nc = lutptr[s_[pos_++]];
                if(nc == int8_t(-1)) {min = ENCODE_OVERFLOW; goto windowed_loop_start;}
                min = (mul * min) | nc;
                ent.push(nc);
                ++filled;
            }
            if(likely(filled == sp_.k_)) {
                if(rht == DNA || rht == DNA2 || rht == PROTEIN_3BIT) min &= mask;
                else {
                    CONST_IF(sizeof(KmerT) <= 8) {
                        assert(div.mod(min) == min % mask);
                        min = div.mod(min);
                    } else {
                        min %= mask;
                    }
                }
                if((kmer = qmap_.next_value(min, min / (ent.value() + .001))) != ENCODE_OVERFLOW) func(kmer);
                --filled;
            }
        }
        if(qmap_.partially_full())
            func(max_in_queue().el_);
    }
    template<typename Functor>
    INLINE void for_each_canon_unspaced_windowed_entropy_(const Functor &func) {
        this->for_each_uncanon_unspaced_windowed_entropy_([&,rht=this->rht](KmerT min) {
            if(rht == DNA) min = canonical_representation(min, sp_.k_);
            return func(min);
        });
    }
    // Utility 'for-each'-like functions.
    template<typename Functor>
    INLINE void for_each_hash(const Functor &func, const char *str, u64 l, unsigned k = 0) {
        s_ = str; l_ = l;
        for_each_hash<Functor>(func, k > 0 ? k: sp_.k_);
    }
    template<typename Functor>
    INLINE void for_each_hash(const Functor &func, unsigned k = 0) const {
        k = k > 0 ? k: sp_.k_;
        if(!sp_.unwindowed()) UNRECOVERABLE_ERROR("Can't for_each_hash for a windowed spacer");
        if(!sp_.unspaced()) UNRECOVERABLE_ERROR("Can't for_each_hash for a spaced spacer");
        if(l_ < k) return;
        size_t i = 0;
        uint64_t fhv=0, rhv=0, hv;
        const char *p, *p2;

        start:
        p = s_ + i;
        while(*p && cstr_lut[*p] < 0) ++p;
        for(;;) {
            p2 = p;
            if(*p2 == 0) return;
            while(*p2 && cstr_lut[*p2] >= 0 and p2 - p < k) ++p2;
            if(*p2 == 0) return;
            if(p2 - p == k) break;
            p = p2 + 1;
        }
        i = p - s_;
        hv = NTC64(s_ + i, k, fhv, rhv);
        func(canonicalize_ ? hv: fhv);
        for(; i < l_ - k; ++i) {
            auto newc = s_[i + k];
            if(cstr_lut[newc] < 0) {
                i += k;
                fhv = rhv = 0;
                goto start;
            }
            hv = NTC64(s_[i], newc, k, fhv, rhv);
            func(canonicalize_ ? hv: fhv);
        }
    }
    template<typename Functor>
    INLINE void for_each_hash(const Functor &func, kseq_t *ks) {
        while(kseq_read(ks) >= 0) assign(ks), for_each_hash<Functor>(func, ks->seq.s, ks->seq.l);
    }
    template<typename Functor>
    void for_each_hash(const Functor &func, gzFile fp, kseq_t *ks=nullptr) {
        bool destroy;
        if(ks == nullptr) ks = kseq_init(fp), destroy = true;
        else            kseq_assign(ks, fp), destroy = false;
        for_each_hash<Functor>(func, ks);
        if(destroy) kseq_destroy(ks);
    }
    template<typename Functor>
    void for_each_hash(const Functor &func, const char *path, kseq_t *ks=nullptr) {
        gzFile fp(gzopen(path, "rb"));
        if(!fp) UNRECOVERABLE_ERROR(ks::sprintf("Could not open file at %s. Abort!\n", path).data());
        gzbuffer(fp, 1<<18);
        for_each_hash<Functor>(func, fp, ks);
        gzclose(fp);
    }
    template<typename Functor>
    INLINE void for_each(const Functor &func, const char *str, u64 l) {
        this->assign(str, l);
        if(!has_next_kmer()) return;
        if(rht != DNA && canonicalize_) {canonicalize_ = false;}
        if(canonicalize_) {
            if(sp_.unwindowed()) {
                 for_each_canon_unwindowed(func);
            } else {
                if(std::is_same<ScoreType, score::Entropy>::value) {
                    if(sp_.unspaced()) for_each_canon_unspaced_windowed_entropy_(func);
                    else               for_each_canon_windowed(func);
                } else for_each_canon_windowed(func);
            }
        } else {
            if(sp_.unspaced()) {
                if(sp_.unwindowed()) for_each_uncanon_unspaced_unwindowed(func);
                else {
                    if(std::is_same<ScoreType, score::Entropy>::value)
                        for_each_uncanon_unspaced_windowed_entropy_(func);
                    else for_each_uncanon_unspaced_windowed(func);
                }
            } else {
                if(canonicalize_) for_each_uncanon_spaced(func);
                //else              for_each_canon_spaced(func); Unless the spaced seed is symmetric, we can't do this
            }
        }
    }
    template<typename Functor>
    INLINE void for_each(const Functor &func, kseq_t *ks) {
        while(kseq_read(ks) >= 0) assign(ks), for_each<Functor>(func, ks->seq.s, ks->seq.l);
    }
    template<typename Functor>
    INLINE void for_each_canon(const Functor &func, kseq_t *ks) {
        if(sp_.unwindowed()) while(kseq_read(ks) >= 0) assign(ks), for_each_canon_unwindowed<Functor>(func);
        else                 while(kseq_read(ks) >= 0) assign(ks), for_each_canon_windowed<Functor>(func);
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
        if(!fp) UNRECOVERABLE_ERROR(ks::sprintf("Could not open file at %s. Abort!\n", path).data());
        gzbuffer(fp, 1<<18);
        for_each_canon<Functor>(func, fp, ks);
        gzclose(fp);
    }
    template<typename Functor>
    void for_each_uncanon(const Functor &func, const char *path, kseq_t *ks=nullptr) {
        gzFile fp(gzopen(path, "rb"));
        if(!fp) UNRECOVERABLE_ERROR(ks::sprintf("Could not open file at %s. Abort!\n", path).data());
        gzbuffer(fp, 1<<18);
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
        const size_t pl = std::strlen(path);
        std::FILE *pfp = 0;
        gzFile fp = 0;
        bool matchxz = pl >= 3 && std::equal(&path[pl - 3], &path[pl], ".xz");
        bool matchbz = pl >= 4 && std::equal(&path[pl - 4], &path[pl], ".bz2");
        bool matchzst = pl >= 4 && std::equal(&path[pl - 4], &path[pl], ".zst");
        if(matchxz || matchbz || matchzst) {
            std::string cmd = std::string(matchxz ? "xz": (matchbz ? "bzip2": "zstd")) + " -dc " + path;
            pfp = ::popen(cmd.data(), "r");
            if(!pfp) UNRECOVERABLE_ERROR(std::string("Failed to open popen call: ") + cmd);
            fp = gzdopen(::fileno(pfp), "rb");
        } else fp = gzopen(path, "rb");
        if(!fp) UNRECOVERABLE_ERROR(ks::sprintf("Could not open file at %s. Abort!\n", path).data());
        gzbuffer(fp, 1<<18);
        if(canonicalize_) for_each_canon<Functor>(func, fp, ks);
        else              for_each_uncanon<Functor>(func, fp, ks);
        gzclose(fp);
        if(pfp) ::pclose(pfp);
    }
    template<typename Functor, typename ContainerType,
             typename=typename std::enable_if<std::is_same<typename ContainerType::value_type::value_type, char>::value ||
                                       std::is_same<typename std::decay<typename ContainerType::value_type>::type, char *>::value
                                             >::type
            >
    void for_each(const Functor &func, const ContainerType &strcon, kseq_t *ks=nullptr) {
        for(const auto &el: strcon) {
            for_each<Functor>(func, get_cstr(el), ks);
        }
    }

    size_t rhmul() const {
        return mul(rht);
    }

    // Encodes a kmer starting at `start` within string `s_`.
    INLINE KmerT kmer(unsigned start) {
        assert(start <= l_ - sp_.c_ + 1);
        if(l_ < sp_.c_)    return ENCODE_OVERFLOW;
        KmerT new_kmer(lutptr[s_[start]]);
        if(new_kmer == ENCODE_OVERFLOW) return ENCODE_OVERFLOW;
        u64 len(sp_.s_.size());
        int8_t nextc;
        auto spaces(sp_.s_.data());
        static constexpr bool isent = std::is_same_v<ScoreType, score::Entropy>;
        CircusEnt *ent_tracker = 0;
        CONST_IF(isent) {
            if(!ent_tracker_)
                ent_tracker_.reset(new CircusEnt(this->k()));
            ent_tracker = ent_tracker_.get();
            ent_tracker->clear();
        }
        if(rht == DNA || rht == PROTEIN || rht == PROTEIN_3BIT || rht == DNA2 || rht == DNAC) {
            const int shift = rht == DNA ? 2: rht == PROTEIN ? 8: (rht == DNA2 || rht == DNAC) ? 1: 3;
            /*std::fprintf(stderr, "Shift %d\n", shift);*/
#define ITER do {\
            start += *spaces++;\
            if((nextc = lutptr[s_[start]]) == int8_t(-1)) {\
                new_kmer = ENCODE_OVERFLOW;\
                goto rnk;\
            }\
            new_kmer = (new_kmer << shift) | nextc;\
            CONST_IF(isent) ent_tracker->push(nextc);\
        } while(0);
            DO_DUFF(len, ITER);
#undef ITER
        } else if(rht == PROTEIN20 || rht == PROTEIN_14 || rht == PROTEIN_6) {
            const size_t mul = rht == PROTEIN20 ? 20: rht == PROTEIN_14 ? 14: 6;
#define ITER do {new_kmer *= mul, start += *spaces++;\
            if((nextc = lutptr[s_[start]]) == int8_t(-1)) {\
                new_kmer = ENCODE_OVERFLOW;\
                goto rnk;\
            }\
            new_kmer = (new_kmer * mul) | nextc;\
            CONST_IF(isent) ent_tracker->push(nextc);\
        } while(0);
            DO_DUFF(len, ITER);
#undef ITER
        }
        rnk:
        return new_kmer;
    }
    // Whether or not an additional kmer is present in the sequence being encoded.
    INLINE int has_next_kmer() const {
        static_assert(std::is_same<decltype((std::int64_t)l_ - sp_.c_ + 1), std::int64_t>::value, "is not same");
        return (pos_ + sp_.c_ - 1) < l_;
    }
    // This fetches our next kmer for our window. It is immediately placed in the qmap_t,
    // which is a tree map containing kmers and scores so we can keep track of the best-scoring
    // kmer in the window.
    INLINE KmerT next_kmer() {
        assert(has_next_kmer());
        return kmer(pos_++);
    }
    // This is the actual point of entry for fetching our minimizers.
    // It wraps encoding and scoring a kmer, updates qmap, and returns the minimizer
    // for the next window.
    static constexpr bool is_entropy = std::is_same_v<ScoreType, score::Entropy>;
    void *getdata() const {
        if constexpr(is_entropy) {
            return ent_tracker_.get();
        } else {
            return data_;
        }
    }
    INLINE KmerT next_minimizer() {
        //if(unlikely(!has_next_kmer())) return ENCODE_OVERFLOW;
        const KmerT k(kmer(pos_++));
        const KmerT kscore(scorer_(k, getdata()));
        return qmap_.next_value(k, kscore);
    }
    INLINE KmerT next_canonicalized_minimizer() {
        assert(has_next_kmer());
        KmerT nk = kmer(pos_++);
        if(rht == DNA) nk = canonical_representation(nk, sp_.k_);
        const KmerT kscore(scorer_(nk, getdata()));
        return qmap_.next_value(nk, kscore);
    }
    auto max_in_queue() const {return qmap_.begin()->first;}

    bool canonicalize() const {return canonicalize_;}
    void canonicalize(bool value) {canonicalize_ = value;}

    auto pos() const {return pos_;}
    void pos(uint64_t v) {pos_ = v;}
    uint32_t k() const {return sp_.k_;}
    size_t n_in_queue() const {return qmap_.n_in_queue();}
};





template<typename IntType, typename HashClass=CyclicHash<IntType>>
struct RollingHasher {
    static_assert(std::is_integral<IntType>::value || sizeof(IntType) > 8, "Must be integral (or by uint128/int128)");
    long long int k_;
    long long int w_;
private:
    InputType enctype_;
public:
    bool canon_;
    HashClass hasher_;
    HashClass rchasher_;
    uint64_t seed1_, seed2_;
    QueueMap<IntType, uint64_t> qmap_;
    const int8_t *lutptr = (const int8_t *)cstr_lut;
    long long int window() const {return w_;}
    InputType hashtype() const {return enctype_;}
    RollingHasher &hashtype(InputType rht) {
        enctype_ = rht; lutptr = rh2lp(rht);
        return *this;
    }
    void window(long long int w) {
        if(w <= k_) w_ = -1;
        else w_ = w;
        if(w_ > 0) {
            qmap_.resize(w_ - k_ + 1);
        }
    }
    static constexpr IntType ENCODE_OVERFLOW = static_cast<IntType>(-1);
    RollingHasher(unsigned k=21, bool canon=false,
                   InputType enc=DNA, long long int wsz = -1, uint64_t seed1=1337, uint64_t seed2=137):
        k_(k), enctype_(enc), canon_(canon), hasher_(k, sizeof(IntType) * CHAR_BIT), rchasher_(k, sizeof(IntType) * CHAR_BIT)
        , seed1_(seed1), seed2_(seed2)
    {
        if(canon_ && enc != InputType::DNA) {
            std::fprintf(stderr, "Note: RollingHasher with Protein alphabet does not support reverse-complementing.\n");
            canon_ = false;
        }
        window(wsz);
        hasher_.seed(seed1, seed2);
        rchasher_.seed(seed2 * seed1, seed2 ^ seed1);
        if(enc == PROTEIN_6_FRAME) throw NotImplementedError("Protein 6-frame not implemented.");
    }
    RollingHasher& operator=(const RollingHasher &o) {
        k_ = o.k_; canon_ = o.canon_; enctype_ = o.enctype_; w_ = o.w_; seed1_ = o.seed1_; seed2_ = o.seed2_;
        window(o.w_);
        return *this;
    }
    RollingHasher(const RollingHasher &o): RollingHasher(o.k_, o.canon_, o.enctype_, o.w_, o.seed1_, o.seed2_) {}
    template<typename Functor>
    void for_each_canon(const Functor &func, const char *s, size_t l) {
        qmap_.reset();
        if(enctype_ != DNA) {
            for_each_uncanon<Functor>(func, s, l);
            return;
        }
        if(l < uint32_t(k_)) return;
        hasher_.reset();
        rchasher_.reset();
        size_t i;
        long long int nf;
        uint8_t v1;
        IntType nextv;
        if(qmap_.size() > 1) {
            auto add_hashes =  [&](auto &hasher) {
            if((nextv = qmap_.next_value(hasher.hashvalue, lex_score(hasher.hashvalue))) != ENCODE_OVERFLOW)
                func(nextv);
            };
            for(i = nf = 0; nf < k_ && i < l; ++i) {
                if((v1 = cstr_lut[s[i]]) == uint8_t(-1)) {
                    fixup_minimizer:
                    if(i + 2 * k_ >= l) goto end;
                    i += k_;
                    nf = 0;
                    hasher_.reset();
                    rchasher_.reset();
                } // Fixme: this ignores both strands when one becomes 'N'-contaminated.
                  // In the future, encode the side that is still valid
                else hasher_.eat(v1), rchasher_.eat(cstr_rc_lut[s[i - nf + k_ - 1]]), ++nf;
            }
            if(nf < k_) goto end; // All failed
            add_hashes(hasher_);
            add_hashes(rchasher_);
            for(;i < l; ++i) {
                if((v1 = cstr_lut[s[i]]) == uint8_t(-1))
                    goto fixup_minimizer;
                hasher_.update(cstr_lut[s[i - k_]], v1);
                rchasher_.reverse_update(cstr_rc_lut[s[i]], cstr_rc_lut[s[i - k_]]);
                add_hashes(hasher_);
                add_hashes(rchasher_);
            }
            end:
            if(qmap_.partially_full())
                func(max_in_queue().el_);
        } else {
            for(i = nf = 0; nf < k_ && i < l; ++i) {
                if((v1 = cstr_lut[s[i]]) == uint8_t(-1)) {
                    fixup:
                    if(i + 2 * k_ >= l) return;
                    i += k_;
                    nf = 0;
                    hasher_.reset();
                    rchasher_.reset();
                } // Fixme: this ignores both strands when one becomes 'N'-contaminated.
                  // In the future, encode the side that is still valid
                else hasher_.eat(v1), rchasher_.eat(cstr_rc_lut[s[i - nf + k_ - 1]]), ++nf;
            }
            if(nf < k_) return; // All failed
            func(std::min(hasher_.hashvalue, rchasher_.hashvalue));
            for(;i < l; ++i) {
                if((v1 = cstr_lut[s[i]]) == uint8_t(-1))
                    goto fixup;
                hasher_.update(cstr_lut[s[i - k_]], v1);
                rchasher_.reverse_update(cstr_rc_lut[s[i]], cstr_rc_lut[s[i - k_]]);
                func(std::min(hasher_.hashvalue, rchasher_.hashvalue));
            }
        }
    }

    template<typename Functor>
    void for_each_uncanon(const Functor &func, const char *s, size_t l) {
        if(l < size_t(k_)) return;
        hasher_.reset();
        qmap_.reset();
        size_t i;
        long long int nf;
        int8_t v1;
        IntType nextv;
        auto use_val = [&](auto v) {
            if(qmap_.size() > 1) {
                if((nextv = qmap_.next_value(v, lex_score(v))) != ENCODE_OVERFLOW) {
                    func(nextv);
                }
            } else if(v != ENCODE_OVERFLOW) func(v);
        };
        for(i = nf = 0; nf < k_ && i < l; ++i) {
            if(unlikely((v1 = lutptr[s[i]]) == int8_t(-1))) {
                //std::fprintf(stderr, "Char %c/%d was missing... %d\n", s[i], s[i], lutptr[s[i]]);
                fixup:
                i += k_; nf = 0; hasher_.reset();
            } else hasher_.eat(v1), ++nf;
        }
        if(nf < k_) return; // All failed
        use_val(hasher_.hashvalue);
        for(;i < l; ++i) {
            if(lutptr[s[i]] == int8_t(-1)) goto fixup;
            //auto ov = hasher_.hashvalue;
            hasher_.update(lutptr[s[i - k_]], lutptr[s[i]]);
            //std::fprintf(stderr, "Updating with new char %c, which is translated to %d, which will replcae old %zu with %zu\n", s[i], int(lutptr[s[i]]), size_t(ov), size_t(hasher_.hashvalue));
            use_val(hasher_.hashvalue);
        }
        if(qmap_.partially_full())
            func(max_in_queue().el_);
    }
    template<typename Functor>
    void for_each_canon(const Functor &func, kseq_t *ks) {
        while(kseq_read(ks) >= 0) {
            for_each_canon<Functor>(func, ks->seq.s, ks->seq.l);
        }
    }
    template<typename Functor>
    void for_each_uncanon(const Functor &func, kseq_t *ks) {
        while(kseq_read(ks) >= 0) {
            for_each_uncanon<Functor>(func, ks->seq.s, ks->seq.l);
        }
    }
    template<typename Functor>
    void for_each_hash(const Functor &func, const char *s, size_t l) {
        if(canon_) for_each_canon<Functor>(func, s, l);
        else       for_each_uncanon<Functor>(func, s, l);
    }
    template<typename Functor>
    void for_each_hash(const Functor &func, gzFile fp, kseq_t *ks=nullptr) {
        if(canon_) for_each_canon<Functor>(func, fp, ks);
        else       for_each_uncanon<Functor>(func, fp, ks);
    }
    template<typename Functor>
    void for_each_hash(const Functor &func, const char *inpath, kseq_t *ks=nullptr) {
        const size_t pl = std::strlen(inpath);
        std::FILE *pfp = 0;
        gzFile fp = 0;
        bool matchxz = pl >= 3 && std::equal(&inpath[pl - 3], &inpath[pl], ".xz");
        bool matchbz = pl >= 4 && std::equal(&inpath[pl - 4], &inpath[pl], ".bz2");
        bool matchzst = pl >= 4 && std::equal(&inpath[pl - 4], &inpath[pl], ".zst");
        if(matchxz || matchbz || matchzst) {
            std::string cmd = std::string(matchxz ? "xz": (matchbz ? "bzip2": "zstd")) + " -dc " + inpath;
            pfp = ::popen(cmd.data(), "r");
            if(!pfp) UNRECOVERABLE_ERROR(std::string("Failed to open popen call: ") + cmd);
            fp = gzdopen(::fileno(pfp), "rb");
        } else fp = gzopen(inpath, "rb");
        if(!fp) UNRECOVERABLE_ERROR(std::string("Could not open file at ") + inpath);
        gzbuffer(fp, 1<<18);
        for_each_hash<Functor>(func, fp, ks);
        gzclose(fp);
        if(pfp) ::pclose(pfp);
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
    void for_each_canon(const Functor &func, gzFile fp, kseq_t *ks=nullptr) {
        bool destroy;
        if(ks == nullptr) ks = kseq_init(fp), destroy = true;
        else            kseq_assign(ks, fp), destroy = false;
        for_each_canon<Functor>(func, ks);
        if(destroy) kseq_destroy(ks);
    }
    template<typename...Args>
    void for_each(Args &&...args) {
        for_each_hash(std::forward<Args>(args)...);
    }
    void reset() {hasher_.reset(); rchasher_.reset();}
    size_t n_in_queue() const {return qmap_.n_in_queue();}
    auto max_in_queue() const {return qmap_.begin()->first;}
    bool canonicalize() const {return canon_;}
    void canonicalize(bool value) {canon_ = value;}
};

template<typename IType>
struct RollingHasherSet {
    std::vector<RollingHasher<IType>> hashers_;
    bool canon_;
    template<typename C>
    RollingHasherSet(const C &c, bool canon=false, InputType enc=DNA, uint64_t seedseed=1337u): canon_(canon) {
        std::mt19937_64 mt(seedseed);
        hashers_.reserve(c.size());
        for(const auto k: c)
            hashers_.emplace_back(k, canon, enc, mt(), mt());
    }
    template<typename Functor>
    void for_each_canon(const Functor &func, const char *s, size_t l) {
        const auto mink = get_mink();
        if(l < mink) return;
        for(auto &h: hashers_) h.reset();
        size_t i = 0;
        long long int nf = 0;
        uint8_t v1;
        for(; nf < mink && i < l; ++i) {
            if((v1 = cstr_lut[s[i]]) == uint8_t(-1)) {
                fixup:
                if(i + 2 * mink >= l) return;
                i += mink;
                nf = 0;
                for(auto &h: hashers_) h.reset();
            } // Fixme: this ignores both strands when one becomes 'N'-contaminated.
              // In the future, encode the side that is still valid
            else {
                for(auto &h: hashers_) h.hasher_.eat(v1), h.rchasher_.eat(cstr_rc_lut[s[h.k_ - i - 1]]);
                ++nf;
            }
        }
        for(size_t hi = 0; hi < hashers_.size(); ++hi) {
            auto &h(hashers_[hi]);
            if(nf >= h.k_)
                func(std::min(h.hasher_.hashvalue, h.rchasher_.hashvalue), hi);
        }
        for(;i < l; ++i) {
            if((v1 = cstr_lut[s[i]]) == uint8_t(-1))
                goto fixup;
            for(size_t hi = 0; hi < hashers_.size(); ++hi) {
                auto &h(hashers_[hi]);
                //h.rchasher_.eat(cstr_rc_lut[s[i - nf + h.k_ - 1]]);
                if(nf >= h.k_) {
                    h.rchasher_.reverse_update(cstr_rc_lut[s[i]], cstr_rc_lut[s[i - h.k_]]);
                    h.hasher_.update(cstr_lut[s[i - h.k_]], v1);
                    func(std::min(h.hasher_.hashvalue, h.rchasher_.hashvalue), hi);
                } else {
                    h.hasher_.eat(v1);
                    h.rchasher_.eat(cstr_rc_lut[s[i + h.k_ - nf - 1]]);
                }
            }
            ++nf;
        }
    }
    uint32_t get_mink() const {return std::accumulate(hashers_.begin(), hashers_.end(), unsigned(-1), [](unsigned x, const auto & y) {return std::min(x, unsigned(y.k_));});}
    template<typename Functor>
    void for_each_uncanon(const Functor &func, const char *s, size_t l) {
        const auto mink = get_mink();
        if(l < mink) return;
        for(auto &h: hashers_) h.reset();
        size_t i = 0;
        long long int nf = 0;
        uint8_t v1;
        for(; nf < mink && i < l; ++i) {
            if((v1 = cstr_lut[static_cast<uint8_t>(s[i])]) == uint8_t(-1)) {
                fixup:
                if(i + 2 * mink >= l) return;
                i += mink;
                nf = 0;
                for(auto &h: hashers_) h.hasher_.reset();
            } // Fixme: this ignores both strands when one becomes 'N'-contaminated.
              // In the future, encode the side that is still valid
            else {
                for(auto &h: hashers_) h.hasher_.eat(v1);
                ++nf;
            }
        }
        for(size_t i = 0; i < hashers_.size(); ++i) {
            auto &h(hashers_[i]);
            if(nf >= h.k_)
                func(h.hasher_.hashvalue, i);
        }
        for(;i < l; ++i) {
            if((v1 = cstr_lut[uint8_t(s[i])]) == uint8_t(-1))
                goto fixup;
            for(size_t hi = 0; hi < hashers_.size(); ++hi) {
                auto &h(hashers_[hi]);
                h.hasher_.eat(v1);
                if(++nf >= h.k_)
                    func(h.hasher_.hashvalue, hi);
            }
        }
    }
    template<typename Functor>
    void for_each_canon(const Functor &func, kseq_t *ks) {
        while(kseq_read(ks) >= 0) {
            for_each_canon<Functor>(func, ks->seq.s, ks->seq.l);
        }
    }
    template<typename Functor>
    void for_each_uncanon(const Functor &func, kseq_t *ks) {
        while(kseq_read(ks) >= 0) {
            for_each_uncanon<Functor>(func, ks->seq.s, ks->seq.l);
        }
    }
    template<typename Functor>
    void for_each_hash(const Functor &func, gzFile fp, kseq_t *ks=nullptr) {
        if(canon_) for_each_canon<Functor>(func, fp, ks);
        else       for_each_uncanon<Functor>(func, fp, ks);
    }
    template<typename Functor>
    void for_each_hash(const Functor &func, const char *inpath, kseq_t *ks=nullptr) {
        gzFile fp = gzopen(inpath, "rb");
        if(!fp) throw file_open_error(inpath);
        gzbuffer(fp, 1<<18);
        for_each_hash<Functor>(func, fp, ks);
        gzclose(fp);
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
    void for_each_canon(const Functor &func, gzFile fp, kseq_t *ks=nullptr) {
        bool destroy;
        if(ks == nullptr) ks = kseq_init(fp), destroy = true;
        else            kseq_assign(ks, fp), destroy = false;
        for_each_canon<Functor>(func, ks);
        if(destroy) kseq_destroy(ks);
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
    enc.for_each([ret](auto x) {int khr; auto it = kh_get(all, ret, x); if(it == ret->n_buckets) {kh_put(all, ret, x, &khr); if(khr < 0) throw std::runtime_error("Error adding to hash table");}}
                 , path.data());
    return ret;
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
                          LOG_WARNING("System error: resource temporarily available. Retry #%i\n", tries + 1);
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


template<typename ScoreType, typename SketchType>
void fill_lmers(SketchType &sketch, const std::string &path, const Spacer &space, bool canonicalize=true,
                void *data=nullptr, kseq_t *ks=nullptr) {
#if USE_HASH_FILLER
    detail::HashFiller<SketchType> hf(sketch);
#endif
    sketch.not_ready();
    Encoder<ScoreType> enc(nullptr, 0, space, data, canonicalize);
#if USE_HASH_FILLER
    enc.for_each([&](u64 min) {hf.add(min);}, path.data(), ks);
#else
    enc.for_each([&](u64 min) {sketch.addh(min);}, path.data(), ks);
#endif
}


#if 0
#endif

} //namespace bns
#endif // _EMP_ENCODER_H__
