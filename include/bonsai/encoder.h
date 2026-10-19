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
#include "polyhash.h"
#include "ntHash/nthash.hpp"
#include "alphabet.h"
#include "rhtraits.h"
#include "sixframe.h"
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
static INLINE u64 lex_score(u64 i) {return FRev64()(i);}
static INLINE u64 lex_score(u64 i, void *) {return FRev64()(i);}
// Entropy-weighted score (lower wins): the k-mer's hash divided by the entropy
// of its characters, so that high-entropy k-mers are preferred. data points to
// a CircusEnt holding the k-mer's characters. Only the top 48 bits of the hash
// are used so that the quotient fits in 64 bits even for zero entropy.
template<typename T>
static INLINE T ent_score(T i, void *data) {
    const double ent = std::max(reinterpret_cast<CircusEnt *>(data)->value(), 0.);
    const u64 h = u64(lex_score(i) >> (sizeof(T) * CHAR_BIT - 48));
    return T(u64(h / (ent + 1e-4)));
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
    QueueMap<KmerT, KmerT> qmap_; // monotone deque holding the window minimum (score, k-mer) pair; see qmap.h.
    const ScoreType  scorer_; // scoring struct
    bool canonicalize_;
    InputType rht = InputType::DNA;
    const int8_t *lutptr = (const int8_t *)DNA4.data();
    size_t nremper = rh2n(rht, sizeof(KmerT));
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
        if(sp_.k_ > nremper) {
            RUNTIME_ERROR(ks::sprintf("k %u exceeds the %zu bases representable in this k-mer type for alphabet %s.",
                                      sp_.k_, nremper, to_string(rht).data()).data());
        }
        if(std::is_same<ScoreType, score::Entropy>::value) {
            ent_tracker_.reset(new CircusEnt(sp_.k_));
        }
        if(!sp_.unspaced() && canonicalize_) {
            canonicalize_ = false;
        }
    }
    Encoder(const Spacer &sp, void *data, bool canonicalize=true): Encoder(nullptr, 0, sp, data, canonicalize) {}
    Encoder(const Spacer &sp, bool canonicalize=true): Encoder(sp, nullptr, canonicalize) {}
    Encoder(const Encoder &o): s_(o.s_), l_(o.l_), sp_(o.sp_), pos_(o.pos_), data_(o.data_), scorer_(o.scorer_), canonicalize_(o.canonicalize_), rht(o.rht), lutptr(o.lutptr), nremper(o.nremper) {
        if(sp_.w_ > sp_.c_)
            qmap_.resize(sp_.w_ - sp_.c_ + 1);
    }
    Encoder(Encoder<ScoreType, KmerT> &&o): s_(o.s_), l_(o.l_), sp_(o.sp_), pos_(o.pos_), data_(o.data_),
            qmap_(std::move(o.qmap_)), scorer_{}, canonicalize_(o.canonicalize_), rht(o.rht), lutptr(o.lutptr), nremper(o.nremper) {
        if(o.ent_tracker_) ent_tracker_.reset(new CircusEnt(*o.ent_tracker_));
    }
    Encoder &operator=(const Encoder<ScoreType, KmerT> &o) {
        s_ = o.s_; l_ = o.l_;
        sp_ = o.sp_;
        pos_ = o.pos_;
        data_ = o.data_;
        qmap_ = o.qmap_;
        canonicalize_ = o.canonicalize_;
        rht = o.rht; lutptr = o.lutptr; nremper = o.nremper;
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
        // A record with fewer k-mers than the window still yields its best one.
        if(qmap_.partially_full())
            func(max_in_queue().el_);
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
        // Windows run over the valid k-mers: a k-mer with an invalid character
        // is skipped rather than entering the window as ENCODE_OVERFLOW. All
        // ones is also poly-T when the k-mer fills KmerT, so validity is
        // rechecked in that case, and a full window emits whatever it holds.
        while(likely(has_next_kmer())) {
            const unsigned start = pos_++;
            const KmerT km = kmer(start);
            if(km == ENCODE_OVERFLOW && !kmer_is_valid(start)) continue;
            const KmerT min = qmap_.next_value(km, scorer_(km, getdata()));
            if(min != ENCODE_OVERFLOW || qmap_.n_in_queue() == qmap_.size())
                func(min);
        }
        // A record with fewer valid k-mers than the window still yields its best one.
        if(qmap_.partially_full())
            func(max_in_queue().el_);
    }
    // Appends one character code to a k-mer for the unspaced rolling encoders.
    // Bit-packed alphabets shift the new code in and mask off the oldest
    // character in finish(). Other alphabets (20, 14 or 6 letters) are encoded
    // in base mul, where OR is not addition, so the oldest character is removed
    // with a modulus by mul^(k-1) before multiplying; the value then stays below
    // mul^k and cannot overflow KmerT.
    struct RollingAppender {
        using DivT = std::conditional_t<(sizeof(KmerT) <= 8), KmerT, uint64_t>;
        bool bitpacked;
        KmerT mul, mask, headmod;
        schism::Schismatic<DivT> div;
        RollingAppender(const Encoder &enc):
            bitpacked(rh_bitpacked(enc.rht)), mul(enc.rhmul()),
            mask(rhmask<KmerT>(enc.rht, enc.sp_.k_)),
            headmod(bitpacked ? KmerT(1): rhmask<KmerT>(enc.rht, enc.sp_.k_ - 1)),
            div(sizeof(KmerT) <= 8 ? DivT(headmod): DivT(1))
        {}
        INLINE KmerT append(KmerT x, int8_t nv) const {
            if(bitpacked) return (x * mul) | KmerT(nv);
            CONST_IF(sizeof(KmerT) <= 8) {
                x = div.mod(x);
            } else {
                x %= headmod;
            }
            return x * mul + KmerT(nv);
        }
        INLINE KmerT finish(KmerT x) const {return bitpacked ? KmerT(x & mask): x;}
    };
    template<typename Functor>
    INLINE void for_each_uncanon_unspaced_unwindowed(const Functor &func) {
        const RollingAppender app(*this);
        KmerT min;
        unsigned filled;
        loop_start:
        min = filled = 0;
        while(likely(pos_ < l_)) {
            while(filled < sp_.k_ && likely(pos_ < l_)) {
                const uint8_t c_at_pos = s_[pos_];
                const int8_t nv = lutptr[c_at_pos];
                ++pos_;
                if(nv == int8_t(-1)) {min = ENCODE_OVERFLOW; goto loop_start;}
                min = app.append(min, nv);
                ++filled;
            }
            if(likely(filled == sp_.k_)) {
                min = app.finish(min);
                func(min);
                --filled;
            }
        }
    }
    template<typename Functor>
    INLINE void for_each_uncanon_unspaced_windowed(const Functor &func) {
        const RollingAppender app(*this);
        KmerT min, kmer;
        unsigned filled;
        windowed_loop_start:
        min = filled = 0;
        while(likely(pos_ < l_)) {
            while(filled < sp_.k_ && likely(pos_ < l_)) {
                const int8_t nv = lutptr[uint8_t(s_[pos_++])];
                if(unlikely(nv == int8_t(-1))) goto windowed_loop_start;
                min = app.append(min, nv);
                ++filled;
            }
            if(likely(filled == sp_.k_)) {
                min = app.finish(min);
                // QueueMap returns all ones for a window that is not full, but all
                // ones is also poly-T when the k-mer fills KmerT, so a full window
                // emits its minimizer whatever its value.
                if((kmer = qmap_.next_value(min, scorer_(min, getdata()))) != ENCODE_OVERFLOW || qmap_.n_in_queue() == qmap_.size())
                    func(kmer);
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
        const RollingAppender app(*this);
        KmerT min, kmer;
        unsigned filled;
        if(!ent_tracker_)
            ent_tracker_.reset(new CircusEnt(this->k()));
        CircusEnt &ent = *ent_tracker_;
        windowed_loop_start:
        ent.clear();
        filled = min = 0;
        while(likely(pos_ < l_)) {
            while(filled < sp_.k_ && likely(pos_ < l_)) {
                const auto nc = lutptr[uint8_t(s_[pos_++])];
                if(nc == int8_t(-1)) {min = ENCODE_OVERFLOW; goto windowed_loop_start;}
                min = app.append(min, nc);
                ent.push(nc);
                ++filled;
            }
            if(likely(filled == sp_.k_)) {
                min = app.finish(min);
                if((kmer = qmap_.next_value(min, ent_score(min, &ent))) != ENCODE_OVERFLOW || qmap_.n_in_queue() == qmap_.size())
                    func(kmer);
                --filled;
            }
        }
        if(qmap_.partially_full())
            func(max_in_queue().el_);
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
        while(*p && cstr_lut[uint8_t(*p)] < 0) ++p;
        for(;;) {
            p2 = p;
            if(*p2 == 0) return;
            while(*p2 && cstr_lut[uint8_t(*p2)] >= 0 and p2 - p < k) ++p2;
            if(*p2 == 0) return;
            if(p2 - p == k) break;
            p = p2 + 1;
        }
        i = p - s_;
        hv = NTC64(s_ + i, k, fhv, rhv);
        func(canonicalize_ ? hv: fhv);
        for(; i < l_ - k; ++i) {
            auto newc = s_[i + k];
            if(cstr_lut[uint8_t(newc)] < 0) {
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
            // Windows are scored on canonical k-mers, as in the file API, so
            // that a sequence and its reverse complement select the same k-mers.
            if(sp_.unwindowed()) for_each_canon_unwindowed(func);
            else                 for_each_canon_windowed(func);
        } else {
            if(sp_.unspaced()) {
                if(sp_.unwindowed()) for_each_uncanon_unspaced_unwindowed(func);
                else {
                    if(std::is_same<ScoreType, score::Entropy>::value)
                        for_each_uncanon_unspaced_windowed_entropy_(func);
                    else for_each_uncanon_unspaced_windowed(func);
                }
            } else {
                // Spaced seeds are never canonicalized: unless the seed is
                // symmetric, a k-mer and its reverse complement use different positions.
                for_each_uncanon_spaced(func);
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
        const bool us = sp_.unspaced(), spu = sp_.unwindowed();
        while(kseq_read(ks) >= 0) {
            assign(ks);
            if(us) {
                if(spu)             for_each_uncanon_unspaced_unwindowed(func);
                else if(is_entropy) for_each_uncanon_unspaced_windowed_entropy_(func);
                else                for_each_uncanon_unspaced_windowed(func);
            } else {
                for_each_uncanon_spaced(func);
            }
        }
    }
    template<typename Functor>
    INLINE void for_each_canon(const Functor &func, gzFile fp, kseq_t *ks=nullptr) {
        bool destroy;
        if(ks == nullptr) ks = kseq_init(fp), destroy = true;
        else            kseq_assign(ks, fp), destroy = false;
        for_each_canon<Functor>(func, ks);
        if(destroy) kseq_destroy(ks);
    }
    template<typename Functor>
    INLINE void for_each_uncanon(const Functor &func, gzFile fp, kseq_t *ks=nullptr) {
        bool destroy;
        if(ks == nullptr) ks = kseq_init(fp), destroy = true;
        else            kseq_assign(ks, fp), destroy = false;
        for_each_uncanon<Functor>(func, ks);
        if(destroy) kseq_destroy(ks);
    }
    template<typename Functor>
    INLINE void for_each_canon(const Functor &func, const char *path, kseq_t *ks=nullptr) {
        gzFile fp(gzopen(path, "rb"));
        if(!fp) UNRECOVERABLE_ERROR(ks::sprintf("Could not open file at %s. Abort!\n", path).data());
        gzbuffer(fp, 1<<18);
        for_each_canon<Functor>(func, fp, ks);
        gzclose(fp);
    }
    template<typename Functor>
    INLINE void for_each_uncanon(const Functor &func, const char *path, kseq_t *ks=nullptr) {
        gzFile fp(gzopen(path, "rb"));
        if(!fp) UNRECOVERABLE_ERROR(ks::sprintf("Could not open file at %s. Abort!\n", path).data());
        gzbuffer(fp, 1<<18);
        for_each_uncanon<Functor>(func, fp, ks);
        gzclose(fp);
    }
    template<typename Functor>
    INLINE void for_each(const Functor &func, gzFile fp, kseq_t *ks=nullptr) {
        bool destroy;
        if(ks == nullptr) ks = kseq_init(fp), destroy = true;
        else            kseq_assign(ks, fp), destroy = false;
        if(canonicalize_) for_each_canon<Functor>(func, ks);
        else              for_each_uncanon<Functor>(func, ks);
        if(destroy) kseq_destroy(ks);
    }
    template<typename Functor>
    INLINE void for_each(const Functor &func, const std::string &path, kseq_t *ks=nullptr) {
        for_each(func, path.data(), ks);
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
        KmerT new_kmer(lutptr[uint8_t(s_[start])]);
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
            ent_tracker->push(int8_t(new_kmer));
        }
        if(rht == DNA || rht == PROTEIN || rht == PROTEIN_3BIT || rht == DNA2 || rht == DNAC) {
            const int shift = rht == DNA ? 2: rht == PROTEIN ? 8: (rht == DNA2 || rht == DNAC) ? 1: 3;
            /*std::fprintf(stderr, "Shift %d\n", shift);*/
#define ITER do {\
            start += *spaces++;\
            if((nextc = lutptr[uint8_t(s_[start])]) == int8_t(-1)) {\
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
#define ITER do {start += *spaces++;\
            if((nextc = lutptr[uint8_t(s_[start])]) == int8_t(-1)) {\
                new_kmer = ENCODE_OVERFLOW;\
                goto rnk;\
            }\
            new_kmer = new_kmer * mul + nextc;\
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
    // ENCODE_OVERFLOW marks a k-mer with an invalid character, but when the
    // k-mer fills KmerT (poly-T at k = 32) it is also a valid encoding.
    INLINE bool kmer_is_valid(unsigned start) const {
        if(lutptr[uint8_t(s_[start])] == int8_t(-1)) return false;
        for(const auto space: sp_.s_)
            if(lutptr[uint8_t(s_[start += space])] == int8_t(-1)) return false;
        return true;
    }
    INLINE KmerT next_canonicalized_minimizer() {
        assert(has_next_kmer());
        const unsigned start = pos_++;
        KmerT nk = kmer(start);
        // Invalid k-mers stay out of the window, as in the unspaced rolling
        // encoders; canonicalizing the marker would turn it into poly-A.
        if(nk == ENCODE_OVERFLOW && !kmer_is_valid(start)) return ENCODE_OVERFLOW;
        if(rht == DNA) nk = canonical_representation(nk, sp_.k_);
        const KmerT kscore(scorer_(nk, getdata()));
        return qmap_.next_value(nk, kscore);
    }
    auto max_in_queue() const {return qmap_.max_in_queue();}

    bool canonicalize() const {return canonicalize_;}
    void canonicalize(bool value) {canonicalize_ = value;}

    auto pos() const {return pos_;}
    void pos(uint64_t v) {pos_ = v;}
    uint32_t k() const {return sp_.k_;}
    size_t n_in_queue() const {return qmap_.n_in_queue();}
};





// HashClass defaults to a polynomial hash modulo 2^61 - 1. A cyclic hash in a
// w-bit word gives positions i and i + w the same rotation, so k-mers longer
// than w collide structurally, for example when they differ by swapping two
// such positions or when their period divides w.
template<typename IntType, typename HashClass=PolyHash<IntType>>
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
        hashtype(enc);
        window(wsz);
        hasher_.seed(seed1, seed2);
        // Both strands must use the same hash function so that a k-mer and its
        // reverse complement map to the same value.
        rchasher_.seed(seed1, seed2);
        if(enc == PROTEIN_6_FRAME) throw NotImplementedError("Protein 6-frame not implemented.");
    }
    RollingHasher& operator=(const RollingHasher &o) {
        k_ = o.k_; canon_ = o.canon_; w_ = o.w_; seed1_ = o.seed1_; seed2_ = o.seed2_;
        hashtype(o.enctype_);
        window(o.w_);
        // The cyclic hashers depend on k, so rebuild them rather than keeping
        // the ones sized for this object's previous k.
        hasher_ = HashClass(k_, sizeof(IntType) * CHAR_BIT);
        rchasher_ = HashClass(k_, sizeof(IntType) * CHAR_BIT);
        hasher_.seed(seed1_, seed2_);
        rchasher_.seed(seed1_, seed2_);
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
            // Each position contributes one canonical hash, the smaller of its
            // two strand hashes, as in the unwindowed branch below.
            auto add_canon_hash = [&]() {
                const IntType v = std::min(hasher_.hashvalue, rchasher_.hashvalue);
                if((nextv = qmap_.next_value(v, lex_score(v))) != ENCODE_OVERFLOW)
                    func(nextv);
            };
            for(i = nf = 0; nf < k_ && i < l; ++i) {
                if((v1 = cstr_lut[uint8_t(s[i])]) == uint8_t(-1)) {
                    fixup_minimizer:
                    // Start a new k-mer after the invalid character.
                    nf = 0;
                    hasher_.reset();
                    rchasher_.reset();
                } // Fixme: this ignores both strands when one becomes 'N'-contaminated.
                  // In the future, encode the side that is still valid
                else hasher_.eat(v1), ++nf;
            }
            if(nf < k_) goto end; // All failed
            // Seed the reverse-strand hasher with the reverse complement of s[i - k_, i).
            for(size_t j = i; j-- > i - k_;) rchasher_.eat(cstr_rc_lut[uint8_t(s[j])]);
            add_canon_hash();
            for(;i < l; ++i) {
                if((v1 = cstr_lut[uint8_t(s[i])]) == uint8_t(-1))
                    goto fixup_minimizer;
                hasher_.update(cstr_lut[uint8_t(s[i - k_])], v1);
                rchasher_.reverse_update(cstr_rc_lut[uint8_t(s[i])], cstr_rc_lut[uint8_t(s[i - k_])]);
                add_canon_hash();
            }
            end:
            if(qmap_.partially_full())
                func(max_in_queue().el_);
        } else {
            for(i = nf = 0; nf < k_ && i < l; ++i) {
                if((v1 = cstr_lut[uint8_t(s[i])]) == uint8_t(-1)) {
                    fixup:
                    // Start a new k-mer after the invalid character.
                    nf = 0;
                    hasher_.reset();
                    rchasher_.reset();
                } // Fixme: this ignores both strands when one becomes 'N'-contaminated.
                  // In the future, encode the side that is still valid
                else hasher_.eat(v1), ++nf;
            }
            if(nf < k_) return; // All failed
            // Seed the reverse-strand hasher with the reverse complement of s[i - k_, i).
            for(size_t j = i; j-- > i - k_;) rchasher_.eat(cstr_rc_lut[uint8_t(s[j])]);
            func(std::min(hasher_.hashvalue, rchasher_.hashvalue));
            for(;i < l; ++i) {
                if((v1 = cstr_lut[uint8_t(s[i])]) == uint8_t(-1))
                    goto fixup;
                hasher_.update(cstr_lut[uint8_t(s[i - k_])], v1);
                rchasher_.reverse_update(cstr_rc_lut[uint8_t(s[i])], cstr_rc_lut[uint8_t(s[i - k_])]);
                func(std::min(hasher_.hashvalue, rchasher_.hashvalue));
            }
        }
    }

    template<typename Functor>
    void for_each_uncanon(const Functor &func, const char *s, size_t l) {
        // Clear the window first, so that a record too short for a k-mer
        // leaves it empty rather than holding the previous record's k-mers.
        qmap_.reset();
        if(l < size_t(k_)) return;
        hasher_.reset();
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
            if(unlikely((v1 = lutptr[uint8_t(s[i])]) == int8_t(-1))) {
                //std::fprintf(stderr, "Char %c/%d was missing... %d\n", s[i], s[i], lutptr[s[i]]);
                fixup:
                nf = 0; hasher_.reset();
            } else hasher_.eat(v1), ++nf;
        }
        if(nf < k_) return; // All failed
        use_val(hasher_.hashvalue);
        for(;i < l; ++i) {
            if(lutptr[uint8_t(s[i])] == int8_t(-1)) goto fixup;
            //auto ov = hasher_.hashvalue;
            hasher_.update(lutptr[uint8_t(s[i - k_])], lutptr[uint8_t(s[i])]);
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
        if(canon_) {
            for_each_canon<Functor>(func, s, l);
        } else       for_each_uncanon<Functor>(func, s, l);
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
    INLINE void for_each_uncanon(const Functor &func, gzFile fp, kseq_t *ks=nullptr) {
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
    auto max_in_queue() const {return qmap_.max_in_queue();}
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
        for(const auto k: c) {
            // The seeds are drawn in separate statements because the evaluation order of
            // function arguments is unspecified and differs between compilers and targets.
            const uint64_t seed1 = mt();
            const uint64_t seed2 = mt();
            hashers_.emplace_back(k, canon, enc, -1, seed1, seed2);
        }
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
            if((v1 = cstr_lut[uint8_t(s[i])]) == uint8_t(-1)) {
                fixup:
                if(i + 2 * mink >= l) return;
                i += mink;
                nf = 0;
                for(auto &h: hashers_) h.reset();
            } // Fixme: this ignores both strands when one becomes 'N'-contaminated.
              // In the future, encode the side that is still valid
            else {
                for(auto &h: hashers_) h.hasher_.eat(v1);
                ++nf;
            }
        }
        // Once a hasher's forward strand holds s[i - k_, i), seed its reverse-strand
        // hasher with the reverse complement of the same k-mer.
        auto seed_rc = [&](auto &h, size_t end) {
            for(size_t j = end; j-- > end - h.k_;) h.rchasher_.eat(cstr_rc_lut[uint8_t(s[j])]);
        };
        for(size_t hi = 0; hi < hashers_.size(); ++hi) {
            auto &h(hashers_[hi]);
            if(nf == h.k_) {
                seed_rc(h, i);
                func(std::min(h.hasher_.hashvalue, h.rchasher_.hashvalue), hi);
            }
        }
        for(;i < l; ++i) {
            if((v1 = cstr_lut[uint8_t(s[i])]) == uint8_t(-1))
                goto fixup;
            for(size_t hi = 0; hi < hashers_.size(); ++hi) {
                auto &h(hashers_[hi]);
                if(nf >= h.k_) {
                    h.rchasher_.reverse_update(cstr_rc_lut[uint8_t(s[i])], cstr_rc_lut[uint8_t(s[i - h.k_])]);
                    h.hasher_.update(cstr_lut[uint8_t(s[i - h.k_])], v1);
                    func(std::min(h.hasher_.hashvalue, h.rchasher_.hashvalue), hi);
                } else {
                    h.hasher_.eat(v1);
                    if(nf + 1 == h.k_) {
                        seed_rc(h, i + 1);
                        func(std::min(h.hasher_.hashvalue, h.rchasher_.hashvalue), hi);
                    }
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
    INLINE void for_each_hash(const Functor &func, gzFile fp, kseq_t *ks=nullptr) {
        if(canon_) for_each_canon<Functor>(func, fp, ks);
        else       for_each_uncanon<Functor>(func, fp, ks);
    }
    template<typename Functor>
    INLINE void for_each_hash(const Functor &func, const char *inpath, kseq_t *ks=nullptr) {
        gzFile fp = gzopen(inpath, "rb");
        if(!fp) throw file_open_error(inpath);
        gzbuffer(fp, 1<<18);
        for_each_hash<Functor>(func, fp, ks);
        gzclose(fp);
    }
    template<typename Functor>
    INLINE void for_each_uncanon(const Functor &func, gzFile fp, kseq_t *ks=nullptr) {
        bool destroy;
        if(ks == nullptr) ks = kseq_init(fp), destroy = true;
        else            kseq_assign(ks, fp), destroy = false;
        for_each_uncanon<Functor>(func, ks);
        if(destroy) kseq_destroy(ks);
    }
    template<typename Functor>
    INLINE void for_each_canon(const Functor &func, gzFile fp, kseq_t *ks=nullptr) {
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
    u64 submitted(0), todo(paths.size());
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


// How a database turns sequence into keys, beyond k/w/spacing: the alphabet
// (DNA, or an amino-acid alphabet such as PROTEIN20 / PROTEIN_14 / PROTEIN_3BIT /
// PROTEIN_6), whether nucleotide input is 6-frame translated first
// (Kraken2-style translated search; sixframe.h), and the hashed-tier key width
// (0 = the default). Build and classify must use the same KeySpec; compact
// databases record it.
struct KeySpec {
    InputType rht = DNA;
    bool translate = false;
    unsigned hash_bits = 0;
};

// KmerKeys<ScoreType, KeyT>: the k-mer -> table key used by the
// classification database (estimate, build, classify), chosen by the
// alphabet's capacity (rhtraits rh2n). With KeyT = u64:
//   EXACT64   k <= rh2n(alphabet, 8) (DNA: k <= 32): the packed k-mer itself
//             from Encoder<ScoreType, u64> (spaced seeds and minimizer
//             windows as before); key_bits = bit width of mul^k - 1 (2k for
//             DNA), so the compact table stores exactly that domain.
//   HASHED64  rh2n(alphabet, 8) < k <= rh2n(alphabet, 16) (DNA 32 < k <= 64):
//             the exact 128-bit k-mer from Encoder<ScoreType, u128> (spaced
//             seeds and minimizer windows supported, the same minimizers as
//             EXACT128) folded to 64 bits by a fixed mixer;
//             k > rh2n(alphabet, 16): a rolling hash (RollingHasher<uint64_t>:
//             a polynomial mod 2^61 - 1, bijectively mixed to 64 bits;
//             canonical for DNA; windows select the minimum hash), unspaced
//             only.
// With KeyT = u128 (compact databases):
//   EXACT128  rh2n(alphabet, 8) < k <= rh2n(alphabet, 16) (DNA: 32 < k <= 64):
//             the packed k-mer from Encoder<ScoreType, u128>; key_bits = bit
//             width of mul^k - 1.
//   HASHED128 k > rh2n(alphabet, 16): RollingHasher<u128> (two independent
//             61-bit lanes, each mixed to 64 bits), unspaced only. (Spaced
//             seeds beyond 128 bits would need a spaced rolling hash.)
// Hashed keys are the low hash_bits bits of the hash (u64: default/max 61;
// u128: default 122, range [65, 128]): a SPECIFICITY/SIZE tradeoff -- each bit
// is one more stored remainder bit per slot, and distinct k-mers merge
// through the hash with ~n^2 / 2^(hash_bits + 1) expected pairs over n
// distinct k-mers (at n = 1e9: 61 bits ~0.2 pairs, 48 bits ~1.8e3, 40 bits
// ~4.5e5). Lookups are exact on the key. Build and classification use the
// same KmerKeys, so they agree by construction.
template<typename ScoreType=score::Lex, typename KeyT=u64>
class KmerKeys {
    static_assert(std::is_same<KeyT, u64>::value || std::is_same<KeyT, u128>::value, "KmerKeys: KeyT is u64 or u128");
    static constexpr bool wide = std::is_same<KeyT, u128>::value;
public:
    enum Tier: uint8_t {EXACT64, HASHED64, EXACT128, HASHED128};
    using key_type = KeyT;
    static constexpr unsigned DEFAULT_HASH_BITS = wide ? 122: 61;
private:
    Tier tier_;
    Spacer sp_;
    bool canon_;
    void *data_;
    InputType rht_;
    std::unique_ptr<Encoder<ScoreType, KeyT>> enc_;
    std::unique_ptr<Encoder<ScoreType, u128>> enc128_;   // HASHED64 for k within 128 bits: fold the exact k-mer
    std::unique_ptr<RollingHasher<KeyT>> rh_;
    bool translate_ = false;
    // fixed 128 -> 64-bit fold (splitmix64 finalizer over both halves)
    static u64 mix64_(u64 z) {
        z = (z ^ (z >> 30)) * 0xbf58476d1ce4e5b9ull;
        z = (z ^ (z >> 27)) * 0x94d049bb133111ebull;
        return z ^ (z >> 31);
    }
    static u64 fold128_(u128 x) {return mix64_(u64(x) ^ mix64_(u64(x >> 64) + 0x9e3779b97f4a7c15ull));}
    unsigned hash_bits_;
    KeyT hmask_;
    static KeyT mask_of(unsigned b) {return b >= sizeof(KeyT) * CHAR_BIT ? ~KeyT(0): ((KeyT(1) << b) - 1);}
    bool exact() const {return tier_ == EXACT64 || tier_ == EXACT128;}
    void make() {
        if(rht_ != DNA) canon_ = false;   // reverse complements are a nucleotide notion
        if(translate_ && rht_ == DNA) RUNTIME_ERROR("KmerKeys: 6-frame translation needs an amino-acid alphabet");
        if(exact()) {
            enc_.reset(new Encoder<ScoreType, KeyT>(nullptr, 0, sp_, data_, canon_));
            if(rht_ != DNA) enc_->hashtype(rht_);
        } else if(!wide && sp_.k_ <= rh2n(rht_, 16)) {
            enc128_.reset(new Encoder<ScoreType, u128>(nullptr, 0, sp_, data_, canon_));
            if(rht_ != DNA) enc128_->hashtype(rht_);
        } else {
            if(!sp_.unspaced())
                RUNTIME_ERROR(ks::sprintf("k %u exceeds the %zu bases a %zu-bit key holds for alphabet %s; the hashed-key tier needs an unspaced seed.",
                                          sp_.k_, rh2n(rht_, sizeof(KeyT)), sizeof(KeyT) * CHAR_BIT, to_string(rht_).data()).data());
            rh_.reset(new RollingHasher<KeyT>(sp_.k_, canon_, rht_, sp_.unwindowed() ? -1: (long long)sp_.w_));
        }
    }
public:
    static Tier tier_for(unsigned k, InputType rht = DNA) {
        if constexpr(wide) return k <= rh2n(rht, 16) ? EXACT128: HASHED128;
        else return k <= rh2n(rht, 8) ? EXACT64: HASHED64;
    }
    // bits the keys span
    static unsigned key_bits_for(unsigned k, InputType rht = DNA, unsigned hash_bits = DEFAULT_HASH_BITS) {
        const Tier t = tier_for(k, rht);
        if(t == HASHED64 || t == HASHED128) return hash_bits;
        unsigned __int128 v = 1;
        for(unsigned i = 0; i < k; ++i) v *= mul(rht);   // k <= rh2n(rht, sizeof(KeyT)): no overflow
        --v;
        unsigned b = 1;
        while(b < 128 && (v >> b)) ++b;
        return b;
    }
    KmerKeys(const Spacer &sp, bool canon, void *data, const KeySpec &spec):
        KmerKeys(sp, canon, data, spec.rht, spec.hash_bits ? spec.hash_bits: DEFAULT_HASH_BITS, spec.translate) {}
    KmerKeys(const Spacer &sp, bool canon, void *data = nullptr, InputType rht = DNA, unsigned hash_bits = DEFAULT_HASH_BITS,
             bool translate = false):
        tier_(tier_for(sp.k_, rht)), sp_(sp), canon_(canon), data_(data), rht_(rht), hash_bits_(hash_bits), hmask_(mask_of(hash_bits)),
        translate_(translate) {
        if(wide ? (hash_bits < 65 || hash_bits > 128): (hash_bits < 16 || hash_bits > 61))
            RUNTIME_ERROR(std::string(wide ? "KmerKeys<u128>: hash_bits must lie in [65, 128]": "KmerKeys<u64>: hash_bits must lie in [16, 61]"));
        make();
    }
    KmerKeys(const KmerKeys &o): tier_(o.tier_), sp_(o.sp_), canon_(o.canon_), data_(o.data_), rht_(o.rht_),
        hash_bits_(o.hash_bits_), hmask_(o.hmask_), translate_(o.translate_) {make();}
    KmerKeys &operator=(const KmerKeys &o) {
        tier_ = o.tier_; sp_ = o.sp_; canon_ = o.canon_; data_ = o.data_; rht_ = o.rht_;
        hash_bits_ = o.hash_bits_; hmask_ = o.hmask_; translate_ = o.translate_; make();
        return *this;
    }
    bool translated() const {return translate_;}
    InputType alphabet() const {return rht_;}
    Tier tier() const {return tier_;}
    unsigned key_bits() const {return key_bits_for(sp_.k_, rht_, hash_bits_);}
    unsigned hash_bits() const {return hash_bits_;}
    const Spacer &spacer() const {return sp_;}
    unsigned span() const {return sp_.c_;}   // positions one key covers
    // keys of one sequence (residues in the key alphabet, or nucleotides
    // to translate when translated())
    template<typename Functor>
    void for_each(const Functor &func, const char *s, size_t l) {
        if(translate_) {
            sixframe::for_each_peptide(s, l, [&](const char *pep, size_t pl, int) {encode_(func, pep, pl);});
            return;
        }
        encode_(func, s, l);
    }
    template<typename Functor>
    void for_each(const Functor &func, const char *path, kseq_t *ks) {
        if(translate_) {   // per record: translate, then encode
            gzFile fp = gzopen(path, "rb");
            if(!fp) RUNTIME_ERROR(std::string("Could not open ") + path);
            const bool own = !ks;
            if(own) ks = kseq_init(fp);
            else kseq_assign(ks, fp);
            while(kseq_read(ks) >= 0) for_each(func, ks->seq.s, ks->seq.l);
            if(own) kseq_destroy(ks);
            gzclose(fp);
            return;
        }
        if(exact()) enc_->for_each([&](KeyT x) {func(x);}, path, ks);
        else if(enc128_) enc128_->for_each([&](u128 x) {func(KeyT(fold128_(x)) & hmask_);}, path, ks);
        else rh_->for_each_hash([&](KeyT h) {func(KeyT(h & hmask_));}, path, ks);
    }
private:
    template<typename Functor>
    void encode_(const Functor &func, const char *s, size_t l) {
        if(exact()) enc_->for_each([&](KeyT x) {func(x);}, s, l);
        else if(enc128_) enc128_->for_each([&](u128 x) {func(KeyT(fold128_(x)) & hmask_);}, s, l);
        else rh_->for_each_hash([&](KeyT h) {func(KeyT(h & hmask_));}, s, l);
    }
public:
};

template<typename ScoreType, typename SketchType>
void fill_lmers(SketchType &sketch, const std::string &path, const Spacer &space, bool canonicalize=true,
                void *data=nullptr, kseq_t *ks=nullptr) {
#if USE_HASH_FILLER
    detail::HashFiller<SketchType> hf(sketch);
#endif
    sketch.not_ready();
    KmerKeys<ScoreType> keys(space, canonicalize, data);   // same keys as the database build (k > 32: hashed)
#if USE_HASH_FILLER
    keys.for_each([&](u64 min) {hf.add(min);}, path.data(), ks);
#else
    keys.for_each([&](u64 min) {sketch.addh(min);}, path.data(), ks);
#endif
}


#if 0
#endif

} //namespace bns
#endif // _EMP_ENCODER_H__
