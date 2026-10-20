#pragma once
#include <array>
#include <atomic>
#include <cmath>
#include <fstream>
#include <unordered_map>
#include "kspp/ks.h"
#include "encoder.h"
#include "feature_min.h"
#include "kmermap.h"
#include "klib/kthread.h"
#include "util.h"
#include "uniq.h"

namespace bns {
using tax_counter = linear::counter<tax_t, u32>; // u32: long reads exceed 65535 hits

static void append_kraken_classification(const tax_counter &hit_counts,
                                  const std::vector<tax_t> &taxa,
                                  const tax_t taxon, const u32 ambig_count, const u32 missing_count,
                                  bseq1_t *bs, kstring_t *bks);
static void append_fastq_classification(const tax_counter &hit_counts,
                                 const std::vector<u32> &taxa,
                                 const tax_t taxon, const u32 ambig_count, const u32 missing_count,
                                 bseq1_t *bs, kstring_t *bks, const int verbose, const int is_paired);
static void append_taxa_runs(tax_t taxon, const std::vector<tax_t> &taxa, kstring_t *bks);


enum output_format: int {
    KRAKEN   = 1,
    FASTQ    = 2,
    EMIT_ALL = 4
};


INLINE void append_taxa_run(const tax_t last_taxa,
                            const u32 taxa_run,
                            ks::string &bks) {
    // U for unclassified (unambiguous but not in database)
    // A for ambiguous: ambiguous nucleotides
    // Actual taxon otherwise.
    switch(last_taxa) {
        case 0:            bks.putc_('U'); break;
        case (tax_t)-1:    bks.putc_('A'); break;
        default:           bks.putuw_(last_taxa); break;
    }
    bks.putc_(':'); bks.putuw_(taxa_run); bks.putc_('\t');
}


inline void append_taxa_runs(tax_t taxon, const std::vector<tax_t> &taxa, ks::string &bks) {
    if(taxon) {
        tax_t last_taxa(taxa[0]);
        unsigned taxa_run(1);
        for(unsigned i(1), end(taxa.size()); i != end; ++i) {
            if(taxa[i] == last_taxa) ++taxa_run;
            else {
                append_taxa_run(last_taxa, taxa_run, bks);
                last_taxa = taxa[i];
                taxa_run  = 1;
            }
        }
        append_taxa_run(last_taxa, taxa_run, bks); // Add last run.
        bks.back() = '\n'; // We add an extra tab. Replace  that with a newline.
        bks.terminate();
    } else bks.putsn("0:0\n", 4);
}

INLINE void append_counts(u32 count, const char character, ks::string &ks) {
    if(count) {
        char buf[] {character, ':'};
        ks.putsn_(buf, 2);
        ks.putuw_(count);
        ks.putc_('\t');
    }
}

// Lineage evidence of a read's call (Perseus-style channels, k-mer hits):
// exact (the call), ancestor (strict ancestors of the call), descendant (strict
// descendants) and conflict (unrelated taxa). hi: the high-confidence call.
struct ReadEvidence {
    u32 exact = 0, ancestor = 0, descendant = 0, conflict = 0;
    tax_t hi = 0;
    double identity = -1.;   // estimated identity to the database (I:), < 0: not computed
    bool emit_lineage = false, emit_hi = false;
};
// Fields, tab-terminated like M:/A:: E:exact,ancestor,descendant,conflict I:identity H:taxon
INLINE void append_evidence(const ReadEvidence *ev, ks::string &ks) {
    if(!ev) return;
    if(ev->emit_lineage) {
        ks.putsn_("E:", 2);
        ks.putuw_(ev->exact); ks.putc_(',');
        ks.putuw_(ev->ancestor); ks.putc_(',');
        ks.putuw_(ev->descendant); ks.putc_(',');
        ks.putuw_(ev->conflict); ks.putc_('\t');
        if(ev->identity >= 0.) {
            char buf[16];
            const int n = std::snprintf(buf, sizeof(buf), "I:%.4f\t", ev->identity);
            ks.putsn_(buf, n);
        }
    }
    if(ev->emit_hi) {
        ks.putsn_("H:", 2);
        ks.putuw_(ev->hi); ks.putc_('\t');
    }
}

inline void append_fastq_classification(const tax_counter &,
                                 const std::vector<tax_t> &taxa,
                                 const tax_t taxon, const u32 ambig_count, const u32 missing_count,
                                 bseq1_t *bs, ks::string &bks, const int verbose, const int is_paired,
                                 const ReadEvidence *ev = nullptr) {
    // Comment start/end as offsets: the appends below can realloc the buffer,
    // so pointers into it would dangle before the paired read copies them.
    bks.puts(bs->name);
    bks.putc_(' ');
    const size_t cms_off(bks.size());
    static const char lut[] {'C', 'U'};
    char tmp[] {lut[taxon == 0], '\t'};
    bks.putsn_(tmp, 2);
    bks.putuw_(taxon);
    bks.putc_('\t');
    bks.putl_(bs->l_seq);
    bks.putc_('\t');
    append_counts(missing_count, 'M', bks);
    append_counts(ambig_count,   'A', bks);
    append_evidence(ev, bks);
    if(verbose) append_taxa_runs(taxon, taxa, bks);
    else        bks.back() = '\n';
    const size_t cme_off(bks.size());
    // And now add the rest of the fastq record
    bks.putsn_(bs->seq, bs->l_seq);
    bks.putsn_("\n+\n", 3);
    bks.putsn_(bs->qual ? bs->qual: bs->seq, bs->l_seq); // Append sequence if it's a fasta record
    bks.putc_('\n');
    if(is_paired) {
        bks.puts((bs + 1)->name);
        bks.putc_(' ');
        bks.putsn_(bks.data() + cms_off, long(cme_off - cms_off)); // Add comment section in.
        bks.putc_('\n');
        bks.putsn_((bs + 1)->seq, (bs + 1)->l_seq);
        bks.putsn_("\n+\n", 3);
        bks.putsn_((bs + 1)->qual ? (bs + 1)->qual: (bs + 1)->seq, (bs + 1)->l_seq);
        bks.putc_('\n');
    }
    bks.terminate();
}



inline void append_kraken_classification(const tax_counter &,
                                  const std::vector<tax_t> &taxa,
                                  const tax_t taxon, const u32 ambig_count, const u32 missing_count,
                                  bseq1_t *bs, ks::string &bks, const ReadEvidence *ev = nullptr) {
    static const char tbl[]{'C', 'U'};
    bks.putc_(tbl[!taxon]);
    bks.putc_('\t');
    bks.puts(bs->name);
    bks.putc_('\t');
    bks.putuw_(taxon);
    bks.putc_('\t');
    bks.putw_(bs->l_seq);
    bks.putc_('\t');
    append_counts(missing_count, 'M', bks);
    append_counts(ambig_count,   'A', bks);
    append_evidence(ev, bks);
    append_taxa_runs(taxon, taxa, bks);
    bks.terminate();
}

// Taxonomic rank levels for identity-capped placement. 0: below species
// (strain, subspecies, serotype, ...), RANK_UNKNOWN: "no rank"/clade/absent.
enum: uint8_t {RANK_BELOW_SPECIES = 0, RANK_SPECIES = 1, RANK_GENUS = 2, RANK_FAMILY = 3, RANK_ORDER = 4,
               RANK_CLASS = 5, RANK_PHYLUM = 6, RANK_KINGDOM = 7, RANK_DOMAIN = 8, RANK_UNKNOWN = 255};
inline uint8_t rank_level(const std::string &r) {
    static const std::unordered_map<std::string, uint8_t> m {
        {"strain", RANK_BELOW_SPECIES}, {"subspecies", RANK_BELOW_SPECIES}, {"serotype", RANK_BELOW_SPECIES},
        {"serogroup", RANK_BELOW_SPECIES}, {"biotype", RANK_BELOW_SPECIES}, {"forma specialis", RANK_BELOW_SPECIES},
        {"varietas", RANK_BELOW_SPECIES}, {"isolate", RANK_BELOW_SPECIES}, {"genotype", RANK_BELOW_SPECIES},
        {"species", RANK_SPECIES}, {"species subgroup", RANK_SPECIES}, {"species group", RANK_SPECIES},
        {"subgenus", RANK_GENUS}, {"genus", RANK_GENUS}, {"subtribe", RANK_GENUS}, {"tribe", RANK_GENUS},
        {"subfamily", RANK_FAMILY}, {"family", RANK_FAMILY}, {"superfamily", RANK_FAMILY},
        {"suborder", RANK_ORDER}, {"order", RANK_ORDER}, {"superorder", RANK_ORDER},
        {"subclass", RANK_CLASS}, {"class", RANK_CLASS}, {"superclass", RANK_CLASS},
        {"subphylum", RANK_PHYLUM}, {"phylum", RANK_PHYLUM},
        {"kingdom", RANK_KINGDOM}, {"superkingdom", RANK_DOMAIN}, {"domain", RANK_DOMAIN}, {"realm", RANK_DOMAIN},
    };
    const auto it = m.find(r);
    return it == m.end() ? uint8_t(RANK_UNKNOWN): it->second;
}
using rank_map_t = std::unordered_map<tax_t, uint8_t>;
// Ranks from the third field of an NCBI nodes.dmp-style file
// ("id\t|\tparent\t|\trank\t|"); files without one give an empty map.
inline rank_map_t build_rank_map(const char *fn) {
    rank_map_t ret;
    std::ifstream is(fn);
    std::string line;
    while(std::getline(is, line)) {
        const size_t a = line.find('|');
        if(a == std::string::npos) continue;
        const size_t b = line.find('|', a + 1);
        if(b == std::string::npos) continue;
        const size_t c = line.find('|', b + 1);
        std::string r = line.substr(b + 1, (c == std::string::npos ? line.size(): c) - b - 1);
        const size_t f = r.find_first_not_of(" \t"), l = r.find_last_not_of(" \t");
        r = f == std::string::npos ? std::string(): r.substr(f, l - f + 1);
        ret[tax_t(std::strtoul(line.data(), nullptr, 10))] = rank_level(r);
    }
    return ret;
}
// Walks `taxon` up until it reaches a node of known rank >= cap; unknown-rank
// nodes are passed through. Stops at the root.
inline tax_t cap_rank(tax_t taxon, uint8_t cap, const khash_t(p) *taxmap, const rank_map_t &ranks) {
    for(unsigned guard = 0; taxon && guard < 4096; ++guard) {
        const auto it = ranks.find(taxon);
        if(it != ranks.end() && it->second != RANK_UNKNOWN && it->second >= cap) return taxon;
        const tax_t p = get_parent(taxmap, taxon);
        if(p == tax_t(-1) || p == 0) return taxon;
        taxon = p;
    }
    return taxon;
}

// DB: khash_t(c) (the original k-mer -> taxon hash) or a compact map
// (kmermap.h: CompactKmerMap with up-to-64-bit keys, CompactKmerMap128 with
// wider ones; looked up in per-read batches). khash and CompactKmerMap
// classify identically for the same keys.
template<typename DB> struct db_key_type {using type = u64;};
template<typename K> struct db_key_type<BasicCompactKmerMap<K>> {using type = K;};

template<typename ScoreType, typename DB=khash_t(c)>
struct ClassifierGeneric {
    using db_type = DB;
    using key_type = typename db_key_type<DB>::type;
    using keys_type = KmerKeys<ScoreType, key_type>;
    const DB *db_;
    const Spacer sp_;
    keys_type enc_;   // k-mer -> key (exact packed k-mer, or a rolling hash beyond the key type's capacity)
    uint32_t          nt_:16;
    uint32_t output_flag_:16;
    mutable std::atomic<u64> classified_[2];
    // Kraken2-style guards against divergent k-mers (defaults: off, the
    // original behavior). confidence_: after the best root-to-leaf path is
    // chosen, walk up until the clade under the call holds >= confidence_ of
    // the read's non-ambiguous k-mers (hits and misses); past the root ->
    // unclassified. min_hit_groups_: a call needs at least this many distinct
    // hitting keys (repeats of one minimizer count once).
    double confidence_ = 0.;
    unsigned min_hit_groups_ = 0;
    // confidence_hi_: a second, stricter call, the confidence walk-up of the
    // call at this threshold (H: field). emit_evidence_: per-read lineage
    // evidence of the call (E: field).
    // Confidence denominators: hits + misses (Kraken2, the default);
    // confidence_hits_only_: hits only (Perseus normalizes by classified
    // k-mers); confidence_quality_: the expected number of error-free k-mer
    // windows given the read's base qualities (expected_intact_windows), which
    // takes precedence and falls back to hits + misses for FASTA. Kraken2's
    // denominator counts error-destroyed k-mers against every clade, so noisy
    // reads fail any threshold; hits-only keeps them but also keeps reads whose
    // only hits are conserved genes of an out-of-database relative; the
    // quality-aware denominator keeps the first and drops the second.
    double confidence_hi_ = 0.;
    bool emit_evidence_ = false;
    bool confidence_hits_only_ = false;
    bool confidence_quality_ = false;
    // Identity-capped placement (ranks_ set): the read's identity to the
    // database is estimated Mash-style from its hits over the expected
    // error-free k-mers, (hits / expected)^(1/k) (expected from base qualities
    // when present, else hits + misses); below species_identity_ the call may be
    // no more specific than a family, below family_identity_ no more specific
    // than an order. Reads from a taxon missing from the database are placed
    // at the deepest rank their divergence supports instead of in a relative's
    // clade (or left unclassified). Defaults calibrated on held-out E. coli
    // reads simulated from real genomes (bench/real_crossfamily): database
    // species >= 0.994 (5th percentile), genus missing ~0.88, family missing
    // ~0.83 (medians).
    const rank_map_t *ranks_ = nullptr;
    double species_identity_ = 0.95, family_identity_ = 0.85;
    // uniq_: KrakenUniq-style per-taxon k-mer counters (uniq.h), one per
    // classification thread; null (the default) collects nothing.
    std::shared_ptr<UniqTable> uniq_;
    public:
    void enable_uniq(unsigned p = 12, UniqSketch kind = UniqSketch::HLL) {
        uniq_ = std::make_shared<UniqTable>(std::max<unsigned>(1, nt_), p, kind);
    }
    UniqTable *uniq() const {return uniq_.get();}
    void set_identity_placement(const rank_map_t *ranks, double species_identity = 0.95, double family_identity = 0.85) {
        if(!(family_identity <= species_identity && species_identity <= 1. && family_identity >= 0.))
            RUNTIME_ERROR("identity thresholds must satisfy 0 <= family <= species <= 1");
        ranks_ = ranks; species_identity_ = species_identity; family_identity_ = family_identity;
    }
    void set_confidence(double c) {
        if(!(c >= 0. && c <= 1.)) RUNTIME_ERROR("confidence must lie in [0, 1]");
        confidence_ = c;
    }
    void set_confidence_hi(double c) {
        if(!(c >= 0. && c <= 1.)) RUNTIME_ERROR("high confidence must lie in [0, 1]");
        confidence_hi_ = c;
    }
    void set_emit_evidence(bool v) {emit_evidence_ = v;}
    void set_confidence_hits_only(bool v) {confidence_hits_only_ = v;}
    void set_confidence_quality(bool v) {confidence_quality_ = v;}
    void set_min_hit_groups(unsigned n) {min_hit_groups_ = n;}
    double confidence() const {return confidence_;}
    unsigned min_hit_groups() const {return min_hit_groups_;}
    void set_emit_all(bool setting) {
        if(setting) output_flag_ |= output_format::EMIT_ALL;
        else        output_flag_ &= (~output_format::EMIT_ALL);
    }
    void set_emit_kraken(bool setting) {
        if(setting) output_flag_ |= output_format::KRAKEN;
        else        output_flag_ &= (~output_format::KRAKEN);
    }
    void set_emit_fastq(bool setting) {
        if(setting) output_flag_ |= output_format::FASTQ;
        else        output_flag_ &= (~output_format::FASTQ);
    }
    INLINE int get_emit_all()    const {return output_flag_ & output_format::EMIT_ALL;}
    INLINE int get_emit_kraken() const {return output_flag_ & output_format::KRAKEN;}
    INLINE int get_emit_fastq()  const {return output_flag_ & output_format::FASTQ;}
    // hash_bits: key width of hashed-tier keys (must match the build; a compact
    // database records it as its key width)
    ClassifierGeneric(const DB *map, const spvec_t &spaces, u8 k, std::uint16_t wsz, int num_threads=16,
                      bool emit_all=true, bool emit_fastq=true, bool emit_kraken=false, bool canonicalize=true,
                      unsigned hash_bits = keys_type::DEFAULT_HASH_BITS):
        ClassifierGeneric(map, spaces, k, wsz, num_threads, emit_all, emit_fastq, emit_kraken, canonicalize,
                          KeySpec{DNA, false, hash_bits}) {}
    // spec: the database's key encoding (alphabet, 6-frame translation, hash width)
    ClassifierGeneric(const DB *map, const spvec_t &spaces, u8 k, std::uint16_t wsz, int num_threads,
                      bool emit_all, bool emit_fastq, bool emit_kraken, bool canonicalize, const KeySpec &spec):
        db_(map),
        sp_(k, wsz, spaces),
        enc_(sp_, canonicalize, nullptr, spec),
        nt_(num_threads > 0 ? (uint16_t)(num_threads): (uint16_t)std::thread::hardware_concurrency()),
        output_flag_(0)
    {
        for(auto &c: classified_) c.store(0);
        set_emit_all(emit_all);
        set_emit_fastq(emit_fastq);
        set_emit_kraken(emit_kraken);
    }
    template<typename Q=DB, typename=typename std::enable_if<std::is_same<Q, khash_t(c)>::value>::type>
    ClassifierGeneric(const char *dbpath, const spvec_t &spaces, u8 k, std::uint16_t wsz, int num_threads=16,
                      bool emit_all=true, bool emit_fastq=true, bool emit_kraken=false, bool canonicalize=true):
        ClassifierGeneric(khash_load<khash_t(c)>(dbpath), spaces, k, wsz, num_threads, emit_all, emit_fastq, emit_kraken, canonicalize) {}
    u64 n_classified()   const {return classified_[0];}
    u64 n_unclassified() const {return classified_[1];}
};

INLINE void append_taxa_run(const tax_t last_taxa,
                            const u32 taxa_run,
                            kstring_t *bks) {
    // U for unclassified (unambiguous but not in database)
    // A for ambiguous: ambiguous nucleotides
    // Actual taxon otherwise.
    switch(last_taxa) {
        case 0:            kputc_('U', bks); break;
        case (tax_t)-1:    kputc_('A', bks); break;
        default:           kputuw_(last_taxa, bks); break;
    }

    kputc_(':', bks); kputuw_(taxa_run, bks); kputc_('\t', bks);
}


INLINE void append_counts(u32 count, const char character, kstring_t *ks) {
    if(count) {
        kputc_(character, ks);
        kputc_(':', ks);
        kputuw_(count, ks);
        kputc_('\t', ks);
    }
}

using Classifier = ClassifierGeneric<score::Lex>;
// Per-thread buffers for the batched (compact map) lookup path.
template<typename KeyT = u64>
struct LookupScratch {
    std::vector<KeyT> hits;          // distinct-hit counting (min_hit_groups)
    std::vector<KeyT> hitkeys;       // key of each hit, parallel to the read's taxa (uniq counting)
    std::vector<KeyT> uniq;          // consecutive-deduplicated keys of a mate
    std::vector<uint32_t> runs;      // run length of each uniq key
    std::vector<KeyT> kmers;
    std::vector<tax_t> taxa;
    std::vector<uint64_t> buf;
};

namespace {
template<typename C>
struct kt_data {
    const C &c_;
    const khash_t(p) *taxmap;
    bseq1_t *bs_;
    const unsigned per_set_;
    const unsigned total_;
    std::atomic<u64> &retstr_size_;
    const int is_paired_;
};
}

// Kraken2's confidence rule: from the call, walk up until the clade rooted at
// the current node holds >= confidence * total of the read's non-ambiguous
// k-mers (hits + misses); returns 0 (unclassified) if even the root falls short.
inline tax_t confident_ancestor(tax_t taxon, const tax_counter &hit_counts, const khash_t(p) *taxmap,
                                double total, double confidence) {
    const double need = confidence * total;
    auto in_clade = [&](tax_t t, tax_t root) {
        for(unsigned guard = 0; t && guard < 4096; ++guard) {
            if(t == root) return true;
            t = get_parent(taxmap, t);
            if(t == tax_t(-1)) return false;
        }
        return false;
    };
    while(taxon) {
        double clade = 0.;
        for(unsigned i = 0; i < hit_counts.size(); ++i)
            if(in_clade(hit_counts.keys()[i], taxon)) clade += hit_counts.vals()[i];
        if(clade >= need) return taxon;
        const tax_t p = get_parent(taxmap, taxon);
        taxon = (p == tax_t(-1)) ? 0: p;
    }
    return 0;
}

// Expected number of error-free span-windows of a read, from its Phred+33
// qualities: sum over windows of prod (1 - p_i). A base outside ACGT/U, or of
// quality 0, counts as certainly wrong (such windows are never emitted). The
// denominator of quality-aware confidence: a read from a genome in the
// database expects about this many hits whatever its error rate, so the
// confidence of the true clade stays near 1 for noisy reads, while an
// out-of-database relative (hits only in conserved genes) stays low.
inline double expected_intact_windows(const char *seq, const char *qual, size_t n, unsigned span) {
    static const auto lg = [] {
        std::array<double, 94> t{};
        t[0] = 0.;
        for(int q = 1; q < 94; ++q) t[q] = std::log1p(-std::pow(10., -q / 10.));
        return t;
    }();
    if(!span || n < span) return 0.;
    auto bad = [&](size_t i) {
        const int q = qual[i] - 33;
        switch(seq[i]) {case 'A': case 'C': case 'G': case 'T': case 'U': case 'a': case 'c': case 'g': case 't': case 'u': break; default: return true;}
        return q <= 0;
    };
    auto term = [&](size_t i) {return lg[std::min(qual[i] - 33, 93)];};
    double e = 0., s = 0.;
    unsigned zeros = 0;
    for(size_t i = 0; i < n; ++i) {
        if(bad(i)) ++zeros; else s += term(i);
        if(i >= span) {
            if(bad(i - span)) --zeros; else s -= term(i - span);
        }
        if(i + 1 >= span && !zeros) e += std::exp(s);
    }
    return e;
}

// Relation of hit taxon t to call x: 0 exact, 1 strict ancestor of x,
// 2 strict descendant of x, 3 unrelated (conflicting).
inline int lineage_relation(tax_t t, tax_t x, const khash_t(p) *taxmap) {
    if(t == x) return 0;
    auto descends = [&](tax_t a, tax_t root) {
        for(unsigned guard = 0; a && guard < 4096; ++guard) {
            a = get_parent(taxmap, a);
            if(a == tax_t(-1)) return false;
            if(a == root) return true;
        }
        return false;
    };
    if(descends(t, x)) return 2;
    if(descends(x, t)) return 1;
    return 3;
}
inline void lineage_evidence(tax_t call, const tax_counter &hit_counts, const khash_t(p) *taxmap, ReadEvidence &ev) {
    ev.exact = ev.ancestor = ev.descendant = ev.conflict = 0;
    if(!call) {
        for(unsigned i = 0; i < hit_counts.size(); ++i) ev.conflict += hit_counts.vals()[i];
        return;
    }
    for(unsigned i = 0; i < hit_counts.size(); ++i) {
        const u32 n = hit_counts.vals()[i];
        switch(lineage_relation(hit_counts.keys()[i], call, taxmap)) {
            case 0: ev.exact += n; break;
            case 1: ev.ancestor += n; break;
            case 2: ev.descendant += n; break;
            default: ev.conflict += n;
        }
    }
}
template<typename ScoreType, typename DB>
unsigned classify_seq(const ClassifierGeneric<ScoreType, DB> &c,
                      typename ClassifierGeneric<ScoreType, DB>::keys_type &enc,
                      const khash_t(p) *taxmap, bseq1_t *bs, const int is_paired, std::vector<tax_t> &taxa,
                      LookupScratch<typename ClassifierGeneric<ScoreType, DB>::key_type> &scratch,
                      UniqBatch *acc = nullptr) {
    LOG_DEBUG("starting classify_seq with bs at pointer = %p\n", static_cast<const void*>(bs));
    khiter_t ki;
    tax_counter hit_counts;
    u32 missing_count(0);
    tax_t taxon(0);
    ks::string bks(bs->sam, bs->l_sam);
    bks.clear();
    taxa.clear();

    // This simplification loses information about the run of congituous labels. Do these matter?
    // K-mers not emitted (ambiguous windows) account for the slack between the
    // number of positions and the number of emitted k-mers.
    const auto n_not_ambig = [&]() {return taxa.size() + missing_count;};
    const auto count_ambig = [&](u64 l, u64 prev_not_ambig) {
        if(enc.translated()) return u32(0);   // nucleotide positions do not map onto peptide k-mers
        const u64 np = l + 1 > enc.span() ? l - enc.span() + 1: 0; // positions holding a c_-span
        const u64 emitted = n_not_ambig() - prev_not_ambig;
        return u32(np > emitted ? np - emitted: 0);
    };
    // khash_t(c): one kh_get per k-mer as it is encoded (the original path).
    // CompactKmerMap: collect the mate's k-mers, look them all up in one
    // prefetched batch, then consume the results in the same order -- the
    // same hits, misses and taxa sequence as the streaming path.
    auto process = [&](const char *seq, size_t len) {
        if constexpr(std::is_same<DB, khash_t(c)>::value) {
            auto fn = [&] (u64 kmer) {
                //If the kmer is missing from our database, just say we don't know what it is.
                if((ki = kh_get(c, c.db_, kmer)) == kh_end(c.db_)) ++missing_count;
                else {
                    taxa.push_back(kh_val(c.db_, ki)), hit_counts.add(kh_val(c.db_, ki));
                    if(acc) scratch.hitkeys.push_back(kmer);
                    if(c.min_hit_groups_) scratch.hits.push_back(kmer);
                }
            };
            enc.for_each(fn, seq, len);
        } else {
            (void)ki;
            // Consecutive windows often share a minimizer: look each run of
            // equal keys up once (Kraken2 does the same), then expand the
            // results in order -- identical hits/misses/taxa to per-key lookups.
            scratch.uniq.clear();
            scratch.runs.clear();
            enc.for_each([&](auto kmer) {
                if(!scratch.uniq.empty() && scratch.uniq.back() == kmer) ++scratch.runs.back();
                else scratch.uniq.push_back(kmer), scratch.runs.push_back(1);
            }, seq, len);
            scratch.taxa.resize(scratch.uniq.size());
            c.db_->get_many(scratch.uniq.begin(), scratch.uniq.end(), scratch.taxa.begin(), scratch.buf);
            for(size_t i = 0; i < scratch.uniq.size(); ++i) {
                const tax_t t = scratch.taxa[i];
                if(t == DB::NOT_FOUND) {missing_count += scratch.runs[i]; continue;}
                for(uint32_t r = 0; r < scratch.runs[i]; ++r) taxa.push_back(t), hit_counts.add(t);
                if(acc) scratch.hitkeys.insert(scratch.hitkeys.end(), scratch.runs[i], scratch.uniq[i]);
                if(c.min_hit_groups_) scratch.hits.push_back(scratch.uniq[i]);
            }
        }
    };
    scratch.hits.clear();
    scratch.hitkeys.clear();
    process(bs->seq, bs->l_seq);
    unsigned ambig_count(count_ambig(bs->l_seq, 0));
    if(is_paired) {
        const auto pre(n_not_ambig());
        process((bs + 1)->seq, (bs + 1)->l_seq);
        ambig_count += count_ambig((bs + 1)->l_seq, pre);
    }

    taxon = resolve_tree(hit_counts, taxmap);
    if(taxon && c.min_hit_groups_) {
        std::sort(scratch.hits.begin(), scratch.hits.end());
        const size_t groups = std::unique(scratch.hits.begin(), scratch.hits.end()) - scratch.hits.begin();
        if(groups < c.min_hit_groups_) taxon = 0;
    }
    double total = double(taxa.size()) + (c.confidence_hits_only_ ? 0.: double(missing_count));
    if(c.confidence_quality_ && !enc.translated() && bs->qual) {
        total = expected_intact_windows(bs->seq, bs->qual, bs->l_seq, enc.span());
        if(is_paired && (bs + 1)->qual) total += expected_intact_windows((bs + 1)->seq, (bs + 1)->qual, (bs + 1)->l_seq, enc.span());
    }
    if(taxon && c.confidence_ > 0.) taxon = confident_ancestor(taxon, hit_counts, taxmap, total, c.confidence_);
    ReadEvidence ev;
    if(!enc.translated() && (c.ranks_ || c.emit_evidence_)) {
        double expected = double(taxa.size() + missing_count);
        if(bs->qual) {
            expected = expected_intact_windows(bs->seq, bs->qual, bs->l_seq, enc.span());
            if(is_paired && (bs + 1)->qual) expected += expected_intact_windows((bs + 1)->seq, (bs + 1)->qual, (bs + 1)->l_seq, enc.span());
        }
        const double ratio = expected > 0. ? std::min(1., double(taxa.size()) / expected): 0.;
        ev.identity = ratio > 0. ? std::pow(ratio, 1. / double(enc.k())): 0.;
        if(taxon && c.ranks_ && ev.identity < c.species_identity_)
            taxon = cap_rank(taxon, ev.identity < c.family_identity_ ? RANK_ORDER: RANK_FAMILY, taxmap, *c.ranks_);
    }
    // Two calls: the call above and its stricter confidence_hi_ walk-up
    // (an ancestor-or-self, or unclassified).
    if((ev.emit_hi = c.confidence_hi_ > 0.) && taxon) ev.hi = confident_ancestor(taxon, hit_counts, taxmap, total, c.confidence_hi_);
    if((ev.emit_lineage = c.emit_evidence_)) lineage_evidence(taxon, hit_counts, taxmap, ev);
    ++c.classified_[!taxon];
    if(acc) {   // hits to the taxon each k-mer maps to; the read to its call (0: unclassified)
        for(size_t i = 0; i < taxa.size(); ++i) acc->add_kmer(taxa[i], uniq_hash(scratch.hitkeys[i]));
        acc->add_read(taxon);
    }
    if(c.get_emit_all() || taxon) {
        switch(c.output_flag_) {
            case EMIT_ALL | FASTQ | KRAKEN: case FASTQ | KRAKEN: case FASTQ: case EMIT_ALL | FASTQ:
                append_fastq_classification(hit_counts, taxa, taxon, ambig_count, missing_count, bs, bks, c.get_emit_kraken(), is_paired, &ev); break;
            case EMIT_ALL | KRAKEN: case KRAKEN:
                append_kraken_classification(hit_counts, taxa, taxon, ambig_count, missing_count, bs, bks, &ev); break;
        }
    }
    LOG_DEBUG("About to return. Len of bks = %zu. len of string: %d\n", bks.size(), std::strlen(bks.data()));
    bs->l_sam = bks.size();
    bs->sam = bks.release();
    return bs->l_sam;
}


template<typename C>
inline void kt_for_helper(void *data_, long index, int tid) {
    kt_data<C> *data((kt_data<C> *)data_);
    size_t retstr_size(0);
    const int inc(!!data->is_paired_ + 1);
    typename C::keys_type enc(data->c_.enc_);
    std::vector<tax_t> taxa;
    LookupScratch<typename C::key_type> scratch;
    UniqBatch *acc = data->c_.uniq() ? &data->c_.uniq()->batch(tid): nullptr;
    // Each task covers [per_set * index, per_set * (index + 1)); the previous
    // (per_set + 1) * index bound left task 0 empty and dropped records.
    for(unsigned i(index * data->per_set_); i < std::min(data->per_set_ * (static_cast<unsigned>(index) + 1), data->total_); retstr_size += classify_seq(data->c_, enc, data->taxmap, data->bs_ + i, data->is_paired_, taxa, scratch, acc), i += inc);
    data->retstr_size_ += retstr_size;
}



using Classifier = ClassifierGeneric<score::Lex>;

template<typename C>
inline void classify_seqs(const C &c, const khash_t(p) *taxmap, bseq1_t *bs,
                          ks::string &cks, const unsigned chunk_size, const unsigned per_set, const int is_paired, ForPool &pool) {
    assert(per_set && ((per_set & (per_set - 1)) == 0));

    std::atomic<u64> retstr_size(0);
    kt_data<C> data{c, taxmap, bs, per_set, chunk_size, retstr_size, is_paired};
    pool.forpool(&kt_for_helper<C>, (void *)&data, chunk_size / per_set + 1);
    cks.resize(retstr_size.load());
    const int inc((is_paired != 0) + 1);
#if !NDEBUG
    for(u32 i(0); i < chunk_size; i += inc) {
        cks.putsn_(bs[i].sam, bs[i].l_sam);
        LOG_DEBUG("sam1: %s. len: %d. chunk size: %u\n", bs[i].sam, bs[i].l_sam, chunk_size);
    }
#else
    for(u32 i(0); i < chunk_size; cks.putsn_(bs[i].sam, bs[i].l_sam), i += inc);
#endif
    cks.terminate();
}

struct del_data {
    bseq1_t *seqs_;
    unsigned per_set_;
    unsigned total_;
};


template<typename C>
inline void process_dataset(const C &c, const khash_t(p) *taxmap, const char *fq1, const char *fq2,
                            std::FILE *out, unsigned chunk_size,
                            unsigned per_set) {
    // TODO: consider reusing buffers for processing large numbers of files.
    int nseq(0), max_nseq(0);
    gzFile ifp1(gzopen(fq1, "rb")), ifp2(fq2 ? gzopen(fq2, "rb"): nullptr);
    kseq_t *ks1(kseq_init(ifp1)), *ks2(ifp2 ? kseq_init(ifp2): nullptr);
    ks::string cks(256u);
    const int fn = fileno(out), is_paired(fq2 != 0);
    del_data dd{nullptr, per_set, chunk_size};
    ForPool pool(c.nt_);
    if((dd.seqs_ = bseq_read(chunk_size, &nseq, (void *)ks1, (void *)ks2)) == nullptr) {
        LOG_WARNING("Could not get any sequences from file, fyi.\n");
        goto fail; // Wheeeeee
    }
    classify_seqs(c, taxmap, dd.seqs_, cks, nseq, per_set, is_paired, pool);
    std::fprintf(stderr, "nseq: %i\n", nseq);
    max_nseq = std::max(max_nseq, nseq);
    while((dd.seqs_ = bseq_realloc_read(chunk_size, &nseq, (void *)ks1, (void *)ks2, dd.seqs_)) && nseq) {
        LOG_INFO("Read %i seqs with chunk size %u\n", nseq, chunk_size);
        max_nseq = std::max(max_nseq, nseq);
        // Classify
        classify_seqs(c, taxmap, dd.seqs_, cks, nseq, per_set, is_paired, pool);
        // Write out
        LOG_DEBUG("Emitting batch. str: %s", cks.data());
        if(cks.size() > (1ull << 16)) {
            cks.write(fn);
            cks.clear();
        }
    }
    cks.write(fn);
    cks.clear();
    // No use parallelizing the frees, there's a global lock on the freeing anyhow.
    for(int i(0); i < max_nseq; bseq_destroy(dd.seqs_ + i++));
    free(dd.seqs_);
    fail:
    // Clean up.
    gzclose(ifp1);
    kseq_destroy(ks1);
    if(ks2) kseq_destroy(ks2);
    if(ifp2) gzclose(ifp2);
}

static void append_fastq_classification(const tax_counter &,
                                        const std::vector<tax_t> &taxa,
                                        const tax_t taxon, const u32 ambig_count, const u32 missing_count,
                                        bseq1_t *bs, kstring_t *bks, const int verbose, const int is_paired) {
    // Comment start/end as offsets: kput* may realloc the buffer, so raw
    // pointers into it would dangle before the paired read copies them.
    kputs(bs->name, bks);
    kputc_(' ', bks);
    const size_t cms_off(bks->l);
    static const char lut[] {'C', 'U'};
    kputc_(lut[taxon == 0], bks);
    kputc_('\t', bks);
    kputuw_(taxon, bks);
    kputc_('\t', bks);
    kputl(bs->l_seq, bks);
    kputc_('\t', bks);
    append_counts(missing_count, 'M', bks);
    append_counts(ambig_count,   'A', bks);
    if(verbose) append_taxa_runs(taxon, taxa, bks);
    else        bks->s[bks->l - 1] = '\n';
    const size_t cme_off(bks->l);
    // And now add the rest of the fastq record
    kputsn_(bs->seq, bs->l_seq, bks);
    kputsn_("\n+\n", 3, bks);
    kputsn_(bs->qual, bs->l_seq, bks);
    kputc_('\n', bks);
    if(is_paired) {
        kputs((bs + 1)->name, bks);
        kputc_(' ', bks);
        kputsn_(bks->s + cms_off, (int)(cme_off - cms_off), bks); // Add comment section in.
        kputc_('\n', bks);
        kputsn_((bs + 1)->seq, (bs + 1)->l_seq, bks);
        kputsn_("\n+\n", 3, bks);
        kputsn_((bs + 1)->qual, (bs + 1)->l_seq, bks);
        kputc_('\n', bks);
    }
    bks->s[bks->l] = 0;
}

static void append_kraken_classification(const tax_counter &,
                                  const std::vector<tax_t> &taxa,
                                  const tax_t taxon, const u32 ambig_count, const u32 missing_count,
                                  bseq1_t *bs, kstring_t *bks) {
    static const char tbl[]{'C', 'U'};
    kputc_(tbl[!taxon], bks);
    kputc_('\t', bks);
    kputs(bs->name, bks);
    kputc_('\t', bks);
    kputuw_(taxon, bks);
    kputc_('\t', bks);
    kputw(bs->l_seq, bks);
    kputc_('\t', bks);
    append_counts(missing_count, 'M', bks);
    append_counts(ambig_count,   'A', bks);
    append_taxa_runs(taxon, taxa, bks);
    bks->s[bks->l] = '\0';
}

static void append_taxa_runs(tax_t taxon, const std::vector<tax_t> &taxa, kstring_t *bks) {
    if(taxon) {
        tax_t last_taxa(taxa[0]);
        unsigned taxa_run(1);
        for(unsigned i(1), end(taxa.size()); i != end; ++i) {
            if(taxa[i] == last_taxa) ++taxa_run;
            else {
                append_taxa_run(last_taxa, taxa_run, bks);
                last_taxa = taxa[i];
                taxa_run  = 1;
            }
        }
        append_taxa_run(last_taxa, taxa_run, bks); // Add last run.
        bks->s[bks->l - 1] = '\n'; // We add an extra tab. Replace  that with a newline.
        bks->s[bks->l] = '\0';
    } else kputsn("0:0\n", 4, bks);
}

} // namespace bns
