#pragma once
#include <atomic>
#include "kspp/ks.h"
#include "encoder.h"
#include "feature_min.h"
#include "kmermap.h"
#include "klib/kthread.h"
#include "util.h"

namespace bns {
using tax_counter = linear::counter<tax_t, u16>;

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

inline void append_fastq_classification(const tax_counter &,
                                 const std::vector<tax_t> &taxa,
                                 const tax_t taxon, const u32 ambig_count, const u32 missing_count,
                                 bseq1_t *bs, ks::string &bks, const int verbose, const int is_paired) {
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
                                  bseq1_t *bs, ks::string &bks) {
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
    append_taxa_runs(taxon, taxa, bks);
    bks.terminate();
}

// DB: khash_t(c) (the original k-mer -> taxon hash) or CompactKmerMap
// (kmermap.h; sketch's compact hash table, looked up in per-read batches).
// Both classify identically.
template<typename ScoreType, typename DB=khash_t(c)>
struct ClassifierGeneric {
    using db_type = DB;
    const DB *db_;
    const Spacer sp_;
    KmerKeys<ScoreType> enc_;   // k-mer -> key (exact 64-bit k-mer, or a rolling hash for k > 32)
    uint32_t          nt_:16;
    uint32_t output_flag_:16;
    mutable std::atomic<u64> classified_[2];
    public:
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
    ClassifierGeneric(const DB *map, const spvec_t &spaces, u8 k, std::uint16_t wsz, int num_threads=16,
                      bool emit_all=true, bool emit_fastq=true, bool emit_kraken=false, bool canonicalize=true):
        db_(map),
        sp_(k, wsz, spaces),
        enc_(sp_, canonicalize),
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
// Per-thread buffers for the batched (CompactKmerMap) lookup path.
struct LookupScratch {
    std::vector<u64> kmers;
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

template<typename ScoreType, typename DB>
unsigned classify_seq(const ClassifierGeneric<ScoreType, DB> &c,
                      KmerKeys<ScoreType> &enc,
                      const khash_t(p) *taxmap, bseq1_t *bs, const int is_paired, std::vector<tax_t> &taxa,
                      LookupScratch &scratch) {
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
                else taxa.push_back(kh_val(c.db_, ki)), hit_counts.add(kh_val(c.db_, ki));
            };
            enc.for_each(fn, seq, len);
        } else {
            (void)ki;
            scratch.kmers.clear();
            enc.for_each([&](u64 kmer) {scratch.kmers.push_back(kmer);}, seq, len);
            scratch.taxa.resize(scratch.kmers.size());
            c.db_->get_many(scratch.kmers.begin(), scratch.kmers.end(), scratch.taxa.begin(), scratch.buf);
            for(const tax_t t: scratch.taxa) {
                if(t == DB::NOT_FOUND) ++missing_count;
                else taxa.push_back(t), hit_counts.add(t);
            }
        }
    };
    process(bs->seq, bs->l_seq);
    unsigned ambig_count(count_ambig(bs->l_seq, 0));
    if(is_paired) {
        const auto pre(n_not_ambig());
        process((bs + 1)->seq, (bs + 1)->l_seq);
        ambig_count += count_ambig((bs + 1)->l_seq, pre);
    }

    ++c.classified_[!(taxon = resolve_tree(hit_counts, taxmap))];
    if(c.get_emit_all() || taxon) {
        switch(c.output_flag_) {
            case EMIT_ALL | FASTQ | KRAKEN: case FASTQ | KRAKEN: case FASTQ: case EMIT_ALL | FASTQ:
                append_fastq_classification(hit_counts, taxa, taxon, ambig_count, missing_count, bs, bks, c.get_emit_kraken(), is_paired); break;
            case EMIT_ALL | KRAKEN: case KRAKEN:
                append_kraken_classification(hit_counts, taxa, taxon, ambig_count, missing_count, bs, bks); break;
        }
    }
    LOG_DEBUG("About to return. Len of bks = %zu. len of string: %d\n", bks.size(), std::strlen(bks.data()));
    bs->l_sam = bks.size();
    bs->sam = bks.release();
    return bs->l_sam;
}


template<typename C>
inline void kt_for_helper(void *data_, long index, int) {
    kt_data<C> *data((kt_data<C> *)data_);
    size_t retstr_size(0);
    const int inc(!!data->is_paired_ + 1);
    KmerKeys<score::Lex> enc(data->c_.enc_);
    std::vector<tax_t> taxa;
    LookupScratch scratch;
    // Each task covers [per_set * index, per_set * (index + 1)); the previous
    // (per_set + 1) * index bound left task 0 empty and dropped records.
    for(unsigned i(index * data->per_set_); i < std::min(data->per_set_ * (static_cast<unsigned>(index) + 1), data->total_); retstr_size += classify_seq(data->c_, enc, data->taxmap, data->bs_ + i, data->is_paired_, taxa, scratch), i += inc);
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
