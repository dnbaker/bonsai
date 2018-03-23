#ifndef _DB_H__
#define _DB_H__
#include <atomic>
#include "kspp/ks.h"
#include "encoder.h"
#include "feature_min.h"
#include "klib/kthread.h"
#include "util.h"

namespace emp {
using tax_counter = linear::counter<tax_t, u16>;

void append_kraken_classification(const tax_counter &hit_counts,
                                  const std::vector<tax_t> &taxa,
                                  const tax_t taxon, const u32 ambig_count, const u32 missing_count,
                                  bseq1_t *bs, kstring_t *bks);
void append_fastq_classification(const tax_counter &hit_counts,
                                 const std::vector<u32> &taxa,
                                 const tax_t taxon, const u32 ambig_count, const u32 missing_count,
                                 bseq1_t *bs, kstring_t *bks, const int verbose, const int is_paired);
void append_taxa_runs(tax_t taxon, const std::vector<tax_t> &taxa, kstring_t *bks);

enum output_format {
    KRAKEN   = 1,
    FASTQ    = 2,
    EMIT_ALL = 4
};

template<typename ScoreType>
struct ClassifierGeneric {
    khash_t(c) *db_;
    Spacer sp_;
    Encoder<ScoreType> enc_;
    int nt_;
    int output_flag_;
    std::atomic<u64> classified_[2];
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
    INLINE int get_emit_all()    {return output_flag_ & output_format::EMIT_ALL;}
    INLINE int get_emit_kraken() {return output_flag_ & output_format::KRAKEN;}
    INLINE int get_emit_fastq()  {return output_flag_ & output_format::FASTQ;}
    ClassifierGeneric(khash_t(c) *map, spvec_t &spaces, u8 k, std::uint16_t wsz, int num_threads=16,
                      bool emit_all=true, bool emit_fastq=true, bool emit_kraken=false, bool canonicalize=true):
        db_(map),
        sp_(k, wsz, spaces),
        enc_(sp_, canonicalize),
        nt_(num_threads > 0 ? num_threads: 16),
        classified_{0, 0}
    {
        set_emit_all(emit_all);
        set_emit_fastq(emit_fastq);
        set_emit_kraken(emit_kraken);
    }
    ClassifierGeneric(const char *dbpath, spvec_t &spaces, u8 k, std::uint16_t wsz, int num_threads=16,
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
template<typename ScoreType>
unsigned classify_seq(ClassifierGeneric<ScoreType> &c,
                      Encoder<ScoreType> &enc,
                      const khash_t(p) *taxmap, bseq1_t *bs, const int is_paired, std::vector<tax_t> &taxa) {
    u64 kmer;
    khiter_t ki;
    linear::counter<tax_t, u16> hit_counts;
    u32 ambig_count(0), missing_count(0);
    tax_t taxon(0);
    ks::string bks(bs->sam, bs->l_sam);
    bks.clear();
    taxa.clear();

    enc.assign(bs->seq, bs->l_seq);
    while(enc.has_next_kmer()) {
        if((kmer = enc.next_kmer()) == BF) ++ambig_count, taxa.push_back((tax_t)-1);
        // If the kmer is ambiguous, ignore it and move on.
        else {
            if((ki = kh_get(c, c.db_, kmer)) == kh_end(c.db_)) ++missing_count, taxa.push_back(0);
            //If the kmer is missing from our database, just say we don't know what it is.
            else {
                // Otherwise, increment the count.
                taxa.push_back(kh_val(c.db_, ki));
                hit_counts.add(kh_val(c.db_, ki));
            }
        }
    }
    if(is_paired) {
        enc.assign((bs + 1)->seq, (bs + 1)->l_seq);
        while(enc.has_next_kmer()) {
            if((kmer = enc.next_kmer()) == BF) ++ambig_count, taxa.push_back((tax_t)-1);
            // If the kmer is ambiguous, ignore it and move on.
            else {
                if((ki = kh_get(c, c.db_, kmer)) == kh_end(c.db_)) ++missing_count, taxa.push_back(0);
                //If the kmer is missing from our database, just say we don't know what it is.
                else {
                    // Check map for presence of the hit.
                    // If it's not there, insert it with a hint.
                    // Otherwise, increment the count.
                    taxa.push_back(kh_val(c.db_, ki));
                    hit_counts.add(kh_val(c.db_, ki));
                }
            }
        }
    }

    ++c.classified_[!(taxon = resolve_tree(hit_counts, taxmap))];
    if(c.get_emit_all() || taxon) {
        if(c.get_emit_fastq()) {
            append_fastq_classification(hit_counts, taxa, taxon, ambig_count, missing_count, bs, kspp2ks(bks), c.get_emit_kraken(), is_paired);
        } else if(c.get_emit_kraken()) {
            append_kraken_classification(hit_counts, taxa, taxon, ambig_count, missing_count, bs, kspp2ks(bks));
        }
    }
    bs->sam   = bks.release();
    return bs->l_sam = bks.size();
}

namespace {
struct kt_data {
    ClassifierGeneric<score::Lex> &c_;
    khash_t(p) *taxmap;
    bseq1_t *bs_;
    const unsigned per_set_;
    const unsigned total_;
    std::atomic<u64> &retstr_size_;
    const int is_paired_;
};
}
void kt_for_helper(void *data_, long index, int tid);

inline void classify_seqs(Classifier &c, khash_t(p) *taxmap, bseq1_t *bs,
                          kstring_t *cks, const unsigned chunk_size, const unsigned per_set, const int is_paired) {
    assert(per_set && ((per_set & (per_set - 1)) == 0));

    std::atomic<u64> retstr_size(0);
    kt_data data{c, taxmap, bs, per_set, chunk_size, retstr_size, is_paired};
    kt_for(c.nt_, &kt_for_helper, (void *)&data, chunk_size / per_set + 1);
    ks_resize(cks, retstr_size.load());
    const int inc(!!is_paired + 1);
    for(unsigned i(0); i < chunk_size; i += inc) kputsn_(bs[i].sam, bs[i].l_sam, cks);
    cks->s[cks->l] = 0;
}

struct del_data {
    bseq1_t *seqs_;
    unsigned per_set_;
    unsigned total_;
};

#if !NDEBUG
#define DBKS(ks) do {\
    LOG_DEBUG("String: %s. Max: %zu. Len: %zu.\n", (ks) ? (ks)->s ? (ks)->s: "Unset string": "nullptr", (ks)->m, (ks)->l);\
    } while(0)
#endif

inline void process_dataset(Classifier &c, khash_t(p) *taxmap, const char *fq1, const char *fq2,
                     std::FILE *out, unsigned chunk_size,
                     unsigned per_set) {
    // TODO: consider reusing buffers for processing large numbers of files.
    int nseq(0), max_nseq(0);
    gzFile ifp1(gzopen(fq1, "rb")), ifp2(fq2 ? gzopen(fq2, "rb"): nullptr);
    kseq_t *ks1(kseq_init(ifp1)), *ks2(ifp2 ? kseq_init(ifp2): nullptr);
    del_data dd{nullptr, per_set, chunk_size};
    ks::string cks(256u);
    const int fn = fileno(out);
    const int is_paired(fq2 != 0);
    if((dd.seqs_ = bseq_read(chunk_size, &nseq, (void *)ks1, (void *)ks2)) == nullptr) {
        LOG_WARNING("Could not get any sequences from file, fyi.\n");
        goto fail; // Wheeeeee
    }
    max_nseq = std::max(max_nseq, nseq);
    while((dd.seqs_ = bseq_realloc_read(chunk_size, &nseq, (void *)ks1, (void *)ks2, dd.seqs_))) {
        LOG_INFO("Read %i seqs with chunk size %u\n", nseq, chunk_size);
        max_nseq = std::max(max_nseq, nseq);
        // Classify
        classify_seqs(c, taxmap, dd.seqs_, kspp2ks(cks), nseq, per_set, is_paired);
        // Write out
        cks.write(fn);
        cks.clear();
    }
    // No use parallelizing the frees, there's a global lock on the freeing anyhow.
    for(int i(0); i < max_nseq; bseq_destroy(dd.seqs_ + i++));
    free(dd.seqs_);
    fail:
    // Clean up.
    gzclose(ifp1); gzclose(ifp2);
    kseq_destroy(ks1); kseq_destroy(ks2);
}


} // namespace emp

#endif // #ifndef _DB_H__

