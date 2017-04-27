#ifndef _DB_H__
#define _DB_H__

#include <atomic>
#include <zlib.h>
#include "encoder.h"
#include "feature_min.h"
#include "klib/kthread.h"
#include "lib/ks.h"
#include "lib/util.h"

namespace emp {
void append_kraken_classification(const std::map<tax_t, std::uint32_t> &hit_counts,
                                  const std::vector<tax_t> &taxa,
                                  const tax_t taxon, const std::uint32_t ambig_count, const std::uint32_t missing_count,
                                  bseq1_t *bs, kstring_t *bks);
void append_fastq_classification(const std::map<tax_t, std::uint32_t> &hit_counts,
                                 const std::vector<std::uint32_t> &taxa,
                                 const tax_t taxon, const std::uint32_t ambig_count, const std::uint32_t missing_count,
                                 bseq1_t *bs, kstring_t *bks, const int verbose, const int is_paired);
void append_taxa_runs(tax_t taxon, const std::vector<tax_t> &taxa, kstring_t *bks);

enum output_format {
    KRAKEN   = 1,
    FASTQ    = 2,
    EMIT_ALL = 4
};

template<std::uint64_t (*score)(std::uint64_t, void *)=lex_score>
struct ClassifierGeneric {
    khash_t(c) *db_;
    Spacer sp_;
    Encoder<score> enc_;
    int nt_;
    int output_flag_;
    std::atomic<std::uint64_t> n_classified_;
    std::atomic<std::uint64_t> n_unclassified_;
    public:
    void set_emit_all(bool setting);
    void set_emit_kraken(bool setting);
    void set_emit_fastq(bool setting);
    INLINE int get_emit_all()    {return output_flag_ & output_format::EMIT_ALL;}
    INLINE int get_emit_kraken() {return output_flag_ & output_format::KRAKEN;}
    INLINE int get_emit_fastq()  {return output_flag_ & output_format::FASTQ;}
    ClassifierGeneric(khash_t(c) *map, spvec_t &spaces, std::uint8_t k, std::uint16_t wsz, int num_threads=16,
                      bool emit_all=true, bool emit_fastq=true, bool emit_kraken=false):
        db_(map),
        sp_(k, wsz, spaces),
        enc_(sp_),
        nt_(num_threads > 0 ? num_threads: 16),
        n_classified_(0),
        n_unclassified_(0)
    {
        set_emit_all(emit_all);
        set_emit_fastq(emit_fastq);
        set_emit_kraken(emit_kraken);
    }
    ClassifierGeneric(const char *dbpath, spvec_t &spaces, std::uint8_t k, std::uint16_t wsz, int num_threads=16,
                      bool emit_all=true, bool emit_fastq=true, bool emit_kraken=false):
        ClassifierGeneric(khash_load<khash_t(c)>(dbpath), spaces, k, wsz, num_threads, emit_all, emit_fastq, emit_kraken) {
    }
    ~ClassifierGeneric() {
    }
};

template<std::uint64_t (*score)(std::uint64_t, void *)>
void ClassifierGeneric<score>::set_emit_all(bool setting) {
    if(setting) output_flag_ |= output_format::EMIT_ALL;
    else        output_flag_ &= (~output_format::EMIT_ALL);
}

template<std::uint64_t (*score)(std::uint64_t, void *)>
void ClassifierGeneric<score>::set_emit_kraken(bool setting) {
    if(setting) output_flag_ |= output_format::KRAKEN;
    else        output_flag_ &= (~output_format::KRAKEN);
}

template<std::uint64_t (*score)(std::uint64_t, void *)>
void ClassifierGeneric<score>::set_emit_fastq(bool setting) {
    if(setting) output_flag_ |= output_format::FASTQ;
    else        output_flag_ &= (~output_format::FASTQ);
}

INLINE void append_taxa_run(const tax_t last_taxa,
                            const std::uint32_t taxa_run,
                            kstring_t *bks) {
    // U for unclassified (unambiguous but not in database)
    // A for ambiguous: ambiguous nucleotides
    // Actual taxon otherwise.
    switch(last_taxa) {
        case 0:            kputc_('U', bks); break;
        case (tax_t)-1:    kputc_('A', bks); break;
        default:           kputuw_(last_taxa, bks); break;
    }

    kputc_(':', bks); kputuw_(taxa_run, bks); kputc('\t', bks);
}


INLINE void append_counts(std::uint32_t count, const char character, kstring_t *ks) {
    if(count) {
        kputc_(character, ks);
        kputc_(':', ks);
        kputuw_(count, ks);
        kputc('\t', ks);
    }
}

using Classifier = ClassifierGeneric<lex_score>;
template<std::uint64_t (*score)(std::uint64_t, void *)>
unsigned classify_seq(ClassifierGeneric<score> &c, Encoder<score> &enc, khash_t(p) *taxmap, bseq1_t *bs, const int is_paired) {
    std::uint64_t kmer;
    khiter_t ki;
    std::map<tax_t, std::uint32_t> hit_counts;
    std::map<tax_t, std::uint32_t>::iterator it;
    std::uint32_t ambig_count(0), missing_count(0);
    tax_t taxon(0);
    ks::KString bks(bs->sam);
    const std::size_t reslen(bs->l_seq - c.sp_.c_ + 1);
    std::vector<tax_t> taxa;
    taxa.reserve(reslen > 0 ? reslen: 0); // Reserve memory without initializing

    enc.assign(bs->seq, bs->l_seq);
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
                if((it = hit_counts.lower_bound(kh_val(c.db_, ki))) == hit_counts.end() ||
                        it->first != kh_val(c.db_, ki)) {
                    hit_counts.emplace_hint(it, kh_val(c.db_, ki), 1);
                } else ++it->second;
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
                    if((it = hit_counts.lower_bound(kh_val(c.db_, ki))) == hit_counts.end() ||
                            it->first != kh_val(c.db_, ki))
                        hit_counts.emplace_hint(it, kh_val(c.db_, ki), 1);
                    else ++it->second;
                    taxa.push_back(0);
                }
            }
        }
    }

    if((taxon = resolve_tree(hit_counts, taxmap))) ++c.n_classified_;
    else                                           ++c.n_unclassified_;
    if(c.get_emit_all() || taxon) {
        if(c.get_emit_fastq()) {
            append_fastq_classification(hit_counts, taxa, taxon, ambig_count, missing_count, bs, bks, c.get_emit_kraken(), is_paired);
        } else if(c.get_emit_kraken()) {
            append_kraken_classification(hit_counts, taxa, taxon, ambig_count, missing_count, bs, bks);
        }
        // Else just increase quantitation... But that requires some kind of lock-free hash table for counting.
        // Use Jellyfish?
    }
    bs->sam   = bks.release();
    return bs->l_sam = bks.size();
}

namespace {
struct kt_data {
    ClassifierGeneric<lex_score> &c_;
    khash_t(p) *taxmap;
    bseq1_t *bs_;
    const unsigned per_set_;
    const unsigned total_;
    std::atomic<std::uint64_t> &retstr_size_;
    const int is_paired_;
};
}
void kt_for_helper(void *data_, long index, int tid);

inline void classify_seqs(Classifier &c, khash_t(p) *taxmap, bseq1_t *bs,
                          kstring_t *cks, const unsigned chunk_size, const unsigned per_set, const int is_paired) {
    assert(per_set && ((per_set & (per_set - 1)) == 0));

    std::atomic<std::uint64_t> retstr_size(0);
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

void kt_del_helper(void *data_, long index, int tid);

#if !NDEBUG
#define DBKS(ks) do {\
    LOG_DEBUG("String: %s. Max: %zu. Len: %zu.\n", (ks) ? (ks)->s ? (ks)->s: "Unset string": "nullptr", (ks)->m, (ks)->l);\
    } while(0)
#endif

inline void process_dataset(Classifier &c, khash_t(p) *taxmap, const char *fq1, const char *fq2,
                     std::FILE *out, unsigned chunk_size,
                     unsigned per_set) {
    int nseq(0);
    gzFile ifp1(gzopen(fq1, "rb")), ifp2(fq2 ? gzopen(fq2, "rb"): nullptr);
    kseq_t *ks1(kseq_init(ifp1)), *ks2(ifp2 ? kseq_init(ifp2): nullptr);
    del_data dd{nullptr, per_set, chunk_size};
    ks::KString cks(256u);
    const int is_paired(!!fq2);
    while((dd.seqs_ = bseq_read(chunk_size, &nseq, (void *)ks1, (void *)ks2))) {
        LOG_INFO("Read %i seqs with chunk size %u\n", nseq, chunk_size);
        // Classify
        classify_seqs(c, taxmap, dd.seqs_, cks, nseq, per_set, is_paired);
        // Write out
        std::fwrite(cks.data(), 1, cks.size(), out);
        // Delete
        kt_for(c.nt_, &kt_del_helper, (void *)&dd, nseq / per_set + 1);
        cks.clear();
        free(dd.seqs_);
    }
    // Clean up.
    gzclose(ifp1); gzclose(ifp2);
    kseq_destroy(ks1); kseq_destroy(ks2);
}


} // namespace emp

#endif // #ifndef _DB_H__

