#include "classifier.h"

namespace kpg {


void kt_for_helper(void *data_, long index, int tid) {
    kt_data *data((kt_data *)data_);
    size_t retstr_size(0);
    const int inc(!!data->is_paired_ + 1);
    Encoder<lex_score> enc(data->c_.enc_);
    for(unsigned i(index * data->per_set_), end(std::min(i + data->per_set_, data->total_));
            i < end; i += inc) {
        retstr_size += classify_seq(data->c_, enc, data->taxmap, data->bs_ + i, data->is_paired_);
    }
    data->retstr_size_ += retstr_size;
}

void append_fastq_classification(const std::map<uint32_t, uint32_t> &hit_counts,
                                 const std::vector<uint32_t> &taxa,
                                 const uint32_t taxon, const uint32_t ambig_count, const uint32_t missing_count,
                                 bseq1_t *bs, kstring_t *bks, const int verbose, const int is_paired) {
    char *cms, *cme; // comment start, comment end -- used for using comment in both output reads.
    kputs(bs->name, bks);
    kputc(' ', bks);
    cms = bks->s + bks->l;
    kputc((taxon <= 0) * ('U' - 'C') + 'C', bks);
    kputc('\t', bks);
    kputuw(taxon, bks);
    kputc('\t', bks);
    kputl(bs->l_seq, bks);
    kputc('\t', bks);
    append_counts(missing_count, 'M', bks);
    append_counts(ambig_count,   'A', bks);
    if(verbose) append_taxa_runs(taxon, taxa, bks);
    else        bks->s[bks->l - 1] = '\n';
    cme = bks->s + bks->l;
    // And now add the rest of the fastq record
    kputsn(bs->seq, bs->l_seq, bks);
    kputsn("\n+\n", 3, bks);
    kputsn(bs->qual, bs->l_seq, bks);
    kputc('\n', bks);
    if(is_paired) {
        kputs((bs + 1)->name, bks);
        kputc(' ', bks);
        kputsn(cms, (int)(cme - cms), bks); // Add comment section in.
        kputc('\n', bks);
        kputsn((bs + 1)->seq, (bs + 1)->l_seq, bks);
        kputsn("\n+\n", 3, bks);
        kputsn((bs + 1)->qual, (bs + 1)->l_seq, bks);
        kputc('\n', bks);
    }
}

void append_kraken_classification(const std::map<uint32_t, uint32_t> &hit_counts,
                                  const std::vector<uint32_t> &taxa,
                                  const uint32_t taxon, const uint32_t ambig_count, const uint32_t missing_count,
                                  bseq1_t *bs, kstring_t *bks) {
    kputc((taxon < 1) * ('U' - 'C') + 'C', bks);
    //kputc(taxon ? 'C': 'U', bks); Equivalent to above but without a branch.
    //Probably meaningless in practice, but....
    kputc('\t', bks);
    kputs(bs->name, bks);
    kputc('\t', bks);
    kputuw(taxon, bks);
    kputc('\t', bks);
    kputw(bs->l_seq, bks);
    kputc('\t', bks);
    append_counts(missing_count, 'M', bks);
    append_counts(ambig_count,   'A', bks);
    append_taxa_runs(taxon, taxa, bks);
}

void append_taxa_runs(uint32_t taxon, const std::vector<uint32_t> &taxa, kstring_t *bks) {
    if(taxon) {
        uint32_t last_taxa(taxa[0]), taxa_run(1);
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
    } else kputsn("0:0\n", 4, bks);
}

void kt_del_helper(void *data_, long index, int tid) {
    del_data *data((del_data *)data_);
    for(unsigned i(index * data->per_set_), end(std::min(i + data->per_set_, data->total_));
            i < end; ++i)
        bseq_destroy(data->seqs_ + i);
}

} // namespace kpg
