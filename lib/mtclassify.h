#ifndef MT_CLASSIFY_H
#define MT_CLASSIFY_H

#include "lib/db.h"
#include "klib/kthread.h"

namespace emp {
void kt_for_helper(void *data_, long index, int tid);

namespace {
struct kt_data {
    Classifier &c_;
    khash_t(p) *taxmap;
    bseq1_t *bs_;
    unsigned per_set_;
    unsigned total_;
    std::atomic<uint64_t> &retstr_size_;
};
}

INLINE void classify_seqs(Classifier &c, khash_t(p) *taxmap, bseq1_t *bs,
                          unsigned chunk_size, kstring_t *cks, unsigned per_set=32u) {
    // Uses bitmath to convert a power of two into its log2
    // Produces bogus results on other numbers!
    static const unsigned ilog2ps(__builtin_ctz(per_set));
    assert(per_set && ((per_set & (per_set - 1)) == 0));

    std::atomic<uint64_t> retstr_size(0);
    kt_data data{c, taxmap, bs, per_set, chunk_size, retstr_size};
    kt_for(c.nt_, &kt_for_helper, (void *)&data, (chunk_size >> ilog2ps) + 1);
    ks_resize(cks, retstr_size);
    for(unsigned i(0); i < chunk_size; ++i) kputsn(bs[i].sam, bs[i].l_sam, cks);
}


void kt_for_helper(void *data_, long index, int tid) {
    kt_data *data((kt_data *)data_);
    size_t retstr_size(0);
    for(unsigned i(index * data->per_set_), end(std::min(i + data->per_set_, data->total_));
            i < end; ++i)
        retstr_size += classify_seq(data->c_, data->taxmap, data->bs_ + i);
    data->retstr_size_ += retstr_size;
}

} // namespace emp

#endif
