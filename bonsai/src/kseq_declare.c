#include "kseq_declare.h"
#include <stdio.h>

#ifdef __cplusplus
extern "C" {
#endif
static inline void trim_readno(kstring_t *s)
{
    if (s->l > 2 && s->s[s->l-2] == '/' && isdigit(s->s[s->l-1]))
        s->l -= 2, s->s[s->l] = 0;
}


bseq1_t *bseq_read(int chunk_size, int *n_, void *ks1_, void *ks2_)
{
    kseq_t *ks = (kseq_t*)ks1_, *ks2 = (kseq_t*)ks2_;
    int m, n, size;
    m = n = size = 0;
    bseq1_t *seqs = 0;
    while (kseq_read(ks) >= 0) {
        if (ks2 && kseq_read(ks2) < 0) { // the 2nd file has fewer reads
            fprintf(stderr, "[W::%s] the 2nd file has fewer sequences.\n", __func__);
            break;
        }
        if (n >= m) {
            m = m? m<<1 : 4096;
            seqs = realloc(seqs, m * sizeof(bseq1_t));
        }
        trim_readno(&ks->name);
        kseq2bseq1(ks, seqs + n);
        seqs[n].id = n;
        size += seqs[n++].l_seq;
        if (ks2) {
            trim_readno(&ks2->name);
            kseq2bseq1(ks2, seqs + n);
            seqs[n].id = n;
            size += seqs[n++].l_seq;
        }
        if (size >= chunk_size && (n&1) == 0) break;
    }
    if (size == 0) { // test if the 2nd file is finished
        if (ks2 && kseq_read(ks2) >= 0)
            fprintf(stderr, "[W::%s] the 1st file has fewer sequences.\n", __func__);
    }
    *n_ = n;
    return seqs;
}

#ifdef __cplusplus
}
#endif
