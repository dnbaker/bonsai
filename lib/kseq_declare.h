#ifndef _KSEQ_DECLARE_
#define _KSEQ_DECLARE_
#include <zlib.h>
#include "klib/kseq.h"
#include <stdlib.h>
#include <stdint.h>

KSEQ_INIT(gzFile, gzread)

#ifndef INLINE
#  if __GNUC__ || __clang__
#  define INLINE __attribute__((always_inline)) inline
#  else
#  define INLINE inline
#  endif
#endif

#ifdef __cplusplus
extern "C" {
#endif

typedef struct {
    int l_seq, id;
    unsigned l_sam;
    char *name, *comment, *seq, *qual, *sam;
} bseq1_t;

static INLINE char *ksdup(const kstring_t *ks) {
    if(ks->l == 0) return NULL;
    char *ret = (char *)malloc(ks->l + 1);
    memcpy(ret, ks->s, ks->l + 1);
    return ret;
}

// From bwa -- batch parsing.
static inline void kseq2bseq1(const kseq_t *ks, bseq1_t *s)
{ // TODO: it would be better to allocate one chunk of memory, but probably it does not matter in practice
    s->name    = ksdup(&ks->name);
    s->comment = ksdup(&ks->comment);
    s->seq     = ksdup(&ks->seq);
    s->qual    = ksdup(&ks->qual);
    s->l_seq   = ks->seq.l;
    s->sam     = NULL;
    s->l_sam   = 0;
    s->id      = 0;
}

bseq1_t *bseq_read(int chunk_size, int *n_, void *ks1_, void *ks2_);

static inline void bseq_destroy(bseq1_t *bs) {
    free(bs->name);
    free(bs->comment);
    free(bs->seq);
    free(bs->qual);
    free(bs->sam);
}

#ifdef __cplusplus
}
#endif

#endif
