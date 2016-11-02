#ifndef _KSEQ_DECLARE_
#define _KSEQ_DECLARE_
#include <zlib.h>
#include "htslib/kseq.h"
#include <cstdlib>
#include <cstdint>

KSEQ_INIT(gzFile, gzread)

extern "C" {

typedef struct {
    int l_seq, id;
    char *name, *comment, *seq, *qual, *sam;
} bseq1_t;

// From bwa -- batch parsing.
static inline void kseq2bseq1(const kseq_t *ks, bseq1_t *s)
{ // TODO: it would be better to allocate one chunk of memory, but probably it does not matter in practice
	s->name = strdup(ks->name.s);
	s->comment = ks->comment.l? strdup(ks->comment.s) : 0;
	s->seq = strdup(ks->seq.s);
	s->qual = ks->qual.l? strdup(ks->qual.s) : 0;
	s->l_seq = strlen(s->seq);
}

bseq1_t *bseq_read(int chunk_size, int *n_, void *ks1_, void *ks2_);
}

#endif
