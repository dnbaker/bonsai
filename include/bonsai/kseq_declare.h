#ifndef _KSEQ_DECLARE__
#define _KSEQ_DECLARE__
#if ZWRAP_USE_ZSTD
#  include "zstd_zlibwrapper.h"
#else
#  include <zlib.h>
#endif
#include "klib/kseq.h"
#include <stdlib.h>
#include <stdint.h>
#include <stdio.h>
#include <assert.h>

#define KSEQ_INIT3(SCOPE, type_t, __read, SIZE)		\
	KSTREAM_INIT(type_t, __read, 16384)			\
	__KSEQ_TYPE(type_t)							\
	__KSEQ_BASIC(SCOPE, type_t)					\
	__KSEQ_READ(SCOPE)

#define KSEQ_INIT_SIZE(type_t, __read, SIZE) KSEQ_INIT3(static, type_t, __read, SIZE)

#ifndef KSTREAM_SIZE
#  define KSTREAM_SIZE (65536u)
#endif

KSEQ_INIT_SIZE(gzFile, gzread, KSTREAM_SIZE)

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
static inline void kseq2bseq1(const kseq_t *ks, bseq1_t *s) // one chunk
{
    s->name = (char *)malloc(ks->name.l + ks->comment.l + ks->seq.l + ks->qual.l + 4);
    memcpy(s->name, ks->name.s, ks->name.l + 1);
    char *start = s->name + ks->name.l + 2;
    if(ks->comment.l == 0) s->comment = NULL;
    else {
        s->comment = start;
        memcpy(s->comment, ks->comment.s, ks->comment.l + 1);
        start += ks->comment.l + 1;
    }
    s->seq = start;
    memcpy(s->seq, ks->seq.s, ks->seq.l + 1);
    start += ks->seq.l + 1;
    if(ks->qual.l == 0) s->qual = NULL;
    else s->qual = start, memcpy(s->qual, ks->qual.s, ks->qual.l + 1);
    s->l_seq   = ks->seq.l;
    s->sam     = NULL;
    s->l_sam   = s->id = 0;
}

static inline void rekseq2bseq1(const kseq_t *ks, bseq1_t *s) // one chunk
{
    if(s->name == NULL) {
        assert(!s->sam && !s->l_sam);
        kseq2bseq1(ks, s);
        return;
    }
    s->name = (char *)realloc(s->name, ks->name.l + ks->comment.l + ks->seq.l + ks->qual.l + 4);
    memcpy(s->name, ks->name.s, ks->name.l + 1);
    char *start = s->name + ks->name.l + 2;
    if(ks->comment.l == 0) s->comment = NULL;
    else {
        s->comment = start;
        memcpy(s->comment, ks->comment.s, ks->comment.l + 1);
        start += ks->comment.l + 1;
    }
    s->seq = start;
    memcpy(s->seq, ks->seq.s, ks->seq.l + 1);
    start += ks->seq.l + 1;
    if(ks->qual.l == 0) s->qual = NULL;
    else s->qual = start, memcpy(s->qual, ks->qual.s, ks->qual.l + 1);
    s->l_seq   = ks->seq.l;
    //s->sam     = NULL;
    //s->l_sam   = s->id = 0;
}

static inline void bseq_destroy(bseq1_t *bs) {
    free(bs->name);
    free(bs->sam);
}

static inline void trim_readno(kstring_t *s)
{
    if (s->l > 2 && s->s[s->l-2] == '/' && isdigit(s->s[s->l-1]))
        s->l -= 2, s->s[s->l] = 0;
}

static bseq1_t *bseq_read(int chunk_size, int *n_, void *ks1_, void *ks2_)
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
            seqs = (bseq1_t *)realloc(seqs, m * sizeof(bseq1_t));
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


static bseq1_t *bseq_realloc_read(int chunk_size, int *n_, void *ks1_, void *ks2_, bseq1_t *seqs) {
    if(!seqs) return bseq_read(chunk_size, n_, ks1_, ks2_);
    int n = 0, size = 0;
    kseq_t *ks = (kseq_t *)ks1_, *ks2 = (kseq_t *)ks2_;
    while (kseq_read(ks) >= 0) {
        if (ks2 && kseq_read(ks2) < 0) { // the 2nd file has fewer reads
            fprintf(stderr, "[W::%s] the 2nd file has fewer sequences.\n", __func__);
            break;
        }
        trim_readno(&ks->name);
        rekseq2bseq1(ks, seqs + n);
        seqs[n].id = n;
        size += seqs[n++].l_seq;
        if (ks2) {
            trim_readno(&ks2->name);
            rekseq2bseq1(ks2, seqs + n);
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

static INLINE void kseq_assign(kseq_t *ks, gzFile fp) {
    if(!ks->f) {
        ks->f = (kstream_t*)calloc(1, sizeof(kstream_t));
        ks->f->buf = (unsigned char*)malloc(KSTREAM_SIZE);
    } else {
        ks->f->is_eof = ks->f->begin = ks->f->end = 0;
    }
    ks->f->f = fp;
}

static inline kseq_t kseq_init_stack() {
    kseq_t ret;
    memset(&ret, 0, sizeof(ret));
    return ret;
}

static INLINE void kseq_destroy_stack(kseq_t &ks) {
    free(ks.name.s); free(ks.comment.s); free(ks.seq.s); free(ks.qual.s);
    ks_destroy(ks.f);
    memset(&ks, 0, sizeof(ks));
}


#ifdef __cplusplus
}
#endif

#endif // #ifndef _KSEQ_DECLARE__
