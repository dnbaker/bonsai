#ifndef __DB_IO_H__
#define __DB_IO_H__


/* CVS Version Info
   $Id: $
*/

/*! @file db_io.h
    @brief
    @author Aaron Carass
    @created Sat 19 Nov 2016 12:57:32 PM EST
    @edited Tue 22 Nov 2016 02:56:07 PM EST


*/


#include "lib/spacer.h"
#include <cinttypes>

#define __fr(item, fp) fread(&(item), 1, sizeof(item), fp)
#define __fw(item, fp) fwrite(&(item), 1, sizeof(item), fp)

namespace emp {


template <typename T>
struct Database {

    unsigned k_, w_;
    T *db_;
    int owns_hash_;
    spvec_t s_;

    Database(const char *fn): owns_hash_(1) {
        FILE *fp(fopen(fn, "rb"));
        if (fp) {
            __fr(k_, fp);
            __fr(w_, fp);
            s_ = spvec_t(k_ - 1);
            LOG_DEBUG("reading %zu bytes from file for vector, with %zu reserved\n", s_.size(), s_.capacity());
            fread(s_.data(), s_.size(), sizeof(uint8_t), fp);
    #if !NDEBUG
            LOG_DEBUG("vec size: %u\n", s_.size());
            for(auto i: s_) fprintf(stderr, "Value in vector is %u\n", (unsigned)i);
    #endif
            db_ = (T *)calloc(1, sizeof(T));
            __fr(*db_,   fp);
    #if !NDEBUG
            print_khash(db_);
    #endif

            LOG_DEBUG("Doing flags. N buckets: %zu. Elements to read: %zu\n", db_->n_buckets, __ac_fsize(db_->n_buckets));
            db_->flags = (uint32_t *)malloc(sizeof(*db_->flags) * __ac_fsize(db_->n_buckets));
            if(!db_->flags) LOG_EXIT("Could not allocate memory for flags.\n");
            fread(db_->flags, __ac_fsize(db_->n_buckets), sizeof(*(db_->flags)), fp);

            LOG_DEBUG("Doing keys\n");
            typedef typename std::remove_pointer<decltype(db_->keys)>::type keytype_t;
            db_->keys = (keytype_t *)malloc(sizeof(*db_->keys) * db_->n_buckets);
            fread(db_->keys, db_->n_buckets, sizeof(*db_->keys), fp);

            LOG_DEBUG("Doing vals\n");
            typedef typename std::remove_pointer<decltype(db_->vals)>::type valtype_t;
            db_->vals = (valtype_t *)malloc(sizeof(*db_->keys) * db_->n_buckets);
            fread(db_->vals, db_->n_buckets, sizeof(*db_->vals), fp);
        } else LOG_EXIT("Could not open %s for reading.\n", fn);
        LOG_DEBUG("Read database!\n");
        fclose(fp);
    }
    Database(unsigned k, unsigned w, spvec_t &s, unsigned owns=1, T *db=nullptr):
        k_(k), w_(w), db_(db), owns_hash_(owns), s_(s)
    {
    }
    Database(Spacer sp, unsigned owns=1, T *db=nullptr):
        Database(sp.k_, sp.w_, sp.s_, owns, db)
    {
    }

    template<typename O>
    Database(Database<O> &other, unsigned owns=0):
        k_(other.k_),
        w_(other.w_),
        db_(nullptr),
        owns_hash_(owns),
        s_(other.s_)
    {
    }

    ~Database() {
        if(owns_hash_) khash_destroy(db_);
    }
    void write(const char *fn) {
        FILE *ofp(fopen(fn, "wb"));
        if(!ofp) LOG_EXIT("Could not open %s for reading.\n", fn);
        for(khiter_t ki(0); ki != kh_end(db_); ++ki)
            if(!kh_exist(db_, ki))
                kh_key(db_, ki) = kh_val(db_, ki) = 0;
        __fw(k_, ofp);
        __fw(w_, ofp);
        fwrite(s_.data(), s_.size(), sizeof(uint8_t), ofp);
        fwrite(db_, 1, sizeof(*db_), ofp);
        fwrite(db_->flags, __ac_fsize(db_->n_buckets), sizeof(*db_->flags), ofp);
        fwrite(db_->keys, db_->n_buckets, sizeof(*db_->keys), ofp);
        fwrite(db_->vals, db_->n_buckets, sizeof(*db_->vals), ofp);
        fclose(ofp);
    }
};


} /* emp namespace */

#undef __fr
#undef __fw


#endif /* ifndef __DB_IO_H__ */

