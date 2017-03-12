#ifndef _DATABASE_H__
#define _DATABASE_H__

#include "lib/spacer.h"
#include <cinttypes>

#define __fr(item, fp) std::fread(&(item), 1, sizeof(item), fp)
#define __fw(item, fp) std::fwrite(&(item), 1, sizeof(item), fp)

namespace emp {


template <typename T>
struct Database {

    unsigned k_, w_;
    T *db_;
    int owns_hash_;
    spvec_t s_;

    Database(const char *fn): owns_hash_(1) {
        std::FILE *fp(std::fopen(fn, "rb"));
        if (fp) {
            __fr(k_, fp);
            __fr(w_, fp);
            s_ = spvec_t(k_ - 1);
            LOG_DEBUG("reading %zu bytes from file for vector, with %zu reserved\n", s_.size(), s_.capacity());
            std::fread(s_.data(), s_.size(), sizeof(uint8_t), fp);
    #if !NDEBUG
            LOG_DEBUG("vec size: %u\n", s_.size());
            for(auto i: s_) fprintf(stderr, "Value in vector is %u\n", (unsigned)i);
    #endif
            db_ = khash_load_impl<T>(fp);
        } else LOG_EXIT("Could not open %s for reading.\n", fn);
        LOG_DEBUG("Read database!\n");
        std::fclose(fp);
    }
    Database(unsigned k, unsigned w, const spvec_t &s, unsigned owns=1, T *db=nullptr):
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
        std::FILE *ofp(std::fopen(fn, "wb"));
        if(!ofp) LOG_EXIT("Could not open %s for reading.\n", fn);
        for(khiter_t ki(0); ki != kh_end(db_); ++ki)
            if(!kh_exist(db_, ki))
                kh_key(db_, ki) = kh_val(db_, ki) = 0;
        __fw(k_, ofp);
        __fw(w_, ofp);
        std::fwrite(s_.data(), s_.size(), sizeof(uint8_t), ofp);
        std::fwrite(db_, 1, sizeof(*db_), ofp);
        std::fwrite(db_->flags, __ac_fsize(db_->n_buckets), sizeof(*db_->flags), ofp);
        std::fwrite(db_->keys, db_->n_buckets, sizeof(*db_->keys), ofp);
        std::fwrite(db_->vals, db_->n_buckets, sizeof(*db_->vals), ofp);
        std::fclose(ofp);
    }

    template<typename Q=T>
    typename std::enable_if<std::is_same<khash_t(c), Q>::value, std::uint32_t>::type
    get_lca(std::uint64_t kmer) {
        khiter_t ki;
        return ((ki = kh_get(c, db_, kmer)) == kh_end(db_)) ? kh_val(db_, ki)
                                                            : -1u;
    }
};


} /* emp namespace */

#undef __fr
#undef __fw


#endif /* ifndef _DATABASE_H__ */

