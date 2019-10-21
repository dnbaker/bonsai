#ifndef _DATABASE_H__
#define _DATABASE_H__

#include "encoder.h"
#include "util.h"
#include <cinttypes>
#include <forward_list>
#include <unordered_set>

#define __fr(item, fp) if(std::fread(&(item), 1, sizeof(item), fp) != sizeof(item)) throw std::runtime_error("Error: Could not read " #item);
#define __fw(item, fp) if(std::fwrite(&(item), 1, sizeof(item), fp) != sizeof(item)) throw std::runtime_error("Error: Could not write " # item);

namespace bns {


template <typename T>
struct Database {

    unsigned k_, w_;
    T       *db_;
    int      owns_hash_;
    spvec_t  s_;
    Spacer  *sp_;

    Spacer *make_sp() {
        //std::fprintf(stderr, "Making sp with spacer = %s\n", str(s_).data());
        Spacer *ret(new Spacer(k_, (uint16_t)w_, s_));
        for(auto &i: ret->s_) --i;
        //std::fprintf(stderr, "Current sp string: %s\n", str(ret->s_).data());
        return ret;
    }

    Database(const char *fn): owns_hash_(1), sp_(nullptr) {
        int filetype(0);
        {
            std::string fns = fn;
            std::string gzsuf   = ".gz";
            std::string zstdsuf = ".zst";
            if(std::equal(std::crbegin(gzsuf), std::crend(gzsuf), std::crbegin(fns))) filetype = 1;
            else if(std::equal(std::crbegin(gzsuf), std::crend(gzsuf), std::crbegin(fns))) filetype = 2;
        }
        std::FILE *fp = filetype ? popen((std::string(filetype == 1 ? "gzip -dc " : "zstd -qdc ") + fn).data(), "rb"): std::fopen(fn, "rb");
        if (fp) {
            __fr(k_, fp);
            __fr(w_, fp);
            s_ = spvec_t(k_ - 1);
            LOG_DEBUG("reading %zu bytes from file for vector, with %zu reserved\n", s_.size(), s_.capacity());
            if(std::fread(s_.data(), s_.size(), sizeof(uint8_t), fp) != s_.size() * sizeof(uint8_t))
                throw std::runtime_error("Error: Could not read spacing from file");
            db_ = khash_load_impl<T>(fp);
        } else LOG_EXIT("Could not open %s for reading.\n", fn);
        sp_ = make_sp();
        assert(sp_);
        LOG_DEBUG("Read database!\n");
        std::fclose(fp);
    }
    Database(unsigned k, unsigned w, const spvec_t &s, unsigned owns=1, T *db=nullptr):
        k_(k), w_(w), db_(db), owns_hash_(owns), s_(s), sp_(make_sp())
    {
    }
    Database(Spacer sp, unsigned owns=1, T *db=nullptr):
        Database(sp.k_, sp.w_, sub1(sp.s_), owns, db)
    {
    }

    template<typename O>
    Database(Database<O> &other, unsigned owns=0):
        k_(other.k_),
        w_(other.w_),
        db_(nullptr),
        owns_hash_(owns),
        s_(other.s_),
        sp_(make_sp())
    {
    }

    ~Database() {
        if(owns_hash_) khash_destroy(db_);
        if(sp_)        delete sp_;
    }
    void write(const char *fn, bool write_gz=false) const {
        // TODO: add compression/work with zlib.
        if(write_gz) {
            gzFile ofp = gzopen(fn, "wb");
            if(!ofp) LOG_EXIT("Could not open %s for writing.\n", fn);
#define gzw(_x, ofp) if(gzwrite(ofp, static_cast<const void *>(&_x), sizeof(_x)) != sizeof(_x)) throw std::runtime_error("Error writing to file")
            gzw(k_, ofp);
            gzw(w_, ofp);
            gzwrite(ofp, static_cast<const void *>(s_.data()), s_.size() * sizeof(s_[0]));
            khash_write_impl<T>(db_, ofp);
            gzclose(ofp);
            return;
#undef gzw
        } // else
        std::FILE *ofp(std::fopen(fn, "wb"));
        if(!ofp) LOG_EXIT("Could not open %s for writing.\n", fn);
        __fw(k_, ofp);
        __fw(w_, ofp);
        if(std::fwrite(s_.data(), s_.size(), sizeof(uint8_t), ofp) != s_.size()) throw std::runtime_error("Error writing database");
        khash_write_impl<T>(db_, ofp);
        std::fclose(ofp);
    }

    template<typename Q=T>
    typename std::enable_if<std::is_same<khash_t(c), Q>::value, u32>::type
    get_lca(u64 kmer) {
        khiter_t ki;
        return ((ki = kh_get(c, db_, kmer)) == kh_end(db_)) ? kh_val(db_, ki)
                                                            : -1u;
    }
};

} /* bns namespace */

#undef __fr
#undef __fw


#endif /* ifndef _DATABASE_H__ */

