#ifndef _DATABASE_H__
#define _DATABASE_H__

#include "lib/encoder.h"
#include "lib/util.h"
#include <cinttypes>
#include <forward_list>
#include <unordered_set>

#define __fr(item, fp) std::fread(&(item), 1, sizeof(item), fp)
#define __fw(item, fp) std::fwrite(&(item), 1, sizeof(item), fp)

namespace emp {


template <typename T>
struct Database {

    unsigned k_, w_;
    T       *db_;
    int      owns_hash_;
    spvec_t  s_;
    Spacer  *sp_;

    Spacer *make_sp() {
        Spacer *ret(new Spacer(k_, w_, s_));
        for(auto &i: ret->s_) --i;
        return ret;
    }

    Database(const char *fn): owns_hash_(1), sp_(nullptr) {
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
        sp_ = make_sp();
        LOG_DEBUG("Read database!\n");
        std::fclose(fp);
    }
    Database(unsigned k, unsigned w, const spvec_t &s, unsigned owns=1, T *db=nullptr):
        k_(k), w_(w), db_(db), owns_hash_(owns), s_(s), sp_(make_sp())
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
        s_(other.s_),
        sp_(make_sp())
    {
    }

    ~Database() {
        if(owns_hash_) khash_destroy(db_);
        if(sp_)        delete sp_;
    }
    void write(const char *fn) {
        std::FILE *ofp(std::fopen(fn, "wb"));
        if(!ofp) LOG_EXIT("Could not open %s for reading.\n", fn);
        __fw(k_, ofp);
        __fw(w_, ofp);
        std::fwrite(s_.data(), s_.size(), sizeof(uint8_t), ofp);
        khash_write_impl<T>(db_, ofp);
        std::fclose(ofp);
#if DATABASE_TEST
        Database<T> test(fn);
        assert(kh_size(test.db_) == kh_size(db_));
        for(khiter_t ki(0); ki != kh_end(test.db_); ++ki) {
            if(kh_key(db_, ki) != kh_key(test.db_, ki))
                std::fprintf(stderr, "key mismatch at %" PRIu64 ". key 1: %" PRIu64 ", key 2: %" PRIu64 ".\n", ki, (std::uint64_t)kh_key(db_, ki), (std::uint64_t)kh_key(test.db_, ki));
            if(kh_val(db_, ki) != kh_val(test.db_, ki))
                std::fprintf(stderr, "val mismatch at %" PRIu64 ". val 1: %" PRIu64 ", val 2: %" PRIu64 ".\n", ki, (std::uint64_t)kh_val(db_, ki), (std::uint64_t)kh_val(test.db_, ki));
        }
#endif
    }

    template<typename Q=T>
    typename std::enable_if<std::is_same<khash_t(c), Q>::value, std::uint32_t>::type
    get_lca(std::uint64_t kmer) {
        khiter_t ki;
        return ((ki = kh_get(c, db_, kmer)) == kh_end(db_)) ? kh_val(db_, ki)
                                                            : -1u;
    }
};

template<typename T>
struct val_data_t {
    const std::unordered_map<tax_t, std::forward_list<std::string>> &tx_;
    const std::vector<tax_t>                                          v_;
    const Database<T>                                               &db_;
};

template<typename T>
void val_helper(void *data_, long index, int tid) {
    val_data_t<T> &data(*static_cast<val_data_t<T> *>(data_));
    const tax_t tax(data.v_[index]);
    if(!has_key(tax, data.tx_)) return;
    std::unordered_set<std::uint64_t> in;
    std::unordered_set<std::uint64_t> failing_kmers;
    Encoder<lex_score> enc(nullptr, 0, *data.db_.sp_, nullptr);
    for(const auto &pair: data.tx_) {
        for(const auto &path: pair.second) {
            LOG_DEBUG("Validator opening path at %s\n", path.data());
            gzFile fp(gzopen(path.data(), "rb"));
            kseq_t *ks(kseq_init(fp));
            std::uint64_t kmer;
            while(kseq_read(ks) >= 0) {
                enc.assign(ks);
                while(enc.has_next_kmer()) if((kmer = enc.next_kmer()) != BF) in.insert(kmer);
            }
        }
    }
    LOG_DEBUG("Validator loaded all kmers from end genome. Now scanning full database for items assigned to tax. (tax: %u. tid; %i)\n", tax, tid);
    for(std::uint64_t i(0); i < kh_size(data.db_.db_); ++i) if(kh_val(data.db_.db_, i) == tax) {
        if((i & 0xFFFFF) == 0) LOG_DEBUG("Processed %" PRIu64 " of %" PRIu64"\n.", i, kh_size(data.db_.db_));
        if(in.find(kh_key(data.db_.db_, i)) == in.end()) {
            failing_kmers.insert(kh_key(data.db_.db_, i));
            LOG_DEBUG("Missing kmer %s (%" PRIu64") from genome that was assigned as lca (%u)\n", data.db_.sp_->to_string(kh_key(data.db_.db_, i)).data(), kh_key(data.db_.db_, i), tax);
        }
    }
    if(failing_kmers.size()) {
        LOG_DEBUG("%zu kmers failed.\n", failing_kmers.size());
        for(const auto kmer: failing_kmers)
            LOG_DEBUG("Failed kmer %s/%" PRIu64"\n",
                      data.db_.sp_->to_string(kh_key(data.db_.db_, kmer)).data(),
                      kh_key(data.db_.db_, kmer));
    } else LOG_DEBUG("Tax %u validated.\n", tax);
}


template<typename T>
void validate_db(const Database<T> &db, std::unordered_set<tax_t> &used_lcas, std::unordered_map<tax_t, std::forward_list<std::string>> &tx2g, int num_threads=-1) {
    if(num_threads < 0) num_threads = std::thread::hardware_concurrency();
    std::unordered_set<std::uint64_t> failing_kmers;
    std::vector<tax_t> lcas(used_lcas.begin(), used_lcas.end());
    val_data_t<T> data{tx2g, lcas, db};
    kt_for(num_threads, &val_helper<T>, &data, lcas.size());
}


} /* emp namespace */

#undef __fr
#undef __fw


#endif /* ifndef _DATABASE_H__ */

