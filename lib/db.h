#ifndef __DB_H__
#define __DB_H__

#include "feature_min.h"
#include "encoder.h"
#include "jellyfish/hash_counter.hpp"
//#include "growt/data-structures/definitions.h"

typedef jellyfish::cooperative::hash_counter<uint64_t> jfhash_t;

namespace kpg {


template<uint64_t (*score)(uint64_t, void *)>
void build_db(khash_t(64) *hash, std::vector<std::string> &paths, const std::string &outpath) {
    // Build database, write it to disk at *file*.
}

jfhash_t load_table(const std::string &path) {
    // Somehow build hash.
    // Make file header, then dumper.
    // Then load the array of keys and values into memory
    // Then build the hash table
    return jfhash_t(1 << 20, 62, 64, 16, 126);
}

template <uint64_t (*score)(uint64_t, void *)>
class Classifier {
    // Further, this needs a vector of Encoder objects
    static const size_t reprobe_limit = 126;

    const size_t default_nelem_;
    const unsigned nt_;
    const unsigned k_;
    std::vector<jfhash_t> hashes_;
    std::vector<Encoder<score>> encoders_;
    public:
    Classifier(size_t nelem, uint8_t k, int num_threads):
        default_nelem_(nelem),
        nt_(num_threads > 0 ? num_threads: 16),
        k_(k)
    {}
    ~Classifier() {
    }
    void add_encoder(Spacer &sp, void *data) {
        encoders_.emplace_back(nullptr, 0, sp, data);
    }
    void add_hash(size_t nelem=0) {
        hashes_.emplace_back(nelem ? nelem: default_nelem_, (k_ << 1), 64, nt_, reprobe_limit);
    }
};
#if 0
template<uint64_t (*score)(uint64_t, void *)>
int mindb_helper(const char *path, const Spacer &sp, void *data, chm_t &ret) {
    uint64_t kmer;
    khash_t(64) *td_map((khash_t(64) *)data);
    Encoder<score> enc(0, 0, sp, data);
    gzFile fp(gzopen(path, "rb"));
    kseq_t *ks(kseq_init(fp));
    khint_t ki;
    if(likely(kseq_read(ks) >= 0)) {
        enc.assign(ks);
        while(enc.has_next_kmer()) {
            if((kmer = enc.next_minimizer()) != BF) {
                if(unlikely((ki = kh_get(64, td_map, kmer)) == kh_end(td_map)))
                    LOG_EXIT("kmer not in taxdepth hash. This should never happen.\n");
                ret[kmer] = (uint32_t)kh_val(td_map, ki);
            }
        }
    }
    gzclose(fp);
    kseq_destroy(ks);
    return 1;
}

template<uint64_t (*score)(uint64_t, void *)>
void build_minimized_database(khash_t(64) *td_map, const Spacer &sp, std::vector<std::string> &paths, chm_t &ret, int num_threads) {
    if(num_threads < 0) num_threads = 16;
    ret.clear();
    std::vector<std::future<int>> futures;
    size_t i(0);
    while(i < (unsigned)num_threads && i < paths.size()) {
        futures.emplace_back(std::async(
          std::launch::async, mindb_helper<score>, paths[i++].data(), sp, (void *)td_map, std::ref(ret)));
    }
    while(futures.size()) {
        for(auto f(futures.begin()); f != futures.end(); ++f) {
            if(is_ready(*f)) {
                futures.erase(f);
                if(i < paths.size())
                    futures.emplace_back(std::async(
                        std::launch::async, mindb_helper<score>, paths[i++].data(), sp, (void *)td_map, std::ref(ret)));
                break;
            }
        }
    }
}
#endif

} // namespace kpg

#endif // #ifndef db
