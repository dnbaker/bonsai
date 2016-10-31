#ifndef __DB_H__
#define __DB_H__

#include "feature_min.h"
#include "cuckoohash_map.hh"
#include "city_hasher.hh"

namespace kpg {
typedef cuckoohash_map<uint64_t, uint32_t, CityHasher<uint64_t>> chm_t;

int mindb_helper(std::string &path, const Spacer &sp, void *data, chm_t &ret);

template<uint64_t (*score)(uint64_t, void *)>
int mindb_helper(std::string &path, const Spacer &sp, void *data, chm_t &ret) {
    uint64_t kmer;
    khash_t(64) *td_map((khash_t(64) *)data);
    Encoder<score> enc(0, 0, sp, data);
    gzFile fp(gzopen(path.data(), "rb"));
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
    while(i < num_threads && i < paths.size()) {
        futures.emplace_back(std::async(
          std::launch::async, mindb_helper<score>, paths[i++], sp, (void *)td_map, ret));
    }
    while(futures.size()) {
        for(auto f(futures.begin()); f != futures.end(); ++f) {
            if(is_ready(*f)) {
                futures.erase(f);
                if(i < paths.size())
                    futures.emplace_back(std::async(
                        std::launch::async, mindb_helper<score>, paths[i++], sp, (void *)td_map, ret));
                break;
            }
        }
    }
}

} // namespace kpg

#endif // #ifndef db
