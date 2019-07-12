#ifndef KMER_COUNTER_H__
#define KMER_COUNTER_H__
#include <cstdint>
#include <vector>
#include <cstdlib>
#include <type_traits>
#include "khash64.h"
#include "bonsai/include/encoder.h"

namespace kmerc {

KHASH_MAP_INIT_INT64(i16, uint16_t)

template<typename C, typename IT=uint64_t, typename ArgType>
std::vector<khash_t(i16)> build_kmer_counts(const C &kmer_sizes, ArgType fp, bool canon=false, size_t presize=0) {
    static_assert(std::is_same<ArgType, gzFile>::value  || std::is_same<ArgType, char *>::value || std::is_same<ArgType, const char *>::value, "Must be gzFile, char *, or const char *");
    bns::RollingHasherSet<IT> rhs(kmer_sizes, canon);
    using T = khash_t(i16);
    std::vector<T> kmer_maps(kmer_sizes.size());
    std::memset(&kmer_maps[0], 0, sizeof(kmer_maps[0]) * kmer_sizes.size());
    if(presize)
        for(auto &x: kmer_maps)
            kh_resize(i16, &x, presize);
    rhs.for_each_hash([&kmer_maps](IT hashvalue, size_t idx){
        auto map_ptr = &kmer_maps[idx];
        if((idx = kh_get(i16, map_ptr, hashvalue)) != kh_end(map_ptr)) {
            auto &val =  map_ptr->vals[idx];
            val += (val != std::numeric_limits<std::decay_t<decltype(*map_ptr->vals)>>::max());
        } else {
            int khr;
            int rc;
            idx = kh_put(i16, map_ptr, hashvalue, &khr);
            if(khr < 0) {
                std::fprintf(stderr, "Error: rc of %d, khr is %d\n", rc, khr);
                throw std::runtime_error("Failed to insert.");
            }
            LOG_ASSERT(idx < map_ptr->n_buckets);
            map_ptr->vals[idx] = 1u;
        }
    }, fp);
    return kmer_maps;
}

enum DumpFlags: int {
    WRITE_SHS = 1,
    WRITE_KVMAP = 2
};

template<typename C, typename IT=uint64_t, typename ArgType>
void dump_maps(const char *prefix, const C &kmer_sizes, ArgType fp, bool canon=false, size_t presize=0, int flag=WRITE_SHS | WRITE_KVMAP) {
    auto maps = build_kmer_counts(kmer_sizes, fp, canon, presize);
    std::vector<IT> buf;
    std::vector<uint16_t> buf16;
    for(size_t kszidx = 0; kszidx < kmer_sizes.size(); ++kszidx) {
        auto k = kmer_sizes[kszidx];
        auto &map = maps[kszidx];
        buf16.resize(kh_size(&map));
        buf.resize(kh_size(&map));
        size_t used = 0;
        for(size_t i = 0; i < map.n_buckets; ++i) {
            if(kh_exist(&map, i))
                buf[used] = map.keys[i], buf16[used] = map.vals[i], ++used;
        }
        const uint64_t count = used;
        gzFile fp;
        if(flag & WRITE_KVMAP) {
            fp = gzopen((std::string(prefix) + "." + std::to_string(k) + ".bin").data(), "wb");
            if(!fp) throw std::runtime_error("Could not open file.");
            gzwrite(fp, &count, sizeof(count));
            gzwrite(fp, buf.data(), buf.size() * sizeof(buf[0]));
            gzwrite(fp, buf16.data(), buf16.size() * sizeof(buf16[0]));
            gzclose(fp);
        }
        if(flag & WRITE_SHS) {
            std::sort(buf.data(), buf.data() + buf.size());
            fp = gzopen((std::string(prefix) + "." + std::to_string(k) + ".shs").data(), "wb");
            gzwrite(fp, &count, sizeof(count));
            gzwrite(fp, buf.data(), buf.size() * sizeof(buf[0]));
            gzclose(fp);
        }
        std::free(map.keys);
        std::free(map.vals);
        std::free(map.flags);
    }
}


} // kmerc

#endif
