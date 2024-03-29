#ifndef KMER_COUNTER_H__
#define KMER_COUNTER_H__
#include <cstdint>
#include <vector>
#include <cstdlib>
#include <type_traits>
#include "khash64.h"
#include "./encoder.h"
#include "sketch/common.h"

namespace kmerc {

KHASH_MAP_INIT_INT64(i16, uint16_t)
KHASH_SET_INIT_INT64(i16s)

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
            idx = kh_put(i16, map_ptr, hashvalue, &khr);
            if(khr < 0) {
                std::fprintf(stderr, "Error: khr is %d\n", khr);
                throw std::runtime_error("Failed to insert.");
            }
            LOG_ASSERT(idx < map_ptr->n_buckets);
            map_ptr->vals[idx] = 1u;
        }
    }, fp);
    return kmer_maps;
}

template<typename C, typename IT=uint64_t, typename ArgType,typename Allocator=sketch::Allocator<IT>>
std::vector<std::vector<IT, Allocator>> build_kmer_sets(const C &kmer_sizes, ArgType fp, bool canon=false, size_t presize=0) {
    static_assert(std::is_same<ArgType, gzFile>::value  || std::is_same<ArgType, char *>::value || std::is_same<ArgType, const char *>::value, "Must be gzFile, char *, or const char *");
    bns::RollingHasherSet<IT> rhs(kmer_sizes, canon);
    using T = std::vector<IT, Allocator>;
    std::vector<T> kmer_sets(kmer_sizes.size());
    if(presize)
        for(auto &x: kmer_sets)
            x.reserve(presize);
    rhs.for_each_hash([&kmer_sets](IT hashvalue, size_t idx){kmer_sets[idx].push_back(hashvalue);}, fp);
    OMP_PRAGMA("omp parallel for")
    for(size_t i = 0; i < kmer_sizes.size(); ++i) {
        auto &v = kmer_sets[i];
        std::sort(v.begin(), v.end()); // TODO: provide support for other sorting methods
        v.erase(std::unique(v.begin(), v.end()), v.end());
    }
    return kmer_sets;
}

enum DumpFlags: int {
    WRITE_SHS = 1,
    WRITE_KVMAP = 2
};

struct WriteFail {};
struct OpenFail {};

template<typename C, typename IT, typename ArgType>
void dump_shs(const char *prefix, const C &kmer_sizes, ArgType cfp, bool canon, size_t presize=0) {
    auto shsets = build_kmer_sets(kmer_sizes, cfp, canon, presize);
    std::atomic<int> ret;
    ret.store(0);
    //#pragma omp parallel for
    for(size_t i = 0; i < kmer_sizes.size(); ++i) {
        auto &vec = shsets[i];
        auto k = kmer_sizes[i];
        std::string fn = std::string(prefix) + "." + std::to_string(k) + ".shs";
        gzFile fp = gzopen(fn.data(), "wb");
        if(!fp) throw std::runtime_error(std::string("Could not open file at ") + fn + " for writing");
        uint64_t nelem = vec.size();
        if(gzwrite(fp, &nelem, sizeof(nelem)) != sizeof(nelem)) ret.store(1 << (i % 64));
        ssize_t nb = sizeof(vec[0]) * vec.size();
        if(gzwrite(fp, vec.data(), nb) != nb) ret.store(1 << (i % 64));
    }
    if(ret) {
        throw WriteFail{};//std::runtime_error("Failed to write");
    }
}

template<typename C, typename IT=uint64_t, typename ArgType>
void dump_maps(const char *prefix, const C &kmer_sizes, ArgType fp, bool canon=false, size_t presize=0, int flag=WRITE_SHS | WRITE_KVMAP) {
    if(flag == WRITE_SHS) {
        dump_shs<C, IT, ArgType>(prefix, kmer_sizes, fp, canon);
    }
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
