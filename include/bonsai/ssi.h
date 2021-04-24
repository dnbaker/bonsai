#ifndef BONSAI_SETSKETCHINDEX_H__
#define BONSAI_SETSKETCHINDEX_H__
#include "sketch/ssi.h"
#include <cstdint>
#include <zlib.h>
#include <unistd.h>

namespace bns {

namespace lsh {
using std::int64_t;

using sketch::SetSketchIndex;
static inline int write_database(const ska::flat_hash_map<uint64_t, std::vector<uint32_t>> &map, std::FILE *ofp, int k) {
    int ret = 0;
    const int fd = ::fileno(ofp);
    if(::write(fd, &k, 4) != 4) throw std::runtime_error("Failed to write k");
    const uint64_t total_vals = map.size();
    const uint64_t total_ids = std::accumulate(map.begin(), map.end(), uint64_t(0), [](uint64_t x, auto &pair) {return x + pair.second.size();});
    std::vector<uint32_t> nids_per_kmer(total_vals);
    {
        auto mit = map.begin();
        for(auto vit = nids_per_kmer.begin(), vie = nids_per_kmer.end(); vit != vie; ++mit, ++vit)
            *vit = mit->second.size();
    }
    static_assert(sizeof(map.begin()->first) == sizeof(uint64_t), "sanity check");
    uint64_t arr [] {total_vals, total_ids};
    int64_t npknb;
    if(::write(fd, arr, sizeof(arr)) != sizeof(arr))
        goto write_error;
    npknb = sizeof(uint32_t) * nids_per_kmer.size();
    if(::write(fd, nids_per_kmer.data(), npknb) != npknb) goto write_error;
    // Write keys
    for(const auto &pair: map)
        if(::write(fd, &pair.first, sizeof(uint64_t)) != int64_t(sizeof(uint64_t)))
            goto write_error;
    for(const auto &pair: map) {
        if(::write(fd, pair.second.data(), sizeof(uint32_t) * pair.second.size()) != int64_t(sizeof(uint32_t) * pair.second.size()))
            goto write_error;
    }
    if(0) {
         write_error: {
             std::fprintf(stderr, "Failed to write to file");
             ret = 1;
         }
    }
    return ret;
}
static inline int write_database(const ska::flat_hash_map<uint64_t, std::vector<uint32_t>> &map, gzFile ofp, int k) {
    int ret = 0;
    const uint64_t total_vals = map.size();
    const uint64_t total_ids = std::accumulate(map.begin(), map.end(), uint64_t(0), [](uint64_t x, auto &pair) {return x + pair.second.size();});
    if(gzwrite(ofp, &k, 4) != 4) throw std::runtime_error("Failed to write k");
    std::vector<uint32_t> nids_per_kmer(total_vals);
    {
        auto mit = map.begin();
        for(auto vit = nids_per_kmer.begin(), vie = nids_per_kmer.end(); vit != vie; ++mit, ++vit)
            *vit = mit->second.size();
    }
    static_assert(sizeof(map.begin()->first) == sizeof(uint64_t), "sanity check");
    uint64_t arr [] {total_vals, total_ids};
    int64_t npknb;
    if(gzwrite(ofp, arr, sizeof(arr)) != sizeof(arr))
        goto write_error;
    npknb = sizeof(uint32_t) * nids_per_kmer.size();
    if(gzwrite(ofp, nids_per_kmer.data(), npknb) != npknb) goto write_error;
        goto write_error;
    // Write keys
    for(const auto &pair: map)
        if(gzwrite(ofp, &pair.first, sizeof(uint64_t)) != sizeof(uint64_t))
            goto write_error;
    for(const auto &pair: map) {
        if(gzwrite(ofp, pair.second.data(), sizeof(uint32_t) * pair.second.size()) != int64_t(sizeof(uint32_t) * pair.second.size()))
            goto write_error;
    }
    if(0) {
         write_error: {
             std::fprintf(stderr, "Failed to write to file");
             ret = 1;
         }
    }
    return ret;
}
// You can use this database for inverted screen sweeps
static inline
std::pair<ska::flat_hash_map<uint64_t, std::vector<uint32_t>>, int>
read_database(std::FILE *fp) {
    auto timestart = std::chrono::high_resolution_clock::now();
    int k;
    uint64_t arr[2];
    uint64_t fret;
    if((fret = std::fread(&k, 4, 1, fp)) != 1) throw std::runtime_error("Failed to read k from database");
    if((fret = std::fread(arr, 16, 1, fp)) != 1) {throw std::runtime_error(std::string("Expected to read 1 16-byte block and got ") + std::to_string(fret));}
    std::unique_ptr<uint32_t[]> data(new uint32_t[arr[0]]);
    if((fret = std::fread(data.get(), sizeof(uint32_t), arr[0], fp)) != arr[0]) {
        throw std::runtime_error(std::string("Expected to read ") + std::to_string(arr[0]) + "-block of u32s and got " + std::to_string(fret));
    }
    std::unique_ptr<uint64_t[]> keys(new uint64_t[arr[0]]);
    if(std::fread(keys.get(), sizeof(uint64_t), arr[0], fp) != arr[0]) throw 3;
    ska::flat_hash_map<uint64_t, std::vector<uint32_t>> map;
    std::vector<uint32_t> buffer;
    map.reserve(arr[0]);
    size_t total_ids_read = 0;
    for(size_t i = 0; i < arr[0]; ++i) {
        const auto nids = data[i];
        buffer.resize(nids);
        if((fret = std::fread(buffer.data(), sizeof(uint32_t), nids, fp)) != nids) {
            throw std::runtime_error(std::string("Expected to read ") + std::to_string(nids) + "-block of u32s and got " + std::to_string(fret));
        }
        total_ids_read += nids;
        map.emplace(keys[i], buffer);
    }
    std::fprintf(stderr, "Time to deserialize: %gms\n", std::chrono::duration<double, std::milli>(std::chrono::high_resolution_clock::now() - timestart).count());
    return std::make_pair(map, k);
}
static inline ska::flat_hash_map<uint64_t, std::vector<uint32_t>> &operator+=(ska::flat_hash_map<uint64_t, std::vector<uint32_t>> &lhs, const ska::flat_hash_map<uint64_t, std::vector<uint32_t>> &rhs)
{
    for(auto &&pair: rhs) {
        auto it = lhs.find(pair.first);
        if(it == lhs.end()) it = lhs.emplace(pair).first;
        else it->second.insert(it->second.end(), pair.second.begin(), pair.second.end());
    }
    return lhs;
}

static inline
std::pair<ska::flat_hash_map<uint64_t, std::vector<uint32_t>>, int>
read_database(gzFile fp) {
    auto timestart = std::chrono::high_resolution_clock::now();
    uint64_t arr[2];
    int k;
    gzread(fp, &k, 4);
    gzread(fp, arr, sizeof(arr));
    std::unique_ptr<uint32_t[]> data(new uint32_t[arr[0]]);
    gzread(fp, data.get(), sizeof(uint32_t) * arr[0]);
    std::unique_ptr<uint64_t[]> keys(new uint64_t[arr[0]]);
    gzread(fp, data.get(), sizeof(uint64_t) * arr[0]);
    ska::flat_hash_map<uint64_t, std::vector<uint32_t>> map;
    std::vector<uint32_t> buffer;
    map.reserve(arr[0]);
    size_t total_ids_read = 0;
    for(size_t i = 0; i < arr[0]; ++i) {
        const auto nids = data[i];
        buffer.resize(nids);
        gzread(fp, buffer.data(), sizeof(uint32_t) * nids);
        total_ids_read += nids;
        map.emplace(keys[i], buffer);
    }
    std::fprintf(stderr, "Time to deserialize: %gms\n", std::chrono::duration<double, std::milli>(std::chrono::high_resolution_clock::now() - timestart).count());
    return std::make_pair(map, k);
}
static inline
std::pair<ska::flat_hash_map<uint64_t, std::vector<uint32_t>>, int>
read_database(std::string s) {
    gzFile fp = gzopen(s.data(), "rb");
    if(!fp) throw std::runtime_error(std::string("Could not read from file ") + s);
    std::pair<ska::flat_hash_map<uint64_t, std::vector<uint32_t>>, int> ret = read_database(fp);
    gzclose(fp);
    return ret;
}

} // namespace lsh

} // namespace bns

#endif /* BONSAI_SETSKETCHINDEX_H__ */
