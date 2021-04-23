#include <getopt.h>
#include <unistd.h>
#include <sys/stat.h>

#include <cstdio>
#include <cstdint>
#include <vector>
#include <chrono>
#include <string>
#include <iostream>
#include <fstream>
#include <numeric>
#include <memory>

#ifdef _OPENMP
#include <omp.h>
#endif
#include "flat_hash_map/flat_hash_map.hpp"

using std::uint64_t;
using std::uint32_t;

int usage() {
    return 1;
}

// Step 1: load k-mer files
// Step 2: invert matrix
ska::flat_hash_map<uint64_t, std::vector<uint32_t>>
read_file(gzFile fp) {
    auto timestart = std::chrono::high_resolution_clock::now();
    uint64_t arr[2];
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
        //if(i % 256 == 0) std::fprintf(stderr, "%zu/%zu, read %zu\n", i + 1, size_t(arr[1]), total_ids_read);
        const auto nids = data[i];
        buffer.resize(nids);
        gzread(fp, buffer.data(), sizeof(uint32_t) * nids);
        total_ids_read += nids;
        map.emplace(keys[i], buffer);
    }
    std::fprintf(stderr, "Time to deserialize: %gms\n", std::chrono::duration<double, std::milli>(std::chrono::high_resolution_clock::now() - timestart).count());
    return map;
}

int main(int argc, char **argv) {
    bool basename = false;
    int ret = 0;
    std::vector<std::string> names;
    std::FILE *ofp = stdout;
    for(int c;(c = getopt(argc, argv, "o:F:bh")) >= 0;) switch(c) {
        case 'b': basename = true; break;
        case 'h': case '?': return usage();
        case 'o': if(!(ofp = std::fopen(optarg, "w"))) {std::fprintf(stderr, "Could not open file at %s\n", optarg); std::abort();}
                  break;
        case 'F': {
            std::ifstream ifs(optarg);
            for(std::string s; std::getline(ifs, s);)
                names.push_back(s);
        }
    }
    const size_t nfiles = names.size();
    if(argc == optind) throw 1;
    auto map = load_map(argv[optind]);
    ska::flat_hash_map<uint64_t, uint32_t> counter;
    counter.reserve(map.size());
    for(const auto &pair: map) counter.emplace(pair.first, 0u);
    std::fprintf(stderr, "map size %zu and total number of ids %zu\n", map.size(), std::accumulate(map.begin(), map.end(), size_t(0), [](size_t x, auto &y) {return x + y.second.size();}));
    if(names.empty()) return usage();
    for(const auto &fn: names) {
        
    }
    return 0;
}
