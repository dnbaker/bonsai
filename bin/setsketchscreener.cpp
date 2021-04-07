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
ska::flat_hash_map<uint64_t, std::vector<uint32_t>> load_map(std::string path) {
    std::FILE *fp = std::fopen(path.data(), "rb");
    uint64_t arr[2];
    std::fread(arr, sizeof(arr), 1, fp);
    uint64_t totalkmers = arr[0], totalids = arr[1];
    std::unique_ptr<uint32_t[]> ids, idn;
    std::unique_ptr<uint64_t[]> kmers(new uint64_t[totalkmers]);
    idn.reset(new uint32_t[totalkmers]);

    if(std::fread(idn.get(), sizeof(uint32_t), totalkmers, fp) != totalkmers) throw 1;
    if(std::fread(kmers.get(), sizeof(uint64_t), totalkmers, fp) != totalkmers) throw 1;

    ids.reset(new uint32_t[totalids]);
    
    if(std::fread(ids.get(), sizeof(uint32_t), totalkmers, fp) != totalkmers) throw 1;
    ska::flat_hash_map<uint64_t, std::vector<uint32_t>> ret;
    ret.reserve(totalkmers);
    size_t idstart = 0;
    for(size_t i = 0; i < totalkmers; ++i) {
        ret.emplace(kmers[i], std::vector<uint32_t>(&ids[idstart], &ids[idstart + idn[i]]));
        idstart += idn[i];
    }
    return ret;
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
    std::fprintf(stderr, "map size %zu and total number of ids %zu\n", map.size(), std::accumulate(map.begin(), map.end(), size_t(0), [](size_t x, auto &y) {return x + y.second.size();}));
    if(names.empty()) return usage();
    return 0;
}
