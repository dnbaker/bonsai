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

#ifdef _OPENMP
#include <omp.h>
#endif
#include "flat_hash_map/flat_hash_map.hpp"

using std::uint64_t;
using std::uint32_t;

int usage() {
    std::fprintf(stderr, "setsketchfolder <flags> sketch1.kmers sketch2.kmers...\n"
                         "-b\tTrim folder path from filenames\n"
                         "-F\tRead paths from <file> instesad of positional arguments\n");
    return 1;
}

// Step 1: load k-mer files
// Step 2: invert matrix

int main(int argc, char **argv) {
    bool basename = false;
    int ret = 0;
    std::vector<std::string> names;
    uint64_t kmer_size_in_bytes = 8;
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
    for(char **p = argv + optind;p < argv + argc; names.emplace_back(*p++));
    if(names.empty()) return usage();
    const size_t nfiles = names.size();
    std::vector<uint64_t *>data(nfiles);
    std::vector<uint64_t> sizes(nfiles);
    // Load files
    auto startload = std::chrono::high_resolution_clock::now();
    for(const auto &n: names) std::fprintf(stderr, "name: %s\n", n.data());
#ifdef _OPENMP
    #pragma omp parallel for
#endif
    for(size_t i = 0; i < nfiles; ++i) {
        std::FILE *fp;
        struct stat st;
        if((fp = std::fopen(names[i].data(), "rb")) == nullptr) {
            std::fprintf(stderr, "Failed to open file at %s\n", names[i].data());
            std::abort();
        }
        const int fd = ::fileno(fp);
        ::fstat(fd, &st);
        if(st.st_size % kmer_size_in_bytes) {
            std::fprintf(stderr, "File size is not a multiple of 8; incorrectly-sized data");
            std::abort();
        }
        auto ptr = (uint64_t *)std::malloc(st.st_size);
        auto nelem = st.st_size / kmer_size_in_bytes;
        if(std::fread(ptr, kmer_size_in_bytes, nelem, fp) != nelem) {
            std::fprintf(stderr, "Error reading from file %s\n", names[i].data());
            std::abort();
        }
        data[i] = ptr;
        sizes[i] = nelem;
    }
    if(basename) {
        std::transform(names.begin(), names.end(), names.begin(),
        [](std::string x) {
            auto pos = x.find_last_of('/');
            if(pos != std::string::npos) return std::string(&x[pos + 1], &x[x.size()]);
            return x;
        });
    }
    const size_t total_samples = std::accumulate(sizes.begin(), sizes.end(), size_t(0));
    auto start = std::chrono::high_resolution_clock::now();
    std::fprintf(stderr, "Sketches for %zu sequence files created in %fms\n", nfiles, std::chrono::duration<double, std::milli>(start - startload).count());
    // Build map
    ska::flat_hash_map<uint64_t, std::vector<uint32_t>> main_map;
    main_map.reserve(total_samples >> 2); // At most two resizes, perhaps fewer if the collection is highly redundant
    for(uint32_t id = 0; id < nfiles; ++id) {
        uint64_t *const eptr = data[id] + sizes[id];
        for(uint64_t *ptr = data[id]; ptr < eptr; ++ptr) {
            const uint64_t kmer = *ptr++;
            auto it = main_map.find(kmer);
            if(it == main_map.end())
                main_map.emplace(kmer, std::vector<uint32_t>{id}).first;
            else
                it->second.push_back(id);
        }
    }
    auto stop = std::chrono::high_resolution_clock::now();
    std::fprintf(stderr, "Index for %zu sequence files created in %fms, %zu (%%%g possible) total k-mers in the index\n", nfiles, std::chrono::duration<double, std::milli>(stop - start).count(),
                 main_map.size(), 100. * main_map.size() / total_samples);
    uint64_t total_vals  = main_map.size();
    std::vector<uint32_t> nids_per_kmer(total_vals);
    {
        auto mit = main_map.begin();
        for(auto vit = nids_per_kmer.begin(), vie = nids_per_kmer.end(); vit != vie; ++mit, ++vit)
            *vit = mit->second.size();
    }
    const uint64_t total_ids = std::accumulate(main_map.begin(), main_map.end(), uint64_t(0), [](uint64_t x, auto &pair) {return x + pair.second.size();});
    {
        uint64_t arr [] {total_vals, total_ids};
        if(std::fwrite(arr, 16, 1, ofp) != 1)
            goto write_error;
        if(std::fwrite(nids_per_kmer.data(), sizeof(uint32_t), nids_per_kmer.size(), ofp) != nids_per_kmer.size())
            goto write_error;
        // Write keys
        for(const auto &pair: main_map)
            if(std::fwrite(&pair.first, sizeof(pair.first), 1, ofp) != 1)
                goto write_error; 
        for(const auto &pair: main_map)
            if(std::fwrite(pair.second.data(), sizeof(pair.second[0]), pair.second.size(), ofp) != pair.second.size())
                goto write_error; 
        if(0) {
             write_error: {
                 std::fprintf(stderr, "Failed to write to file");
                 ret = 1;
             }
        }
    }
    if(ofp != stdout) std::fclose(ofp);
    for(auto p: data) std::free(p);
    return ret;
}
