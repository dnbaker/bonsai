#define NO_BLAZE

#include <getopt.h>
#include <unistd.h>
#include <sys/stat.h>

#include <cstdio>
#include <memory>
#include <cstdint>
#include <vector>
#include <chrono>
#include <string>
#include <iostream>
#include <fstream>
#include <numeric>
#include <cassert>
#include <thread>

#ifdef _OPENMP
#include <omp.h>
#endif
#include <zlib.h>
#include "bonsai/ssi.h"

#if _OPENMP
#define OMP_ELSE(x, y) x
#define OMP_ONLY(...) __VA_ARGS__
#else
# define OMP_ELSE(x, y) y
#define OMP_ONLY(...)
#endif

using std::uint64_t;
using std::uint32_t;
using namespace bns::lsh;

int usage() {
    std::fprintf(stderr, "setsketchroller <flags> sketch1.kmers sketch2.kmers...\n"
                         "-b\tTrim folder path from filenames\n"
                         "-F\tRead paths from <file> instesad of positional arguments\n"
                         "-k\tProvide k for database.\n"
                         "-o\tWrite database to file instead of stdout\n");
    return 1;
}

template<typename T>
void par_reduce(T *x, size_t n) {
    const unsigned int ln = static_cast<int>(std::ceil(std::log2(n)));
    for(size_t i = 0; i < ln; ++i) {
        const size_t step_size = 1 << i;
        const size_t sweep_size = (i + 1) << 1;
        const size_t nsweeps = (n + (sweep_size - 1)) / sweep_size;
        OMP_ONLY(_Pragma("omp parallel for"))
        for(size_t j = 0; j < nsweeps; ++j) {
            const auto lh = j * sweep_size, rh = lh + step_size;
            if(rh < n)
                x[lh] += x[rh];
        }
    }
}

// Step 1: load k-mer files
// Step 2: invert matrix

int main(int argc, char **argv) {
    bool basename = false;
    int ret = 0, nthreads = 1;
    std::vector<std::string> names, qnames;
    uint64_t kmer_size_in_bytes = 8;
    std::FILE *ofp = stdout;
    std::string ofname;
    int k = -1;
    for(int c;(c = getopt(argc, argv, "k:p:o:F:bh")) >= 0;) switch(c) {
        case 'b': basename = true; break;
        case 'h': case '?': return usage();
        case 'o': ofname = optarg; if(!(ofp = std::fopen(optarg, "w"))) {std::fprintf(stderr, "Could not open file at %s\n", optarg); std::abort();}
                  break;
        case 'p': nthreads = std::atoi(optarg); break;
        case 'k': k = std::atoi(optarg); break;
        case 'F': {
            std::ifstream ifs(optarg);
            for(std::string s; std::getline(ifs, s);)
                names.push_back(s);
            break;
        }
        case 'q': {
            qnames.push_back(optarg);
            break;
        }
    }
    if(nthreads < 0) nthreads = std::thread::hardware_concurrency();
    else if(nthreads == 0) nthreads = 1;
    omp_set_num_threads(nthreads);
    for(char **p = argv + optind;p < argv + argc; names.emplace_back(*p++));
    if(names.empty()) return usage();
    const size_t nfiles = names.size();
    std::vector<uint64_t *>data(nfiles);
    std::vector<uint64_t> sizes(nfiles);
    // Load files
    auto startload = std::chrono::high_resolution_clock::now();
    OMP_ONLY(_Pragma("omp parallel for"))
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
    std::vector<ska::flat_hash_map<uint64_t, std::vector<uint32_t>>> maps;
    while(maps.size() < uint32_t(nthreads)) {
        maps.emplace_back();
    }
    assert(maps.size());
    auto &main_map = maps.front();
    OMP_ONLY(_Pragma("omp parallel for"))
    for(uint32_t id = 0; id < nfiles; ++id) {
        const int tid = OMP_ELSE(omp_get_thread_num(), 0);
        auto &map = maps[tid];
        uint64_t *const eptr = data[id] + sizes[id];
        for(uint64_t *ptr = data[id]; ptr < eptr; ++ptr) {
            const uint64_t kmer = *ptr++;
            auto it = map.find(kmer);
            if(it == map.end())
                map.emplace(kmer, std::vector<uint32_t>{id}).first;
            else
                it->second.push_back(id);
        }
    }
    par_reduce(maps.data(), maps.size());
    auto stop = std::chrono::high_resolution_clock::now();
    std::fprintf(stderr, "Index for %zu sequence files created in %fms, %zu (%%%g possible) total k-mers in the index\n", nfiles, std::chrono::duration<double, std::milli>(stop - start).count(),
                 main_map.size(), 100. * main_map.size() / total_samples);
    if(k < 0) std::fprintf(stderr, "Warning: k not provided. This will need to be passed at query time for accurate processing.\n");
    ret = bns::lsh::write_database(main_map, ofp, k);
    if(ofp != stdout) std::fclose(ofp);
    for(auto p: data) std::free(p);
#if !defined(NDEBUG)
    if(ofp != stdout) {
        std::FILE *tfp = std::fopen(ofname.data(), "rb");
        auto matcpypair = bns::lsh::read_database(tfp);
        assert(matcpypair.second == k);
        auto matcpy = std::move(matcpypair.first);
        std::fclose(tfp);
        std::fprintf(stderr, "Sizes: %zu/%zu\n", matcpy.size(), main_map.size());
        assert(matcpy == main_map);
    }
#endif
    return ret;
}

