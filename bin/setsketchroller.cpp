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
#include "flat_hash_map/flat_hash_map.hpp"

using std::uint64_t;
using std::uint32_t;

int usage() {
    std::fprintf(stderr, "setsketchfolder <flags> sketch1.kmers sketch2.kmers...\n"
                         "-b\tTrim folder path from filenames\n"
                         "-F\tRead paths from <file> instesad of positional arguments\n");
    return 1;
}

template<typename T>
void par_reduce(T *x, size_t n) {
    const unsigned int ln = static_cast<int>(std::ceil(std::log2(n)));
    for(size_t i = 0; i < ln; ++i) {
        const size_t step_size = 1 << i;
        const size_t sweep_size = (i + 1) << 1;
        const size_t nsweeps = (n + (sweep_size - 1)) / sweep_size;
        #if _OPENMP
        #pragma omp parallel for
        #endif
        for(size_t j = 0; j < nsweeps; ++j) {
            const auto lh = j * sweep_size, rh = lh + step_size;
            if(rh < n)
                x[lh] += x[rh];
        }
    }
}

ska::flat_hash_map<uint64_t, std::vector<uint32_t>> &operator+=(ska::flat_hash_map<uint64_t, std::vector<uint32_t>> &lhs, const ska::flat_hash_map<uint64_t, std::vector<uint32_t>> &rhs) {
    for(auto &&pair: rhs) {
        auto it = lhs.find(pair.first);
        if(it == lhs.end()) it = lhs.emplace(pair).first;
        else it->second.insert(it->second.end(), pair.second.begin(), pair.second.end());
    }
    return lhs;
}


// You can use this database for inverted screen sweeps
ska::flat_hash_map<uint64_t, std::vector<uint32_t>>
read_file(std::FILE *fp) {
    auto timestart = std::chrono::high_resolution_clock::now();
    uint64_t arr[2];
    uint64_t fret;
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
        //if(i % 256 == 0) std::fprintf(stderr, "%zu/%zu, read %zu\n", i + 1, size_t(arr[1]), total_ids_read);
        const auto nids = data[i];
        buffer.resize(nids);
        if((fret = std::fread(buffer.data(), sizeof(uint32_t), nids, fp)) != nids) {
            throw std::runtime_error(std::string("Expected to read ") + std::to_string(nids) + "-block of u32s and got " + std::to_string(fret));
        }
        total_ids_read += nids;
        map[keys[i]] = buffer;
        map.emplace(keys[i], buffer);
    }
    std::fprintf(stderr, "Time to deserialize: %gms\n", std::chrono::duration<double, std::milli>(std::chrono::high_resolution_clock::now() - timestart).count());
    return map;
}
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
// Step 1: load k-mer files
// Step 2: invert matrix

int main(int argc, char **argv) {
    bool basename = false;
    int ret = 0, nthreads = 1;
    std::vector<std::string> names;
    uint64_t kmer_size_in_bytes = 8;
    std::FILE *ofp = stdout;
    std::string ofname;
    for(int c;(c = getopt(argc, argv, "p:o:F:bh")) >= 0;) switch(c) {
        case 'b': basename = true; break;
        case 'h': case '?': return usage();
        case 'o': ofname = optarg; if(!(ofp = std::fopen(optarg, "w"))) {std::fprintf(stderr, "Could not open file at %s\n", optarg); std::abort();}
                  break;
        case 'p': nthreads = std::atoi(optarg); break;
        case 'F': {
            std::ifstream ifs(optarg);
            for(std::string s; std::getline(ifs, s);)
                names.push_back(s);
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
    std::vector<ska::flat_hash_map<uint64_t, std::vector<uint32_t>>> maps;
    while(maps.size() < uint32_t(nthreads)) {
        maps.emplace_back();
    }
    assert(maps.size());
    auto &main_map = maps.front();
#if _OPENMP
    #pragma omp parallel for
#endif
    for(uint32_t id = 0; id < nfiles; ++id) {
        const int tid =
#if _OPENMP
                omp_get_thread_num();
#else
                0;
#endif
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
        //std::fprintf(stderr, "After %u/%zu, %zu total keys\n", id + 1, nfiles, map.size());
    }
    par_reduce(maps.data(), maps.size());
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
            if(std::fwrite(&pair.first, sizeof(uint64_t), 1, ofp) != 1)
                goto write_error; 
        static_assert(sizeof(main_map.begin()->first) == sizeof(uint64_t), "sanity check");
        for(const auto &pair: main_map) {
            if(std::fwrite(pair.second.data(), sizeof(uint32_t), pair.second.size(), ofp) != pair.second.size())
                goto write_error; 
        }
        if(0) {
             write_error: {
                 std::fprintf(stderr, "Failed to write to file");
                 ret = 1;
             }
        }
    }
    if(ofp != stdout) std::fclose(ofp);
    for(auto p: data) std::free(p);
#if !defined(NDEBUG)
    if(ofp != stdout) {
        std::FILE *tfp = std::fopen(ofname.data(), "rb");
        auto matcpy = read_file(tfp);
        std::fclose(tfp);
        std::fprintf(stderr, "Sizes: %zu/%zu\n", matcpy.size(), main_map.size());
        assert(matcpy == main_map);
    }
#endif
    return ret;
}
