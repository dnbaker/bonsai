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
#include "bonsai/encoder.h"

#if _OPENMP
#define OMP_ELSE(x, y) x
#define OMP_ONLY(...) __VA_ARGS__
#else
# define OMP_ELSE(x, y) y
#define OMP_ONLY(...)
#endif

using namespace bns::lsh;
using namespace bns;

#include "hll/include/flat_hash_map/flat_hash_map.hpp"
using MapT = ska::flat_hash_map<uint64_t, uint32_t>;

using std::uint64_t;
using std::uint32_t;

int usage() {
    std::fprintf(stderr, "Usage: setsketchscreener <opts> setsketch.db [input sequence files..]\n"
                         "-c: cyclic hash\n"
                         "-N: use cyclic hash\n"
                         "-C: Do not canonicalize\n"
                         "-P: Enable protein encoding (implies cylic)\n"
    );
    return 1;
}

template<typename Func>
void for_each_kmer(Encoder<> &enc, RollingHasher<uint64_t> &rolling_hasher, const std::string &path, const int htype, const Func &func, kseq_t *kseq=static_cast<kseq_t*>(nullptr)) {
    if(htype == 0) {
        enc.for_each(func, path.data(), kseq);
    } else if(htype == 1) {
        rolling_hasher.for_each_hash(func, path.data(), kseq);
    } else if(htype == 2) {
        enc.for_each_hash(func, path.data(), kseq);
    } else {
        std::fprintf(stderr, "Error: this should never happen. htype should be [0, 1, 2]\n");
        std::exit(EXIT_FAILURE);
    }
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

MapT process(bns::Encoder<> *encoders, bns::RollingHasher<uint64_t> *rencoders,
             kseq_t *kseqs, std::vector<MapT> &matchcounts, int htype, const std::string &fn, int nthreads,
                std::vector<std::pair<char *, size_t>> &seqbuffers);

int main(int argc, char **argv) {
    bool basename = false;
    int ret = 0;
    std::vector<std::string> names;
    std::FILE *ofp = stdout;
    std::string kmerparsetype = "bns";
    bool enable_protein = false;
    bool canon = true;
    int k = -1, nthreads = 1;
    const size_t initsize = 2000000;
    for(int c;(c = getopt(argc, argv, "pPCcNk:o:F:bh")) >= 0;) switch(c) {
        case 'p': nthreads = std::atoi(optarg); break;
        case 'b': basename = true; break;
        case 'c': kmerparsetype = "cyclic"; break;
        case 'C': canon = true; break;
        case 'h': case '?': return usage();
        case 'o': if(!(ofp = std::fopen(optarg, "w"))) {std::fprintf(stderr, "Could not open file at %s\n", optarg); std::abort();}
                  break;
        case 'k': k = std::atoi(optarg); break;
        case 'P': enable_protein = true; kmerparsetype = "cyclic"; break;
        case 'N': kmerparsetype = "nthash"; break;
        case 'F': {
            std::ifstream ifs(optarg);
            for(std::string s; std::getline(ifs, s);)
                names.push_back(s);
        }
    }
    if(nthreads <= 0) nthreads = 1;
    const int htype = kmerparsetype == "bns" ? 0: kmerparsetype == "cyclic"? 1: 2;
    if(argc == optind) throw 1;
    else if(argc == optind + 1) {
        if(names.empty()) throw std::runtime_error("no query paths provided");
    } else
        names.insert(names.end(), argv + optind + 1, argv + argc);
    auto mapk = read_database(argv[optind]);
    auto map = std::move(mapk.first);
    auto loaded_k = mapk.second;
    if(k < 0) {
        if(loaded_k <= 0)
            throw std::invalid_argument("k must be provided by command-line or at db construction.");
    } else k = loaded_k;
    MapT counter;
    counter.reserve(map.size());
    for(const auto &pair: map) counter.emplace(pair.first, 0u);
    std::fprintf(stderr, "map size %zu and total number of ids %zu\n", map.size(), std::accumulate(map.begin(), map.end(), size_t(0), [](size_t x, auto &y) {return x + y.second.size();}));
    if(names.empty()) return usage();
    const RollingHashingType alphabet = enable_protein ? RollingHashingType::PROTEIN: RollingHashingType::DNA;
    const size_t memtoalloc = (sizeof(bns::RollingHasher<uint64_t>) + sizeof(bns::Encoder<>)) * nthreads + 63;
    std::unique_ptr<uint8_t> mem(new uint8_t[memtoalloc]);
    uint8_t *memp = mem.get(), *ap = reinterpret_cast<uint8_t *>(reinterpret_cast<uint64_t>(memp) + (uint64_t)memp % 64 ? int(64 - (uint64_t)memp % 64): 0);
    bns::RollingHasher<uint64_t> *rencoders = (bns::RollingHasher<uint64_t> *)ap;
    bns::Encoder<> *encoders = (bns::Encoder<> *)(&rencoders[nthreads]);
    kseq_t *kseqs = static_cast<kseq_t *>(std::calloc(nthreads, sizeof(kseq_t)));
    std::vector<MapT> matchcounts(nthreads);
    std::vector<std::pair<char *, size_t>> seqbuffers;
    OMP_PFOR
    for(int idx = 0; idx < nthreads; ++idx) {
        auto &back = kseqs[idx];
        back.seq.m = initsize;
        back.seq.s = static_cast<char *>(std::malloc(initsize));
        new (encoders + idx) Encoder<>(k, canon);
        new (rencoders + idx) RollingHasher<uint64_t>(k, canon, alphabet);
    }
    std::vector<MapT> results;
    results.reserve(names.size());
    for(const auto &fn: names) {
        results.push_back(process(encoders, rencoders, kseqs, matchcounts, htype, fn, nthreads, seqbuffers));
        //for_each_kmer(encoder, renc, fn.data(), htype,
    }
    OMP_PFOR
    for(int i = 0; i < nthreads; ++i) {
        kseq_destroy_stack(kseqs[i]);
        encoders[i].~Encoder<>();
        rencoders[i].~RollingHasher<uint64_t>();
    }
    std::free(kseqs);
    return 0;
}

MapT process(bns::Encoder<> *encoders, bns::RollingHasher<uint64_t> *rencoders,
             kseq_t *kseqs, std::vector<MapT> &matchcounts, int htype, const std::string &fn, int nthreads,
                std::vector<std::pair<char *, size_t>> &seqbuffers)
{
#if 0
    for(auto &m: matchcounts) m.clear();
    for(auto &pair: seqbuffers) std::free(seqbuffers.first), seqbuffers.first = 0;
    gzFile ifp = gzopen(fn.data(), "rb");
    kseq_assign(kseqs, ifp);
    for(;kseq_read(kseqs) >= 0;) {
        std::pair<char *, size_t> tup;
        tup.first = (char *)std::malloc(kseqs->seq.l + 1);
        std::memcpy(tup.first, kseqs->seq.s, kseqs->seq.l + 1);
        tup.second = kseqs->seq.l;
    }
    OMP_PFOR
    for(size_t i = 0; i < seqbuffers.size(); ++i) {
    }
    gzclose(ifp);
    par_reduce(matchcounts.data(), matchcounts.size());
    return matchcounts.front();  // Copy is implied
#endif
    return MapT();
}
