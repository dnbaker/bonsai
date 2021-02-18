#include "bonsai/encoder.h"
#include "bonsai/util.h"
#include "kseq_declare.h"
#include "hll/flat_hash_map/flat_hash_map.hpp"
#include <getopt.h>
#ifdef _OPENMP
#include <omp.h>
#endif

using namespace bns;

#ifndef VALUE_TYPE
#define VALUE_TYPE uint32_t
#endif
using CType = ska::flat_hash_map<uint64_t, VALUE_TYPE>;

void update_kmerc(CType &kmerc, const std::string &path, int k, bool canon, const int htype, kseq_t *kseq=static_cast<kseq_t*>(nullptr), RollingHashingType rht=RollingHashingType::DNA) {
    Encoder<> enc(k, canon);
    RollingHasher<uint64_t> rolling_hasher(k, canon, rht);
    auto update_fn = [&kmerc](uint64_t x) {
        auto it = kmerc.find(x);
        if(it == kmerc.end()) kmerc.emplace(x, 1u);
        else ++it->second;
    };
    if(htype == 0) {
        enc.for_each(update_fn, path.data(), kseq);
    } else if(htype == 1) {
        rolling_hasher.for_each_hash(update_fn, path.data(), kseq);
    } else if(htype == 2) {
        enc.for_each_hash(update_fn, path.data(), kseq);
    } else {
        std::fprintf(stderr, "Warning: this should never happen\n");
    }
}

template<typename T, typename Func>
void par_reduce(T *x, size_t n, const Func &func) {
    const unsigned int ln = static_cast<int>(std::ceil(std::log2(n)));
    for(size_t i = 0; i < ln; ++i) {
        const size_t step_size = 1 << i;
        const size_t sweep_size = (i + 1) << 1;
        const size_t nsweeps = (n + (sweep_size - 1)) / sweep_size;
        OMP_PFOR
        for(size_t j = 0; j < nsweeps; ++j) {
            const auto lh = j * sweep_size, rh = lh + step_size;
            if(rh < n)
                func(x[lh], x[rh]);
        }
    }
}


void usage() {
    std::fprintf(stderr, "Usage: kmercnt <opts> (input paths) \n"
                        " Opts: [-k [31] -o [/dev/stdout] -C [true] -N [use_nthash] -c [use_cyclic_hash] -p [1] -S [sort by hash]\n"
                        "-k: set k-mer length\n"
                        "-C: do not canonicalize\n"
                        "-p: set number of threads\n"
                        "-S: sort by hash instead of count\n"
                        "-F: Load paths from <file> (in addition to positional arguments)\n"
                        "-P: Parse protein k-mers instead of DNA k-mers\n"
                        "-B: emit binary (sparse vector) notation\n"
                        "-b: emit binary stream of [uint64_t, uint64_t] k-mer/count pairs\n"
        );
}

int main(int argc, char **argv) {
    std::string ofile = "/dev/stdout", fpaths;
    std::string kmerparsetype = "bns";
    bool canon = true, sort_by_hash = false, enable_protein = false;
    int binary_output = false;
    int k = 31, nthreads = 1;
    for(int c;(c = getopt(argc, argv, "k:F:o:p:bBcCNh?")) >= 0;) {
        switch(c) {
            case 'k': k = std::atoi(optarg); break;
            case 'o': ofile = optarg; break;
            case 'h': usage(); return EXIT_FAILURE;
            case 'N': kmerparsetype = "nthash"; break;
            case 'C': canon = false; break;
            case 'c': kmerparsetype = "cyclic"; break;
            case 'F': fpaths = optarg; break;
            case 'p': nthreads = std::atoi(optarg); break;
            case 'S': sort_by_hash = true; break;
            case 'P': enable_protein = true; kmerparsetype = "cyclic"; break;
            case 'B': binary_output = true; break;
            case 'b': binary_output = 2; break;
            //case 'S': spacestr = optarg; break;
        }
    }
    std::vector<std::string> infiles(argv + optind, argv + argc);
    if(fpaths.size()) {
        std::ifstream ifs(fpaths);
        for(std::string line;std::getline(ifs, line);) {
            infiles.emplace_back(line);
        }
    }
    if(infiles.empty()) {
        infiles.emplace_back("/dev/stdin");
    }
    if(k > 32 && kmerparsetype == "bns") {
        kmerparsetype = "cyclic";
    }
    nthreads = std::max(nthreads, 1);
#ifdef _OPENMP
    omp_set_num_threads(nthreads);
#endif
    std::fprintf(stderr, "Counting %s %u-mers, %s\n", canon ? "canonical": "stranded", k, &kmerparsetype[0]);
    std::vector<kseq_t> kseqs;
    std::vector<CType> threadkmercs(nthreads);
    for(auto &kmerc: threadkmercs) kmerc.reserve(1<<22); // reserve 4MB to start
    while(std::ptrdiff_t(kseqs.size()) < nthreads) kseqs.emplace_back(kseq_init_stack());
    const size_t initsize = 1ull << 22;
    for(auto &ks: kseqs) {
        ks.seq.m = initsize;
        ks.seq.s = static_cast<char *>(std::malloc(initsize));
    }
    const int htype = kmerparsetype == "bns" ? 0: kmerparsetype == "cyclic"? 1: 2;
#ifdef _OPENMP
    #pragma omp parallel for schedule(dynamic)
#endif
    for(size_t i = 0; i < infiles.size(); ++i) {
        std::fprintf(stderr, "Reading from infile %zu (%s)\n", i, infiles[i].data());
        int tid = 0;
#ifdef _OPENMP
        tid = omp_get_thread_num();
#endif
        assert(tid < threadkmercs.size());
        update_kmerc(threadkmercs.at(tid), infiles[i], k, canon, htype, &kseqs[tid], enable_protein ? bns::RollingHashingType::PROTEIN: RollingHashingType::DNA);
    }
    auto &kmerc = threadkmercs.front();
    par_reduce(threadkmercs.data(), threadkmercs.size(), [](CType &lhs, const CType &rhs) {
        for(const auto &pair: rhs) {
            auto it = lhs.find(pair.first);
            if(it != lhs.end()) it->second += pair.second;
            else lhs.emplace(pair);
        }
    });
    std::vector<std::pair<uint64_t, VALUE_TYPE>> res(kmerc.size());
    std::copy(kmerc.begin(), kmerc.end(), res.begin());
    std::ios_base::sync_with_stdio(false);
    const auto maxc = std::max_element(res.begin(), res.end(), [](auto x, auto y) {return x.second < y.second;})->second;
    std::fprintf(stderr, "maxcount: %zu\n", size_t(maxc));
    if(sort_by_hash) {
        std::sort(res.begin(), res.end());
    } else {
        std::sort(res.begin(), res.end(), [](auto x, auto y) {return x.second > y.second;});
    }
    std::fprintf(stderr, "Counted %zu unique k-mers with max count %zu\n", res.size(), size_t(maxc));
    if(binary_output == 1) {
        if(ofile == "/dev/stdout") throw std::runtime_error("Must provide a path for binary output");
        std::FILE *ofb = std::fopen((ofile + "hash").data(), "wb");
        std::FILE *ofc = std::fopen((ofile + "count").data(), "wb");
        for(const auto &pair: res) {
            uint64_t h = pair.first;
            uint32_t uv;
            uint16_t us;
            uint8_t ub;
            std::fwrite(&h, 1, sizeof(h), ofb);
            if(sizeof(VALUE_TYPE) > 4 && maxc > 0xFFFFFFFFull) {h = pair.second;std::fwrite(&h, 1, sizeof(h), ofc);}
            else if(maxc > 0xFFFFull) {uv = pair.second;std::fwrite(&uv, 1, sizeof(uv), ofc);}
            else if(maxc > 0xFFull) {us = pair.second;std::fwrite(&us, 1, sizeof(us), ofc);}
            else {ub = pair.second; std::fwrite(&ub, 1, sizeof(ub), ofc);}
        }
        std::fclose(ofc); std::fclose(ofb);
    } else if(binary_output == 2) {
        std::FILE *ofp = std::fopen(ofile.data(), "wb");
        uint64_t x[2];
        for(const auto &pair: res) {
            x[0] = pair.first, x[1] = pair.second;
            std::fwrite(x, 1, sizeof(x), ofp);
        }
        std::fclose(ofp);
    } else {
        std::ofstream ofs(ofile);
        ofs << "ID\tCount\n";
        for(const auto &pair: res) {
            ofs << pair.first << '\t' << pair.second << '\n';
        }
    }
}
