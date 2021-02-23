#include "bonsai/encoder.h"
#include "bonsai/util.h"
#include "kseq_declare.h"
#include "hll/flat_hash_map/flat_hash_map.hpp"
#include <getopt.h>
#ifdef _OPENMP
#include <omp.h>
#endif
#include "sketch/setsketch.h"

using namespace bns;

#ifndef VALUE_TYPE
#define VALUE_TYPE uint32_t
#endif

template<typename Sketch>
void update_sketch(Encoder<> &enc, RollingHasher<uint64_t> &rolling_hasher, Sketch &sketch, const std::string &path, const int htype, kseq_t *kseq=static_cast<kseq_t*>(nullptr)) {
    auto update_fn = [&sketch](uint64_t x) {
        sketch.update(x);
    };
    if(htype == 0) {
        enc.for_each(update_fn, path.data(), kseq);
    } else if(htype == 1) {
        rolling_hasher.for_each_hash(update_fn, path.data(), kseq);
    } else if(htype == 2) {
        enc.for_each_hash(update_fn, path.data(), kseq);
    } else {
        std::fprintf(stderr, "Error: this should never happen. htype should be [0, 1, 2]\n");
        std::exit(EXIT_FAILURE);
    }
}

struct PlusEq {
    template<typename T>
    T &operator()(T &lhs, const T &rhs) const {return lhs += rhs;}
};

template<typename T, typename Func=PlusEq>
void par_reduce(T *x, size_t n, const Func &func=Func()) {
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
    std::fprintf(stderr, "Usage: setsketcher <opts> (input paths) \n"
                        " Opts: [-k [31] -o [/dev/stdout] -C [true] -N [use_nthash] -c [use_cyclic_hash] -p [1]] <paths>\n"
                        "-k: set k-mer length\n"
                        "-C: do not canonicalize\n"
                        "-p: set number of threads\n"
                        "-s: save kmers. This yields additional files consisting of sampled k-mers as well. (This can be used to build estimators for GC bias, Shannon entropy, and seeding)\n"
                        "-S: save kmer counts [Implies -s]. This yields additional files describing multiplicities of k-mers.\n"
                        "-z: Set sketch size (default: 4096)\n"
                        "-F: Load paths from <file> (in addition to positional arguments)\n"
                        "-P: Parse protein k-mers instead of DNA k-mers [this implies cyclic, avoiding direct encoding]\n"
                        "-I: set initial buffer size for sequence parsing to [size_t] (4194304 = 4MiB)\n"
        );
}

int main(int argc, char **argv) {
    std::string ofile = "/dev/stdout", fpaths;
    std::string kmerparsetype = "bns";
    bool canon = true, enable_protein = false, save_kmers = 0, save_kmer_counts = 0;
    int k = 31, nthreads = 1;
    size_t initsize = 1ull << 22, sketchsize = 4096;
    for(int c;(c = getopt(argc, argv, "I:k:F:o:p:z:PsScCNh?")) >= 0;) {
        switch(c) {
            case 'k': k = std::atoi(optarg); break;
            case 'o': ofile = optarg; break;
            case 'h': usage(); return EXIT_FAILURE;
            case 'N': kmerparsetype = "nthash"; break;
            case 'C': canon = false; break;
            case 'c': kmerparsetype = "cyclic"; break;
            case 'z': sketchsize = std::strtoull(optarg, nullptr, 10); break;
            case 's': save_kmer_counts = true; [[fallthrough]];
            case 'S': save_kmers = true; break;
            case 'F': fpaths = optarg; break;
            case 'p': nthreads = std::atoi(optarg); break;
            case 'P': enable_protein = true; kmerparsetype = "cyclic"; break;
            case 'I': initsize = std::strtoull(optarg, nullptr, 10); break;
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
    if(enable_protein) {
        if(kmerparsetype != "cyclic") {
            std::fprintf(stderr, "Warning: enable_protein implies cyclic hashing.\n");
            kmerparsetype = "cyclic";
        }
    } else if(k > 32 && kmerparsetype == "bns") {
        std::fprintf(stderr, "Warning: k > 32 implies nthash-based rolling hashing.\n");
        kmerparsetype = "nthash";
    }
    if((save_kmers || save_kmer_counts) && (ofile == "/dev/stdout" || ofile == "-")) {
        std::fprintf(stderr, "Error: -o required if save_kmers or save_kmer_counts is set. (This yields multiple files using -o as a prefix.)");
        std::exit(EXIT_FAILURE);
    }
    nthreads = std::max(nthreads, 1);
#ifdef _OPENMP
    omp_set_num_threads(nthreads);
#endif
    std::fprintf(stderr, "Sketching %s %u-mers, %s\n", canon ? "canonical": "stranded", k, &kmerparsetype[0]);
    std::vector<kseq_t> kseqs;
    std::vector<Encoder<>> encoders;
    std::vector<RollingHasher<uint64_t>> rencoders;
    std::vector<sketch::CSetSketch<double>> sketches;
#define RSV(x) x.reserve(nthreads);
    RSV(kseqs)  RSV(encoders) RSV(rencoders) RSV(sketches)
#undef RSV
    const RollingHashingType rht = enable_protein ? RollingHashingType::PROTEIN: RollingHashingType::DNA;
    while(std::ptrdiff_t(kseqs.size()) < nthreads) {
        kseqs.emplace_back(kseq_init_stack());
        auto &back = kseqs.back();
        back.seq.m = initsize;
        back.seq.s = static_cast<char *>(std::malloc(initsize));
        encoders.emplace_back(k, canon);
        rencoders.emplace_back(k, canon, rht);
        sketches.emplace_back(sketchsize, save_kmers, save_kmer_counts);
    }
    const int htype = kmerparsetype == "bns" ? 0: kmerparsetype == "cyclic"? 1: 2;
#ifdef _OPENMP
    #pragma omp parallel for schedule(dynamic)
#endif
    for(size_t i = 0; i < infiles.size(); ++i) {
        int tid = 0;
#ifdef _OPENMP
        tid = omp_get_thread_num();
#endif
        update_sketch(
                encoders[tid], rencoders[tid], sketches[tid], // Parsing/Sketching prep
                infiles[i], htype, &kseqs[tid]                // Path/Sketch format/buffer
        );
    }
    par_reduce(&sketches[0], nthreads);
    if(!save_kmers) {
        sketches.front().write(ofile);
    } else {
        auto &f = sketches.front();
        f.write(ofile + ".csketch");
        {
            {
                std::ofstream ofs(ofile + ".u64.kmers", std::ios::out | std::ios::binary);
                ofs.write((const char *)f.ids().data(), f.ids().size() * sizeof(f.ids()[0]));
            }
            {
                std::ofstream ofshr(ofile + ".humanreadable.kmers", std::ios::out);
                const auto idp = f.ids().data();
                const auto idcp = f.idcounts().data();
                for(size_t i = 0; i < f.ids().size(); ++i) {
                    uint64_t id = idp[i];
                    if(htype == 0) ofshr << encoders.front().sp_.to_string(id);
                    else           ofshr << id;
                    if(f.idcounts().size()) {
                        ofshr << '\t' << idcp[i];
                    }
                    ofshr << '\n';
                }
            }
            {
                std::ofstream ofs(ofile + ".u32.kmercounts", std::ios::out | std::ios::binary);
                ofs.write((const char *)f.idcounts().data(), f.idcounts().size() * sizeof(f.idcounts()[0]));
            }
        }
    }
    return EXIT_SUCCESS;
}
