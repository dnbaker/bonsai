#include "hll/sketch.h"
#include "bonsai/include/encoder.h"
#include "omp.h"

void usage() {
    std::fprintf(stderr, "rolling_multk_sketch <opts> in.fa\n-P: set prefix\n-p: set number of threads\n-k: Add kmer length\n-C: Do not canonicalize\n-r:  set START,END for kmer range (e.g., -r34,38 will use 34, 35, 36, 37]).\n");
    std::fprintf(stderr, "-S: set log 2 sketch size (14)\n");
    std::exit(1);
}

template<typename C, typename Sketch=sketch::hll_t, typename IT=uint64_t, typename ArgType, typename ... Args>
std::vector<Sketch> build_multk_sketches(const C &kmer_sizes, ArgType fp, bool canon=false, Args &&... args) {
    static_assert(std::is_same<ArgType, gzFile>::value  || std::is_same<ArgType, char *>::value || std::is_same<ArgType, const char *>::value, "Must be gzFile, char *, or const char *");
    bns::RollingHasherSet<IT> rhs(kmer_sizes, canon);
    const size_t nsk = kmer_sizes.size();
    std::vector<Sketch> sketches;
    sketches.reserve(nsk);
    while(sketches.size() < kmer_sizes.size())
        sketches.emplace_back(std::forward<Args>(args)...);
    rhs.for_each_hash([&sketches](IT hashvalue, size_t idx) {
        auto &sketch = sketches[idx];
        sketch.add(hashvalue);
    }, fp);
    return sketches;
}


int main(int argc, char *argv[]) {
    int c;
    std::string prefix;
    int canon = true;
    size_t l2sz = 14;
    std::vector<uint32_t> ks;
    while((c = getopt(argc, argv, "ChP:p:k:r:S:")) >= 0) {
        switch(c) {
            case 'h': usage(); break;
            case 'k': {auto i = std::atoi(optarg); if(i > 0) ks.push_back(i);} break;
            case 'C': canon = false; break;
            case 'P': prefix = optarg; break;
            case 'S': l2sz = std::atoi(optarg); break;
            case 'p': omp_set_num_threads(std::atoi(optarg)); break;
            case 'r':
                {
                    auto s = std::atoi(optarg);
                    auto p = std::strchr(optarg, ',');
                    auto step = p ? std::strchr(p + 1, ',') ? std::atoi(std::strchr(p + 1, ',') + 1): 1: 1;
                    int e = p ? std::atoi(p + 1): s + 1;
                    ks.clear();
                    while(s < e)
                        ks.push_back(s), s += step;
                }
        }
    }
    if(optind == argc) usage();
    if(prefix.empty()) prefix = argv[optind];
    gzFile fp = gzopen(argv[optind], "rb");
    auto sketches = build_multk_sketches(ks, fp, canon, l2sz);
    for(const auto &s: sketches)
        s.write((prefix + "." + std::to_string(ks[&s - &*sketches.begin()]) + ".sketch." + std::to_string(l2sz) + ".hll").data());
}
