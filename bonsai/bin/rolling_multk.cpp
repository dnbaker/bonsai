#include "bonsai/include/kmer_counter.h"
#include "bonsai/include/util.h"
#include <thread>
#include <algorithm>
#include <omp.h>

using namespace kmerc;

void usage() {
    std::fprintf(stderr, "rolling_multk <opts> in.fa\n-P: set prefix\n-p: set number of threads\n-k: Add kmer length\n-C: Do not canonicalize\n-r:  set START,END for kmer range (e.g., -r34,38 will use 34, 35, 36, 37]).\n");
    std::fprintf(stderr, "-K: do not write kvmap\n-S: do not write SHS\n");
    std::exit(1);
}

int main(int argc, char *argv[]) {
    int c;
    std::vector<uint32_t> ks;
    std::string prefix;
    int nt = 1;
    bool canon = true;
    int write_flags = WRITE_KVMAP | WRITE_SHS;
    while((c = getopt(argc, argv, "KSChP:p:k:r:")) >= 0) {
        switch(c) {
            case 'h': usage(); break;
            case 'k': {auto i = std::atoi(optarg); if(i > 0) ks.push_back(i);} break;
            case 'C': canon = false; break;
            case 'P': prefix = optarg; break;
            case 'p': nt = std::atoi(optarg); break;
            case 'K': write_flags &= ~WRITE_KVMAP; break;
            case 'S': write_flags &= ~WRITE_SHS; break;
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
    omp_set_num_threads(nt ? nt: int(std::thread::hardware_concurrency()));
    size_t presize = bns::filesize(argv[optind]);
    if(optind == argc) usage();
    if(prefix.empty()) prefix = argv[optind];
    dump_maps(prefix.data(), ks, argv[optind], canon, presize, write_flags);
    return EXIT_SUCCESS;
}
