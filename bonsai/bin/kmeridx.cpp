#include "bonsai/include/kmeridx.h"
#include <getopt.h>

int main(int argc, char *argv[]) {
    int c;
    bool pos32 = false;
    unsigned k = 31;
    if(argc == 1) goto help;
    while((c = getopt(argc, argv, "k:3:h")) >= 0) {
        switch(c) {
            case '3':  pos32 = true; break;
            case 'k': k = std::atoi(optarg); break;
            case 'h': help: std::fprintf(stderr, "Usage: <executable> [opts] path.fa out.index.gz\nFlags:\n\n-k\tSet kmer [31]\n-3\tUse 32-bit positions [64-bit]\n"); std::exit(0);
        }
    }
    if(k > 32) RUNTIME_ERROR("k must be <= 32");
    bool small_idx = k < 17;
    if(pos32) {
        if(small_idx) {
            bns::KmerIdx<uint32_t, uint32_t> idx(k, argv[optind]);
            idx.write(argv[optind + 1]);
        } else {
            bns::KmerIdx<uint32_t, uint64_t> idx(k, argv[optind]);
            idx.write(argv[optind + 1]);
        }
    }
    else {
        if(small_idx) {
            bns::KmerIdx<uint64_t, uint32_t> idx(k, argv[optind]);
            idx.write(argv[optind + 1]);
        } else {
            bns::KmerIdx<uint64_t, uint64_t> idx(k, argv[optind]);
            idx.write(argv[optind + 1]);
        }
    }
}
