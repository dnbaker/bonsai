#include "bonsai/include/kmeridx.h"
#include <getopt.h>

template<typename T1, typename T2>
void build_write(unsigned k, const char *s, const char *o) {
    bns::KmerIdx<T1, T2>(k, s).write(o);
}

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
    auto fptr = pos32 ? (small_idx ? build_write<uint32_t, uint32_t>: build_write<uint32_t, uint64_t>)
                      : (small_idx ? build_write<uint64_t, uint32_t>: build_write<uint64_t, uint64_t>);
    fptr(k,argv[optind], argv[optind + 1] ? const_cast<const char *>(argv[optind + 1]): "/dev/stdout");
}
