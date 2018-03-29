#include <getopt.h>
#include "include/kfreq.h"

using namespace emp;

void usage(char **argv) {
    std::fprintf(stderr, "Usage: %s [max_kmer_size] [genome1] [genome2] ...\n", *argv);
    std::exit(EXIT_FAILURE);
}

std::string canonicalize(const char *str) {
    if(const char *slash = std::strrchr(str, '/'))
        return slash + 1;
    return str;
}

int main(int argc, char *argv[]) {
    if(argc == 1) usage(argv);
    std::vector<std::string> paths;
    unsigned ks = std::atoi(argv[1]);
    if(ks > 16 || ks < 2) {
        std::fprintf(stderr, "ks: %u. Max supported: 16. Min: 2\n", ks);
        usage(argv);
    }
    for(char **p(argv + 2); *p; paths.emplace_back(*p++));
    kseq_t seq = kseq_init_stack();
    freq::KFC kc(ks);
    for(const auto &path: paths) {
        std::string outpath = canonicalize(path.data()) + ".k" + std::to_string(ks) + ".bin";
        kc.add(path.data(), &seq);
        kc.write(outpath.data());
        kc.clear();
    }
    kseq_destroy_stack(seq);
}
