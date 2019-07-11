#include "hll/sketch.h"
#include "bonsai/include/encoder.h"
#include "omp.h"

void usage() {
    std::fprintf(stderr, "rolling_multk_sketch <opts> in.fa\n-P: set prefix\n-p: set number of threads\n-k: Add kmer length\n-C: Do not canonicalize\n-r:  set START,END for kmer range (e.g., -r34,38 will use 34, 35, 36, 37]).\n");
    std::fprintf(stderr, "-S: set log 2 sketch size (14)\n");
    std::exit(1);
}

template<typename Sketch>
Sketch &fn2sketch(std::string fn, Sketch &sk) {
    gzFile fp = gzopen(fn.data(), "rb");
    uint64_t w;
    if(gzread(fp, &w, sizeof(w)) != sizeof(w)) throw 1;
    size_t nelem = w;
    size_t count = 0;
    while(gzread(fp, &w, w) == sizeof(w)) {
        sk.add(w);
        ++count;
    }
    if(count != nelem) throw 2;
    return sk;
}


int main(int argc, char *argv[]) {
    int c;
    std::string output;
    size_t l2sz = 14;
    while((c = getopt(argc, argv, "S:h")) >= 0) {
        switch(c) {
            case 'h': usage(); break;
            case 'S': l2sz = std::atoi(optarg); break;
        }
    }
    std::vector<std::string> inputs(argv + optind, argv + argc);
    std::vector<sketch::hll_t> sketches;
    sketches.reserve(inputs.size());
    while(sketches.size() < inputs.size())
        sketches.emplace_back(l2sz);
    #pragma omp parallel for
    for(size_t i = 0; i < inputs.size(); ++i) {
        fn2sketch(inputs[i], sketches[i]);
        sketches[i].write(inputs[i] + ".sketch.hll");
    }
}
