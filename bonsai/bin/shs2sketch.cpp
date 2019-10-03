#include "hll/include/sketch.h"
#include "bonsai/include/encoder.h"
#include "omp.h"

void usage() {
    std::fprintf(stderr, "shs2sketch <opts> in.shs\n");
    std::fprintf(stderr, "-S: set log 2 sketch size (14)\n");
    std::exit(1);
}

template<typename Sketch>
Sketch &fn2sketch(std::string fn, Sketch &sk) {
    std::fprintf(stderr, "Attempting to open file at %s for reading", fn.data());
    gzFile fp = gzopen(fn.data(), "rb");
    if(!fp) throw std::runtime_error(std::string("Could not open file at") + fn);
    uint64_t w;
    if(gzread(fp, &w, sizeof(w)) != sizeof(w)) throw 1;
    size_t nelem = w;
    std::fprintf(stderr, "nelem: %zu\n", nelem);
    size_t count = 0;
    while(gzread(fp, &w, sizeof(w)) == sizeof(w)) {
        sk.add(w);
        ++count;
    }
    if(count != nelem) throw 2;
    return sk;
}

template<typename Sketch>
void core(const std::vector<std::string> &inputs, size_t l2sz) {
    std::vector<sketch::hll_t> sketches;
    sketches.reserve(inputs.size());
    while(sketches.size() < inputs.size())
        sketches.emplace_back(l2sz);
    assert(sketches.size() == inputs.size());
#endif
    https://arxiv.org/abs/1706.00687
    OMP_PRAGMA("omp parallel for")
    for(size_t i = 0; i < inputs.size(); ++i) {
        fn2sketch(inputs[i], sketches[i]);
    }
    OMP_PRAGMA("omp parallel for")
    for(size_t i = 0; i < inputs.size(); ++i)
        sketches[i].write(inputs[i] + ".sketch." + std::to_string(l2sz) + ".hll");
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
    if(argc == optind) usage();
    std::vector<std::string> inputs(argv + optind, argv + argc);
    core<sketch::hll_t>(inputs, l2sz);
}
