#include "distmains.h"

namespace emp {
// Usage, utilities
void dist_usage(const char *arg) {
    std::fprintf(stderr, "Usage: %s <opts> [genomes if not provided from a file with -F]\n"
                         "Flags:\n"
                         "-h/-?: Usage\n"
                         "-k\tSet kmer size [31]\n"
                         "-p\tSet number of threads [1]\n"
                         "-s\tadd a spacer of the format <int>x<int>,<int>x<int>,"
                         "..., where the first integer corresponds to the space "
                         "between bases repeated the second integer number of times\n"
                         "-w\tSet window size [max(size of spaced kmer, [parameter])]\n"
                         "-S\tSet sketch size [16, for 2**16 bytes each]\n"
                         "-o\tOutput for genome size estimates [stdout]\n"
                         "-O\tOutput for genome distance matrix [stdout]\n"
                         "-e\tEmit in scientific notation\n"
                         "-F\tGet paths to genomes from file rather than positional arguments\n"
                         "TODO: Add separate sketch_usage.\n"
                , arg);
    std::exit(EXIT_FAILURE);
}

std::string &extend_suffix(std::string &suffix, int wsz, const Spacer &sp, int k, const std::string &spacing) {
    return suffix += ".w" + std::to_string(std::max((int)sp.c_, wsz)) + "." + std::to_string(k) + ".spacing" + spacing;
}

ks::string hll_fname(const char *path, size_t sketch_p, const std::string &suffix) {
    return ks::sprintf("%s.%s%zu.hll", get_cstr(path), (suffix.size() ? std::string(".")  + suffix + ".": std::string(".")).data(), sketch_p);
}

bool has_hll(const char *path, size_t sketch_p, const std::string &suffix) {
    return isfile(hll_fname(path, sketch_p, suffix).data());
}

// Main functions
int sketch_main(int argc, char *argv[]) {
    int wsz(-1), k(31), sketch_size(16), skip_cached(false), co;
    std::string spacing, paths_file, suffix;
    omp_set_num_threads(1);
    while((co = getopt(argc, argv, "F:c:p:x:s:S:k:w:ceh?")) >= 0) {
        switch(co) {
            case 'k': k = std::atoi(optarg); break;
            case 'x': suffix = optarg; break;
            case 'p': omp_set_num_threads(std::atoi(optarg)); break;
            case 's': spacing = optarg; break;
            case 'S': sketch_size = std::atoi(optarg); break;
            case 'w': wsz = std::atoi(optarg); break;
            case 'c': skip_cached = true; break;
            case 'F': paths_file = optarg; break;
            case 'h': case '?': dist_usage(*argv);
        }
    }
    spvec_t sv(spacing.size() ? parse_spacing(spacing.data(), k): spvec_t(k - 1, 0));
    Spacer sp(k, wsz, sv);
    extend_suffix(suffix, wsz, sp, k, spacing);
    std::vector<std::vector<std::string>> ivecs;
    {
        std::vector<std::string> inpaths(paths_file.size() ? get_paths(paths_file.data())
                                                           : std::vector<std::string>(argv + optind, argv + argc));
        for(const auto &el: inpaths) ivecs.emplace_back(std::vector<std::string>{el});
    }
    std::vector<hll::hll_t> hlls;
    while(hlls.size() < (unsigned)omp_get_num_threads()) hlls.emplace_back(sketch_size);
    if(wsz < sp.c_) wsz = sp.c_;
    if(ivecs.size() == 0) {
        std::fprintf(stderr, "No paths. See usage.\n");
        dist_usage(*argv);
    }
    #pragma omp parallel for
    for(size_t i = 0; i < hlls.size(); ++i) {
        const std::string &path(ivecs[i][0]);
        if(skip_cached && has_hll(path.data(), sketch_size, suffix)) continue;
        const auto tnum(omp_get_thread_num());
        hlls[tnum] = make_hll(ivecs[i],
                              k, wsz, sv, nullptr, 1, sketch_size);
        hlls[tnum].write(hll_fname(path.data(), sketch_size, suffix).data());
    }
    return EXIT_SUCCESS;
}

int dist_main(int argc, char *argv[]) {
    int wsz(-1), k(31), sketch_size(16), use_scientific(false), co, cache_sketch(false);
    std::string spacing, paths_file, suffix;
    FILE *ofp(stdout), *pairofp(stdout);
    omp_set_num_threads(1);
    while((co = getopt(argc, argv, "x:F:c:p:o:s:w:O:S:k:Meh?")) >= 0) {
        switch(co) {
            case 'k': k = std::atoi(optarg); break;
            case 'x': suffix = optarg; break;
            case 'p': omp_set_num_threads(std::atoi(optarg)); break;
            case 's': spacing = optarg; break;
            case 'S': sketch_size = std::atoi(optarg); break;
            case 'w': wsz = std::atoi(optarg); break;
            case 'F': paths_file = optarg; break;
            case 'o': ofp = fopen(optarg, "w"); if(ofp == nullptr) LOG_EXIT("Could not open file at %s for writing.\n", optarg); break;
            case 'O': pairofp = fopen(optarg, "w"); if(pairofp == nullptr) LOG_EXIT("Could not open file at %s for writing.\n", optarg); break;
            case 'e': use_scientific = true; break;
            case 'W': cache_sketch = true;  break;
            case 'h': case '?': dist_usage(*argv);
        }
    }
    spvec_t sv(spacing.size() ? parse_spacing(spacing.data(), k): spvec_t(k - 1, 0));
    Spacer sp(k, wsz, sv);
    extend_suffix(suffix, wsz, sp, k, spacing);
    std::vector<std::string> inpaths(paths_file.size() ? get_paths(paths_file.data())
                                                       : std::vector<std::string>(argv + optind, argv + argc));
    std::vector<hll::hll_t> hlls;
    std::vector<std::vector<std::string>> scratch_vv;
    while(scratch_vv.size() < (unsigned)omp_get_num_threads()) scratch_vv.emplace_back(std::vector<std::string>{"empty"}), scratch_vv.back()[0].reserve(256);
    while(hlls.size() < inpaths.size()) hlls.emplace_back(sketch_size);
    if(wsz < sp.c_) wsz = sp.c_;
    if(inpaths.size() == 0) {
        std::fprintf(stderr, "No paths. See usage.\n");
        dist_usage(*argv);
    }
    #pragma omp parallel for
    for(size_t i = 0; i < hlls.size(); ++i) {
        const std::string &path(inpaths[i]);
        if(cache_sketch && has_hll(path.data(), sketch_size, suffix))
            hlls[i].read(path);
        else {
            // By reserving 256 character, we make it probably that no allocation is necessary in this loop.
            std::vector<std::string> &scratch_stringvec(scratch_vv[omp_get_thread_num()]);
            scratch_stringvec[0] = inpaths[i];
#if 0
            hlls[i] = make_hll(scratch_stringvec,
                               k, wsz, sv, nullptr, 1, sketch_size);
#else
            fill_hll(hlls[i], scratch_stringvec, k, wsz, sv, nullptr, 1, sketch_size); // Avoid allocation fights.
#endif
        }
        if(cache_sketch) hlls[i].write(hll_fname(path.data(), sketch_size, suffix).data());
    }
    ks::string str("#Path\tSize (est.)\n");
    assert(str == "#Path\tSize (est.)\n");
    str.resize(1 << 18);
    {
        const int fn(fileno(ofp));
        for(size_t i(0); i < hlls.size(); ++i) {
            str.sprintf("%s\t%lf\n", inpaths[i].data(), hlls[i].report());
            if(str.size() >= 1 << 18) str.write(fn), str.clear();
        }
        str.write(fn), str.clear();
    }
    // TODO: Emit overlaps and symmetric differences.
    if(ofp != stdout) std::fclose(ofp);
    std::vector<double> dists(hlls.size() - 1);
    str.clear();
    str.sprintf("##Names \t");
    for(const auto &path: inpaths) str.sprintf("%s\t", path.data());
    str.back() = '\n';
    str.write(fileno(pairofp)); str.free();
    for(auto &el: hlls) el.sum();
    const char *fmt(use_scientific ? "\t%e": "\t%f");
    std::vector<hll::hll_t> scratch_hlls;
    while(scratch_hlls.size() < (unsigned)omp_get_num_threads()) scratch_hlls.emplace_back(sketch_size);
    for(size_t i = 0; i < hlls.size(); ++i) {
        hll::hll_t &h1(hlls[i]);
        #pragma omp parallel for
        for(size_t j = i + 1; j < hlls.size(); ++j)
            dists[j - i - 1] = jaccard_index(hlls[j], h1, scratch_hlls[omp_get_thread_num()]);
        h1.free();
        str += inpaths[i];
        for(size_t k(0); k < i + 1; ++k) str.putc_('\t'), str.putc_('-');
        for(size_t k(0); k < hlls.size() - i - 1; ++k) str.sprintf(fmt, dists[k]);
        str.putc_('\n');
        if(str.size() >= 1 << 18) str.write(fileno(pairofp)), str.clear();
    }
    str.write(fileno(pairofp)); str.clear();
    if(pairofp != stdout) fclose(pairofp);
    return EXIT_SUCCESS;
}

int setdist_main(int argc, char *argv[]) {
    int wsz(-1), k(31), use_scientific(false), co;
    unsigned bufsize(1 << 18);
    std::string spacing, paths_file;
    FILE *ofp(stdout), *pairofp(stdout);
    omp_set_num_threads(1);
    while((co = getopt(argc, argv, "F:c:p:o:O:S:B:k:Meh?")) >= 0) {
        switch(co) {
            case 'B': std::stringstream(optarg) << bufsize; break;
            case 'k': k = std::atoi(optarg); break;
            case 'p': omp_set_num_threads(std::atoi(optarg)); break;
            case 's': spacing = optarg; break;
            case 'w': wsz = std::atoi(optarg); break;
            case 'F': paths_file = optarg; break;
            case 'o': ofp = fopen(optarg, "w"); break;
            case 'O': pairofp = fopen(optarg, "w"); break;
            case 'e': use_scientific = true; break;
            case 'h': case '?': dist_usage(*argv);
        }
    }
    std::vector<char> rdbuf(bufsize);
    spvec_t sv(spacing.size() ? parse_spacing(spacing.data(), k): spvec_t(k - 1, 0));
    Spacer sp(k, wsz, sv);
    std::vector<std::string> inpaths(paths_file.size() ? get_paths(paths_file.data())
                                                       : std::vector<std::string>(argv + optind, argv + argc));
    std::vector<khash_t(all) *> hashes;
    while(hashes.size() < inpaths.size()) hashes.emplace_back((khash_t(all) *)calloc(sizeof(khash_t(all)), 1));
    const size_t nhashes(hashes.size());
    if(wsz < sp.c_) wsz = sp.c_;
    if(inpaths.size() == 0) {
        std::fprintf(stderr, "No paths. See usage.\n");
        dist_usage(*argv);
    }
    #pragma omp parallel for
    for(size_t i = 0; i < hashes.size(); ++i) {
        const char *path(inpaths[i].data());
        khash_t(all) *hash(hashes[i]);
        fill_set_genome<score::Lex>(path, sp, hash, i, nullptr);
    }
    LOG_DEBUG("Filled genomes. Now analyzing data.\n");
    ks::string str;
    str.sprintf("#Path\tSize (est.)\n");
    {
        const int fn(fileno(ofp));
        for(size_t i(0); i < hashes.size(); ++i) {
            str.sprintf("%s\t%zu\n", inpaths[i].data(), kh_size(hashes[i]));
            if(str.size() > 1 << 17) str.write(fn), str.clear();
        }
        str.write(fn), str.clear();
    }
    // TODO: Emit overlaps and symmetric differences.
    if(ofp != stdout) std::fclose(ofp);
    std::vector<double> dists(nhashes - 1);
    str.clear();
    str.sprintf("##Names \t");
    for(auto &path: inpaths) str.sprintf("%s\t", path.data());
    str.back() = '\n';
    str.write(fileno(pairofp)); str.free();
    setvbuf(pairofp, rdbuf.data(), _IOLBF, rdbuf.size());
    const char *const fmt(use_scientific ? "\t%e": "\t%f");
    for(size_t i = 0; i < hashes.size(); ++i) {
        auto &h1(hashes[i]);
        size_t j;
        #pragma omp parallel for
        for(j = i + 1; j < hashes.size(); ++j)
            dists[j - i - 1] = jaccard_index(hashes[j], h1);
        for(j = 0; j < i + 1; ++j) fputc('\t', pairofp), fputc('-', pairofp);
        for(j = 0; j < hashes.size() - i - 1; ++j)
            fprintf(pairofp, fmt, dists[j]);
        fputc('\n', pairofp);
        khash_destroy(h1), h1 = nullptr;
        // Delete data as soon as we don't need it.
    }
    return EXIT_SUCCESS;
}


}
