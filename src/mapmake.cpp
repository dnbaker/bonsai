#include "lib/feature_min.h"
#include "lib/util.h"

using namespace kpg;

int taxdbb_main(int argc, char *argv[]) {
    int c, taxmap_preparsed(0), use_hll(0), num_threads(-1), mode(score_scheme::LEX);
    unsigned k(31);
    std::string spacing;
    if(argc < 5) {
        usage:
        // TODO update this usage.
        fprintf(stderr, "Usage: %s <flags> <seq2tax.path> <taxmap.path> <out.path> <paths>\nFlags:\n"
                "-k: Set k.\n"
                "-p: Number of threads\n-S: add a spacer of the format"
                "<int>,<int>,<int>, (...), where each integer is the number of spaces"
                "between successive bases included in the seed. There must be precisely k - 1"
                "elements in this list. Use this option multiple times to specify multiple seeds.\n"
                , *argv);
        exit(EXIT_FAILURE);
    }
    while((c = getopt(argc, argv, "S:p:k:tfTHh?")) >= 0) {
        switch(c) {
            case 'h': case '?': goto usage;
            case 'k': k = atoi(optarg); break;
            case 'p': num_threads = atoi(optarg); break;
            case 'S': spacing = optarg; break;
            case 'T': taxmap_preparsed = 1; break;
            case 'H': use_hll = 1; break;
            case 't': mode = score_scheme::TAX_DEPTH; break;
            case 'f': mode = score_scheme::FEATURE_COUNT; break;
        }
    }
    spvec_t sv(encode_spacing(spacing.data()));
    Spacer sp(k, k, sv);
    std::vector<std::string> inpaths(argv + optind + 3, argv + argc);
    khash_t(64) *td(load_khash_map<khash_t(64)>(argv[optind]));
    // Then 
    kh_destroy(64, td);
    return EXIT_SUCCESS;
}

int lca_main(int argc, char *argv[]) {
    int c, taxmap_preparsed(0), use_hll(0);
    unsigned k(31);
    int num_threads(-1);
    std::string spacing;
    if(argc < 5) {
        usage:
        fprintf(stderr, "Usage: %s <flags> <seq2tax.path> <taxmap.path> <out.path> <paths>\nFlags:\n"
                "-k: Set k.\n"
                "-p: Number of threads\n-S: add a spacer of the format"
                "<int>,<int>,<int>, (...), where each integer is the number of spaces"
                "between successive bases included in the seed. There must be precisely k - 1"
                "elements in this list. Use this option multiple times to specify multiple seeds.\n"
                , *argv);
        exit(EXIT_FAILURE);
    }
    while((c = getopt(argc, argv, "S:p:k:tfTHh?")) >= 0) {
        switch(c) {
            case 'h': case '?': goto usage;
            case 'k': k = atoi(optarg); break;
            case 'p': num_threads = atoi(optarg); break;
            case 'S': spacing = optarg; break;
            case 'T': taxmap_preparsed = 1; break;
            case 'H': use_hll = 1; break;
            //case 't': mode = score_scheme::TAX_DEPTH; break;
            //case 'f': mode = score_scheme::FEATURE_COUNT; break;
        }
    }
    khash_t(p) *taxmap(taxmap_preparsed ? load_khash_map<khash_t(p)>(argv[optind + 1])
                                        : build_parent_map(argv[optind + 1]));
    spvec_t sv(encode_spacing(spacing.data()));
    Spacer sp(k, k, sv);
    std::vector<std::string> inpaths(argv + optind + 3, argv + argc);
#if !NDEBUG
    {
        std::string tmp(inpaths[0]);
        for(size_t i(1); i < inpaths.size(); ++i) tmp += ",", tmp += inpaths[i];
        //fprintf(stderr, "Trying to gather kmers from %s.\n", tmp.data());
    }
#endif
    size_t hash_size(use_hll ? estimate_cardinality<lex_score, 24>(inpaths, k, k, sv, nullptr, num_threads): 1 << 16);
    if(use_hll) fprintf(stderr, "Estimated number of elements: %" PRIu64 "\n", hash_size);
    khash_t(64) *out(taxdepth_map<lex_score>(inpaths, taxmap, argv[optind], sp, num_threads, hash_size));
    write_khash_map<khash_t(64)>(out, argv[optind + 2]);
    kh_destroy(p, taxmap);
    kh_destroy(64, out);
    return EXIT_SUCCESS;
}

static std::vector<std::pair<std::string, int (*)(int, char **)>> mains {
    {"lca", lca_main},
    //{"lca2depth", lca2depth_main},
    {"taxdbb", taxdbb_main}
};
int main(int argc, char *argv[]) {
    
    if(argc > 1) for(auto &i: mains) if(i.first == argv[1]) return i.second(argc - 1, argv + 1);
    fprintf(stderr, "No valid subcommand provided. Options: lca, lca2depth\n");
    return EXIT_FAILURE;
}
