#include "lib/feature_min.h"
#include "lib/util.h"
#include "lib/database.h"
#include "lib/classifier.h"
#include "lib/tx.h"

using namespace emp;

int classify_main(int argc, char *argv[]) {
    int co, num_threads(16), emit_kraken(1), emit_fastq(0), emit_all(0), chunk_size(1 << 20), per_set(32);
    FILE *ofp(stdout);
    if(argc < 4) {
        usage:
        fprintf(stderr, "Usage:\n%s <dbpath> <tax_path> <inr1.fq> [Optional: <inr2.fq>]\n"
                        "Flags:\n-o:\tRedirect output to path instead of stdout.\n"
                        "-c:\tSet chunk size. Default: %i\n"
                        "-a:\tEmit all records, not just classified.\n"
                        "-p:\tSet number of threads. Default: 16.\n"
                        "-k:\tEmit kraken-style output.\n"
                        "-K:\tDo not emit kraken-style output.\n"
                        "-f:\tEmit fastq-style output.\n"
                        "-K:\tDo not emit fastq-formatted output.\n"
                        "\nIf -f and -k are set, full kraken output will be contained in the fastq comment field."
                        "\n  Default: kraken-style only output.\n",
                 *argv, 1 << 14);
        exit(EXIT_FAILURE);
    }
    while((co = getopt(argc, argv, "c:p:o:S:afFkKh?")) >= 0) {
        switch(co) {
            case 'h': case '?': goto usage;
            case 'a': emit_all = 1; break;
            case 'c': chunk_size = atoi(optarg); break;
            case 'F': emit_fastq  = 0; break;
            case 'f': emit_fastq  = 1; break;
            case 'K': emit_kraken = 0; break;
            case 'k': emit_kraken = 1; break;
            case 'p': num_threads = atoi(optarg); break;
            case 'o': ofp = fopen(optarg, "w"); break;
            case 'S': per_set = atoi(optarg); break;
        }
    }
    LOG_ASSERT(ofp);
    switch(argc - optind) {
        default: goto usage;
        case 3:  LOG_DEBUG("Processing in single-end mode.\n"); break;
        case 4:  LOG_DEBUG("Processing in paired-end mode.\n"); break;
    }
    Database<khash_t(c)> db(argv[optind]);
    //reportDB<khash_t(c)>(&db, stderr);
    //for(auto &i: db._s) --i; // subtract by one since we'll re-subtract during construction.
    ClassifierGeneric<lex_score> c(db.db_, db.s_, db.k_, db.k_, num_threads,
                                   emit_all, emit_fastq, emit_kraken);
    khash_t(p) *taxmap(build_parent_map(argv[optind + 1]));
    // We can use optind + 3 for both single-end and paired-end mode since the argument at
    // index argc is null when argc - optind == 3.
    process_dataset(c, taxmap, argv[optind + 2], argv[optind + 3],
                    ofp, chunk_size, per_set);
    if(ofp != stdout) fclose(ofp);
    kh_destroy(p, taxmap);
    LOG_INFO("Successfully completed classify!\n");
    return EXIT_SUCCESS;
}

int phase2_main(int argc, char *argv[]) {
    int c, mode(score_scheme::LEX), wsz(-1), num_threads(-1);
    unsigned k(31);
    std::size_t start_size(1<<16);
    std::string spacing, tax_path, seq2taxpath;
    // TODO: update documentation for tax_path and seq2taxpath options.
    if(argc < 4) {
        usage:
        fprintf(stderr, "Usage: %s <flags> [tax_path if lex else <phase1map.path>] <out.path> <paths>\nFlags:\n"
                "-k: Set k.\n"
                "-p: Number of threads\n"
                "-t: Build for taxonomic minimizing\n-f: Build for feature minimizing\n"
                , *argv);
        exit(EXIT_FAILURE);
    }
    while((c = getopt(argc, argv, "w:M:S:p:k:T:tfHh?")) >= 0) {
        switch(c) {
            case 'h': case '?': goto usage;
            case 'k': k = atoi(optarg); break;
            case 'p': num_threads = atoi(optarg); break;
            case 's': start_size = strtoull(optarg, nullptr, 10); break;
            case 't': mode = score_scheme::TAX_DEPTH; break;
            case 'f': mode = score_scheme::FEATURE_COUNT; break;
            case 'w': wsz = atoi(optarg); break;
            case 'T': tax_path = optarg; break;
            case 'M': seq2taxpath = optarg; break;
        }
    }
    if(wsz < 0 || wsz < (int)k) LOG_EXIT("Window size must be set and >= k for phase2.\n");
    const int lex(score_scheme::LEX == mode);
    if(lex) {
        if(seq2taxpath.empty()) LOG_EXIT("seq2taxpath required for lexicographic mode for final database generation.");
        Spacer sp(k, wsz, spvec_t(k - 1, 0));
        Database<khash_t(c)>  phase2_map(sp);
        std::vector<std::string> inpaths(argv + optind + 2, argv + argc);
        khash_t(p) *taxmap(build_parent_map(argv[optind]));
        phase2_map.db_ = lca_map<lex_score>(inpaths, taxmap, seq2taxpath.data(), sp, num_threads);
        phase2_map.write(argv[optind + 1]);
        kh_destroy(p, taxmap);
        return EXIT_SUCCESS;
    }
    Database<khash_t(64)> phase1_map(Database<khash_t(64)>(argv[optind]));
    Spacer sp(k, wsz, phase1_map.s_);
    Database<khash_t(c)>  phase2_map(phase1_map);
    std::vector<std::string> inpaths(argv + optind + 2, argv + argc);
    khash_t(p) *taxmap(tax_path.empty() ? nullptr: build_parent_map(tax_path.data()));
    phase2_map.db_ = minimized_map<hash_score>(inpaths, phase1_map.db_, sp, num_threads, start_size, mode);
    // Write minimized map
    phase2_map.write(argv[optind + 1]);
    if(taxmap) kh_destroy(p, taxmap);
    return EXIT_SUCCESS;
}

int phase1_main(int argc, char *argv[]) {
    int c, taxmap_preparsed(0), use_hll(0), mode(score_scheme::LEX), wsz(-1);

    unsigned k(31);
    int num_threads(-1);
    std::string spacing;
    if(argc < 5) {
        usage:
        fprintf(stderr, "Usage: %s <flags> <seq2tax.path> <taxmap.path> <out.path> <paths>\nFlags:\n"
                "-k: Set k.\n"
                "-p: Number of threads.\n-S: add a spacer of the format "
                "<int>,<int>,<int>, (...), where each integer is the number of spaces"
                "between successive bases included in the seed. There must be precisely k - 1"
                "elements in this list. Use this option multiple times to specify multiple seeds.\n"
                "-s: add a spacer of the format <int>x<int>,<int>x<int>,"
                "..., where the first integer corresponds to the space "
                "between bases repeated the second integer number of times.\n"
                "-t: Build for taxonomic minimizing.\n-f: Build for feature minimizing.\n"
                "-H: Estimate rather than count kmers exactly before building map.\n"
                "-T: Path to taxonomy map to load, if you've preparsed it. Not really worth it, building from scratch is fast.\n"
                "-d: Write out in database format version 1.\n"
                , *argv);
        exit(EXIT_FAILURE);
    }
    if(std::string("lca") == argv[0])
        fprintf(stderr, "[W:%s] lca subcommand has been renamed phase1. "
                        "This has been deprecated and will be removed.\n", __func__);
    while((c = getopt(argc, argv, "s:S:p:k:tfTHh?")) >= 0) {
        switch(c) {
            case 'h': case '?': goto usage;
            case 'k': k = atoi(optarg); break;
            case 'p': num_threads = atoi(optarg); break;
            case 's': spacing = optarg; break;
            case 'S': spacing = optarg; break;
            case 'T': taxmap_preparsed = 1; break;
            case 'H': use_hll = 1; break;
            case 't': mode = score_scheme::TAX_DEPTH; break;
            case 'f': mode = score_scheme::FEATURE_COUNT; break;
            //case 'w': wsz = atoi(optarg); break;
        }
    }
    if(wsz < 0) wsz = k;
    LOG_DEBUG("Loading tax\n");
    khash_t(p) *taxmap(taxmap_preparsed ? khash_load<khash_t(p)>(argv[optind + 1])
                                        : build_parent_map(argv[optind + 1]));
    LOG_DEBUG("Loaded tax\n");
    spvec_t sv(parse_spacing(spacing.data(), k));
    Spacer sp(k, wsz, sv);
    std::vector<std::string> inpaths(argv + optind + 3, argv + argc);
    {
        std::string tmp(inpaths[0]);
        for(std::size_t i(1); i < inpaths.size(); ++i) tmp += ",", tmp += inpaths[i];
        //fprintf(stderr, "Trying to gather kmers from %s.\n", tmp.data());
    }
    std::size_t hash_size(use_hll ? estimate_cardinality<lex_score>(inpaths, k, k, sv, nullptr, num_threads, 24): 1 << 16);
    if(use_hll) fprintf(stderr, "Estimated number of elements: %zu\n", hash_size);

    if(mode == score_scheme::LEX)
        LOG_EXIT("No phase1 required for lexicographic. Use phase2 instead.\n");
    auto mapbuilder(mode == score_scheme::TAX_DEPTH ? taxdepth_map<lex_score>
                                                    : ftct_map<lex_score>);
 
    Database<khash_t(64)> db(sp, 1, mapbuilder(inpaths, taxmap, argv[optind], sp, num_threads, hash_size));
    for(auto &i: db.s_) {
        LOG_DEBUG("Decrementing value %i to %i\n", i, i - 1);
        --i;
    }
    db.write(argv[optind + 2]);

    kh_destroy(p, taxmap);
    return EXIT_SUCCESS;
}

static std::vector<std::pair<std::string, int (*)(int, char **)>> mains {
    {"phase1", phase1_main},
    {"p1",     phase1_main},
    {"phase2", phase2_main},
    {"p2",     phase2_main},
    {"lca", phase1_main},
    {"classify", classify_main}
};
int main(int argc, char *argv[]) {
    
    if(argc > 1) for(auto &i: mains) if(i.first == argv[1]) return i.second(argc - 1, argv + 1);
    fprintf(stderr, "No valid subcommand provided. Options: phase1, phase2, classify\n");
    return EXIT_FAILURE;
}
