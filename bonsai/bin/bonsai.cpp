#include <fstream>
#include <omp.h>
#include "feature_min.h"
#include "util.h"
#include "database.h"
#include "classifier.h"
#include "bitmap.h"
#include "tx.h"
#include "setcmp.h"
#include "flextree.h"
#include "cppitertools/groupby.hpp"
#include "distmains.h"
#include "metamain.h"
#include <sstream>

using namespace emp;

using namespace std::literals;
using std::cerr;
using std::cout;
using std::begin;
using std::end;

int classify_main(int argc, char *argv[]) {
    int co, num_threads(16), emit_kraken(1), emit_fastq(0), emit_all(0), chunk_size(1 << 20), per_set(32);
    bool canonicalize(true);
    std::ios_base::sync_with_stdio(false);
    std::FILE *ofp(stdout);
    if(argc < 4) {
        usage:
        std::fprintf(stderr, "Usage:\n%s <dbpath> <tax_path> <inr1.fq> [Optional: <inr2.fq>]\n"
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
        std::exit(EXIT_FAILURE);
    }
    while((co = getopt(argc, argv, "Cc:p:o:S:afFkKh?")) >= 0) {
        switch(co) {
            case 'h': case '?': goto usage;
            case 'C': canonicalize = false; break;
            case 'a': emit_all = 1; break;
            case 'c': chunk_size = std::atoi(optarg); break;
            case 'F': emit_fastq  = 0; break;
            case 'f': emit_fastq  = 1; break;
            case 'K': emit_kraken = 0; break;
            case 'k': emit_kraken = 1; break;
            case 'p': num_threads = std::atoi(optarg); break;
            case 'o': ofp = std::fopen(optarg, "w"); break;
            case 'S': per_set = std::atoi(optarg); break;
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
    ClassifierGeneric<score::Lex> c(db.db_, db.s_, db.k_, db.k_, num_threads,
                                   emit_all, emit_fastq, emit_kraken, canonicalize);
    khash_t(p) *taxmap(build_parent_map(argv[optind + 1]));
    // We can use optind + 3 for both single-end and paired-end mode since the argument at
    // index argc is null when argc - optind == 3.
    process_dataset(c, taxmap, argv[optind + 2], argv[optind + 3],
                    ofp, chunk_size, per_set);
    if(ofp != stdout) std::fclose(ofp);
    kh_destroy(p, taxmap);
    LOG_INFO("Successfully completed classify!\n");
    return EXIT_SUCCESS;
}

int phase2_main(int argc, char *argv[]) {
    int c, mode(score_scheme::LEX), wsz(-1), num_threads(-1), k(31);
    bool canon(true);
    std::size_t start_size(1<<16);
    std::string spacing, tax_path, seq2taxpath, paths_file;
    std::ios_base::sync_with_stdio(false);
    // TODO: update documentation for tax_path and seq2taxpath options.
    if(argc < 4) {
        usage:
        std::fprintf(stderr, "Usage: %s <flags> [tax_path if lex/ent else <phase1map.path>] <out.path> <paths>\nFlags:\n"
                     "-k: Set k.\n"
                     "-p: Number of threads\n"
                     "-t: Build for taxonomic minimizing\n-f: Build for feature minimizing\n"
                     "-F: Load paths from file provided instead further arguments on the command-line.\n"
                     "-e: Use entropy maximization.\n"
                     "-f: Use feature count minimization.\n"
                     "-t: Use tax depth maximization.\n"
                     "-w: Set window size.\n"
                     "-T: Set tax_path.\n"
                     "-M: Set seq2taxpath.\n"
                     "-S: Set spacing.\n"
                     , *argv);
        std::exit(EXIT_FAILURE);
    }
    while((c = getopt(argc, argv, "Cw:M:S:p:k:T:F:tefHh?")) >= 0) {
        switch(c) {
            case 'C': canon = false; break;
            case 'h': case '?': goto usage;
            case 'k': k = std::atoi(optarg); break;
            case 'p': num_threads = std::atoi(optarg); break;
            case 'S': spacing = optarg; break;
            case 's': start_size = strtoull(optarg, nullptr, 10); break;
            case 't': mode = score_scheme::TAX_DEPTH; break;
            case 'f': mode = score_scheme::FEATURE_COUNT; break;
            case 'w': wsz = std::atoi(optarg); break;
            case 'T': tax_path = optarg; break;
            case 'M': seq2taxpath = optarg; break;
            case 'F': paths_file = optarg; break;
            case 'e': mode = score_scheme::ENTROPY; break;
        }
    }
    if(wsz < 0 || wsz < k) LOG_EXIT("Window size must be set and >= k for phase2.\n");
    spvec_t sv(spacing.size() ? parse_spacing(spacing.data(), k): spvec_t(k - 1, 0));
    std::vector<std::string> inpaths(paths_file.size() ? get_paths(paths_file.data())
                                                       : std::vector<std::string>(argv + optind + 2, argv + argc));
    if(inpaths.empty()) LOG_EXIT("Need input files from command line or file. See usage.\n");
    LOG_DEBUG("Got paths\n");
    if(score_scheme::LEX == mode || score_scheme::ENTROPY) {
        if(seq2taxpath.empty()) LOG_EXIT("seq2taxpath required for lexicographic mode for final database generation.");
        Spacer sp(k, wsz, sv);
        Database<khash_t(c)>  phase2_map(sp);
        // Force using hll so that we can use __sync_bool_compare_and_swap to parallelize.
#pragma message("Add cached/high water-point kseq buffers.")
        std::size_t hash_size(estimate_cardinality<score::Lex>(inpaths, k, k, sp.s_, canon, nullptr, num_threads, 24));
        LOG_DEBUG("Estimated cardinality: %zu\n", hash_size);
        LOG_DEBUG("Parent map bulding from %s\n", argv[optind]);
        khash_t(p) *taxmap(build_parent_map(argv[optind]));
        phase2_map.db_ = score_scheme::LEX == mode ? lca_map<score::Lex>(inpaths, taxmap, seq2taxpath.data(), sp, num_threads, canon, hash_size)
                                                   : lca_map<score::Entropy>(inpaths, taxmap, seq2taxpath.data(), sp, num_threads, canon, hash_size);
        phase2_map.write(argv[optind + 1]);
        kh_destroy(p, taxmap);
        return EXIT_SUCCESS;
    }
    Database<khash_t(64)> phase1_map{Database<khash_t(64)>(argv[optind])};
    Spacer sp(k, wsz, phase1_map.s_);
    Database<khash_t(c)>  phase2_map{phase1_map};
    khash_t(p) *taxmap(tax_path.empty() ? nullptr: build_parent_map(tax_path.data()));
    phase2_map.db_ = minimized_map<score::Hash>(inpaths, phase1_map.db_, sp, num_threads, start_size, canon);
    // Write minimized map
    phase2_map.write(argv[optind + 1]);
    if(taxmap) kh_destroy(p, taxmap);
    return EXIT_SUCCESS;
}


int hll_main(int argc, char *argv[]) {
    int c, wsz(-1), k(31), num_threads(-1), sketch_size(24);
    bool canon(true);
    std::string spacing, paths_file;
    std::ios_base::sync_with_stdio(false);
    if(argc < 2) {
        usage: LOG_EXIT("Usage: %s <opts> <paths>\nFlags:\n"
                        "-k:\tkmer length (Default: 31. Max: 31)\n"
                        "-w:\twindow size (Default: -1)  Must be -1 (ignored) or >= kmer length.\n"
                        "-s:\tspacing (default: none). format: <value>x<times>,<value>x<times>,...\n"
                        "   \tOmitting x<times> indicates 1 occurrence of spacing <value>\n"
                        "-S:\tsketch size (default: 24). (Allocates 2 << [param] bytes of memory per HyperLogLog.\n"
                        "-p:\tnumber of threads.\n"
                        "-F:\tPath to file which contains one path per line\n"
                        , argv[0]);
    }
    while((c = getopt(argc, argv, "Cw:s:S:p:k:tfh?")) >= 0) {
        switch(c) {
            case 'C': canon = false; break;
            case 'h': case '?': goto usage;
            case 'k': k = std::atoi(optarg); break;
            case 'p': num_threads = std::atoi(optarg); break;
            case 's': spacing = optarg; break;
            case 'S': sketch_size = std::atoi(optarg); break;
            case 'w': wsz = std::atoi(optarg); break;
            case 'F': paths_file = optarg; break;
        }
    }
    if(wsz < k) wsz = k;
    std::vector<std::string> inpaths(paths_file.empty() ? get_paths(paths_file.data())
                                                        : std::vector<std::string>(argv + optind, argv + argc));
    spvec_t sv(spacing.empty() ? spvec_t(k - 1, 0): parse_spacing(spacing.data(), k));
#pragma message("Add cached/high water-point kseq buffers.")
    const double est(estimate_cardinality<score::Lex>(inpaths, k, wsz, sv, canon, nullptr, num_threads, sketch_size));
    std::fprintf(stderr, "Estimated number of unique exact matches: %lf\n", est);
    return EXIT_SUCCESS;
}

int phase1_main(int argc, char *argv[]) {
    int c, taxmap_preparsed(0), use_hll(0), mode(score_scheme::LEX), wsz(-1), k(31), num_threads(-1), sketch_size(24);
    bool canon(true);
    std::ios_base::sync_with_stdio(false);
    std::string spacing;

    if(argc < 5) {
        usage:
        std::fprintf(stderr, "Usage: %s <flags> <seq2tax.path> <taxmap.path> <out.path> <paths>\nFlags:\n"
                     "-k: Set k.\n"
                     "-p: Number of threads.\n-S: add a spacer of the format "
                     "<int>,<int>,<int>, (...), where each integer is the number of spaces"
                     "between successive bases included in the seed. There must be precisely k - 1"
                     "elements in this list. Use this option multiple times to specify multiple seeds.\n"
                     "-s: add a spacer of the format <int>x<int>,<int>x<int>,"
                     "..., where the first integer corresponds to the space "
                     "between bases repeated the second integer number of times.\n"
                     "-S: Set HyperLogLog sketch size. For very large cardinalities, this may need to be increased for accuracy.\n"
                     "-t: Build for taxonomic minimizing.\n-f: Build for feature minimizing.\n"
                     "-H: Estimate rather than count kmers exactly before building map.\n"
                     "-T: Path to taxonomy map to load, if you've preparsed it. Not really worth it, building from scratch is fast.\n"
                     "-d: Write out in database format version 1.\n"
                     , *argv);
        std::exit(EXIT_FAILURE);
    }
    if("lca"s == argv[0])
        std::fprintf(stderr, "[W:%s] lca subcommand has been renamed phase1. "
                             "This has been deprecated and will be removed.\n", __func__);
    while((c = getopt(argc, argv, "Cs:S:p:k:tfTHh?")) >= 0) {
        switch(c) {
            case 'C': canon = false; break;
            case 'h': case '?': goto usage;
            case 'k': k = std::atoi(optarg); break;
            case 'p': num_threads = std::atoi(optarg); break;
            case 's': spacing = optarg; break;
            case 'S': sketch_size = std::atoi(optarg); break;
            case 'T': taxmap_preparsed = 1; break;
            case 'H': use_hll = 1; break;
            case 't': mode = score_scheme::TAX_DEPTH; break;
            case 'f': mode = score_scheme::FEATURE_COUNT; break;
            //case 'w': wsz = std::atoi(optarg); break;
        }
    }
    if(wsz < 0) wsz = k;
    khash_t(p) *taxmap(taxmap_preparsed ? khash_load<khash_t(p)>(argv[optind + 1])
                                        : build_parent_map(argv[optind + 1]));
    spvec_t sv(parse_spacing(spacing.data(), k));
    Spacer sp(k, wsz, sv);
    std::vector<std::string> inpaths(argv + optind + 3, argv + argc);
#pragma message("Add cached/high water-point kseq buffers.")
    std::size_t hash_size(use_hll ? estimate_cardinality<score::Lex>(inpaths, k, k, sv, canon, nullptr, num_threads, sketch_size): 1 << 16);
    if(use_hll) LOG_INFO("Estimated number of elements: %zu\n", hash_size);

    if(mode == score_scheme::LEX) LOG_EXIT("No phase1 required for lexicographic. Use phase2 instead.\n");
    auto mapbuilder(mode == score_scheme::TAX_DEPTH ? taxdepth_map<score::Lex>
                                                    : ftct_map<score::Lex>);
    Database<khash_t(64)> db(sp, 1, mapbuilder(inpaths, taxmap, argv[optind], sp, num_threads, canon, hash_size));
    for(auto &i: db.s_) {
        LOG_DEBUG("Decrementing value %i to %i\n", i, i - 1);
        --i;
    }
    db.write(argv[optind + 2]);

    kh_destroy(p, taxmap);
    return EXIT_SUCCESS;
}

int hist_main(int argc, char *argv[]) {
    Database<khash_t(c)> db(argv[1]);
    khash_t(c) *map(db.db_);
    std::FILE *ofp(stdout);
    count::Counter<u32> counter;
    if(argc > 2) ofp = std::fopen(argv[2], "w");
    for(khiter_t ki(0); ki != kh_end(map); ++ki) if(kh_exist(map, ki)) counter.add(kh_val(map, ki));
    auto &cmap(counter.get_map());
    using elcount = std::pair<tax_t, u32>;
    std::vector<elcount> structs;
    for(auto& i: cmap) structs.emplace_back(i.first, i.second);
    SORT(std::begin(structs), std::end(structs), [] (const elcount &a, const elcount &b) {
        return a.second < b.second;
    });
    std::fputs("Name\tCount\n", ofp);
    for(const auto &i: structs) std::fprintf(ofp, "%u\t%u\n", i.first, i.second);
    if(ofp != stdout) std::fclose(ofp);
     return EXIT_SUCCESS;
 }

int err_main(int argc, char *argv[]) {
    std::fputs("No valid subcommand provided. Options: phase1, phase2, classify, hll, metatree, dist, sketch, setdist\n", stderr);
    return EXIT_FAILURE;
}

int zomg_main(int argc, char *argv[]) {
    // This is testing code which will eventually prepare the reformed taxonomy (without and with)
    // new nodes for downstream work and classificaiton.
    std::vector<std::string> paths {"1", "2", "3"};
    TaxonomyReformation tr("Foobar", paths, nullptr);
    return EXIT_FAILURE;
}


int main(int argc, char *argv[]) {
    const static std::unordered_map<std::string, emp::MainFnPtr> md {
        {"phase1",   phase1_main},
        {"p1",       phase1_main},
        {"phase2",   phase2_main},
        {"p2",       phase2_main},
        {"lca",      phase1_main},
        {"hll",      hll_main},
        {"dist",     emp::dist_main},
        {"sketch",   emp::sketch_main},
        {"setdist",  emp::setdist_main},
        {"hist",     hist_main},
        {"metatree", metatree_main},
        {"classify", classify_main}
    };
    return (argc > 1 && md.find(argv[1]) != md.end() ? md.find(argv[1])->second
                                                     : err_main)(argc - 1, argv + 1);
}
