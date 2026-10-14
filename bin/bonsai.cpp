#include <fstream>
#include <sstream>
#include <omp.h>
#include "bonsai/feature_min.h"
#include "bonsai/util.h"
#include "bonsai/database.h"
#include "bonsai/classifier.h"
#include "bonsai/bitmap.h"
#include "bonsai/tx.h"
#include "bonsai/setcmp.h"
#include "bonsai/flextree.h"
#include "bonsai/grafting.h"
#include "bonsai/kraken_export.h"
#include "sketch/hll.h"
// #include "sketch/filterhll.h"

#if USE_CPPITERTOOLS
#include "cppitertools/groupby.hpp"
#endif

using namespace bns;

namespace {
inline double now_s() {
    return std::chrono::duration<double>(std::chrono::steady_clock::now().time_since_epoch()).count();
}
}

using namespace std::literals;
using std::cerr;
using std::cout;
using std::begin;
using std::end;
template<typename SketchType>
inline hll::hll_t &get_hll(SketchType &s);
template<>
inline hll::hll_t &get_hll<hll::hll_t>(hll::hll_t &s) {return s;}
/*
template<>
inline hll::hll_t &get_hll<fhll::fhll_t>(fhll::fhll_t &s) {return s.hll();}
template<>
inline hll::hll_t &get_hll<fhll::pcfhll_t>(fhll::pcfhll_t &s) {return s.hll();}
*/

template<typename SketchType> inline const hll::hll_t &get_hll(const SketchType &s);
template<>
inline const hll::hll_t &get_hll<hll::hll_t>(const hll::hll_t &s) {return s;}
/*
template<>
inline const hll::hll_t &get_hll<fhll::fhll_t>(const fhll::fhll_t &s) {return s.hll();}
*/

template<typename SketchType>
struct est_helper {
    const Spacer                      &sp_;
    const std::vector<std::string> &paths_;
    std::mutex                         &m_;
    const u64                          np_;
    const bool                      canon_;
    void                            *data_;
    std::vector<SketchType>         &hlls_;
    kseq_t                            *ks_;
};
template<typename SketchType, typename ScoreType=score::Lex>
void est_helper_fn(void *data_, long index, int tid) {
    est_helper<SketchType> &h(*(est_helper<SketchType> *)(data_));
    fill_lmers<ScoreType, SketchType>(h.hlls_[tid], h.paths_[index], h.sp_, h.canon_, h.data_, h.ks_ + tid);
}

template<typename SketchType, typename ScoreType=score::Lex>
void fill_sketch(SketchType &ret, const std::vector<std::string> &paths,
              unsigned k, uint16_t w, const spvec_t &spaces, bool canon=true,
              void *data=nullptr, int num_threads=1, u64 np=23, kseq_t *ks=nullptr) {
    // Default to using all available threads if num_threads is negative.
    if(num_threads < 0) {
        num_threads = std::thread::hardware_concurrency();
        LOG_INFO("Number of threads was negative and has been adjusted to all available threads (%i).\n", num_threads);
    }
    const Spacer space(k, w, spaces);
    if(num_threads <= 1) {
        for(u64 i(0); i < paths.size(); fill_lmers<ScoreType, SketchType>(ret, paths[i++], space, canon, data, ks));
    } else {
        std::mutex m;
        KSeqBufferHolder kseqs(num_threads);
        std::vector<SketchType> sketches;
        while(sketches.size() < (unsigned)num_threads) sketches.emplace_back(ret.clone());
        if(ret.size() != sketches.back().size()) throw "a party";
        est_helper<SketchType> helper{space, paths, m, np, canon, data, sketches, kseqs.data()};
        {
            ForPool pool(num_threads);
            pool.forpool(&est_helper_fn<SketchType, ScoreType>, &helper, paths.size());
        }
        auto &rhll = get_hll(ret);
        for(auto &sketch: sketches) rhll += get_hll(sketch);
    }
}

template<typename ScoreType=score::Lex>
hll::hll_t make_hll(const std::vector<std::string> &paths,
                unsigned k, uint16_t w, spvec_t spaces, bool canon=true,
                void *data=nullptr, int num_threads=1, u64 np=23, kseq_t *ks=nullptr, hll::EstimationMethod estim=hll::EstimationMethod::ERTL_MLE, uint16_t jestim=hll::JointEstimationMethod::ERTL_JOINT_MLE) {
    hll::hll_t global(np, estim, (hll::JointEstimationMethod)jestim);
    fill_sketch<hll::hll_t, ScoreType>(global, paths, k, w, spaces, canon, data, num_threads, np, ks);
    return global;
}

template<typename ScoreType=score::Lex>
u64 estimate_cardinality(const std::vector<std::string> &paths,
                            unsigned k, uint16_t w, spvec_t spaces, bool canon,
                            void *data=nullptr, int num_threads=-1, u64 np=23, kseq_t *ks=nullptr, hll::EstimationMethod estim=hll::EstimationMethod::ERTL_MLE) {
    auto tmp(make_hll<ScoreType>(paths, k, w, spaces, canon, data, num_threads, np, ks, estim));
    return tmp.report();
}

int classify_main(int argc, char *argv[]) {
    int co, num_threads(1), emit_kraken(1), emit_fastq(0), emit_all(0), chunk_size(1 << 20), per_set(32);
    bool canonicalize(true);
    std::ios_base::sync_with_stdio(false);
    std::FILE *ofp(stdout);
    if(argc < 4) {
        usage:
        std::fprintf(stderr, "Usage:\n%s <dbpath> <tax_path> <inr1.fq> [Optional: <inr2.fq>]\n"
                             "Flags:\n-o:\tRedirect output to path instead of stdout.\n"
                             "-c:\tSet chunk size. Default: %i\n"
                             "-a:\tEmit all records, not just classified.\n"
                             "-p:\tSet number of threads. [1] (Set -1 to use all threads.)\n"
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
    if(db.k_ > 32) LOG_EXIT("Database k %u is too large for 64-bit encodings.\n", db.k_);
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

bool endswith(const std::string &path, const std::string &suf) {
    return std::equal(std::crbegin(suf), std::crend(suf), std::crbegin(path));
}

int phase2_main(int argc, char *argv[]) {
    int c, mode(score_scheme::LEX), wsz(-1), num_threads(1), k(31);
    bool canon(true);
    WRITE write_fmt = UNCOMPRESSED;
    std::size_t start_size(1<<16);
    std::string spacing, tax_path, seq2taxpath, paths_file;
    std::ios_base::sync_with_stdio(false);
    std::string dbpath;
    // TODO: update documentation for tax_path and seq2taxpath options.
    if(argc < 4) {
        usage:
        std::fprintf(stderr, "Usage: %s <flags> -T <tax.path> -M <seq2tax.path> <out.path> <paths>\n"
                     "(minimized modes instead take <phase1map.path> <out.path> <paths>)\nFlags:\n"
                     "-k: Set k.\n"
                     "-p: Number of threads [1] (set to -1 to use all threads)\n"
                     "-t: Build for taxonomic minimizing\n-f: Build for feature minimizing\n"
                     "-F: Load paths from file provided instead further arguments on the command-line.\n"
                     "-e: Use entropy maximization.\n"
                     "-f: Use feature count minimization.\n"
                     "-t: Use tax depth maximization.\n"
                     "-w: Set window size.\n"
                     "-T: Set tax_path.\n"
                     "-M: Set seq2taxpath.\n"
                     "-S: Set spacing.\n"
                     "-z: Write gzip-compressed.\n"
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
            case 'z': write_fmt = ZLIB; break;
        }
    }
    dbpath = argv[optind];
    if(num_threads < 0) num_threads = std::thread::hardware_concurrency();
    if(wsz < k) wsz = k;
    // The build pipeline encodes k-mers in 64-bit words; k > 32 would alias.
    if(k > 32) LOG_EXIT("k must be <= 32 for this alphabet (64-bit encodings).\n");
#ifdef ZWRAP_USE_ZSTD
    const std::string suf(".zst");
#else
    const std::string suf(".gz");
#endif
    if(endswith(dbpath, suf))     write_fmt = ZLIB;
    if(write_fmt && !endswith(dbpath, ".gz"))
        dbpath += suf, LOG_INFO("Writing gzipped, but without a .gz suffix. Adding it.\n");
    LOG_INFO("db output path: %s\n", dbpath.data());
    spvec_t sv(parse_spacing(spacing.data(), k));
    // The lex/entropy branch takes positional <out.path> <paths...> (see the
    // README); the minimized branch below takes <phase1map> <out.path> <paths...>.
    std::vector<std::string> inpaths(paths_file.size() ? get_paths(paths_file.data())
                                                       : std::vector<std::string>(argv + optind + 1, argv + argc));
    if(inpaths.empty()) LOG_EXIT("Need input files from command line or file. See usage.\n");
    LOG_DEBUG("Got paths\n");
    if(seq2taxpath.empty()) LOG_EXIT("seq2taxpath required for final database generation.");
    if(score_scheme::LEX == mode || score_scheme::ENTROPY) {
        LOG_INFO("Final map will be written to %s\n", dbpath.data());
        Spacer sp(k, wsz, sv);
        Database<khash_t(c)>  phase2_map(sp);
        // Force using hll so that we can use __sync_bool_compare_and_swap to parallelize.
        LOG_INFO("About to estimate cardinality\n");
        std::size_t hash_size(estimate_cardinality<score::Lex>(inpaths, k, k, sv, canon, nullptr, num_threads, 24));
#if !NDEBUG
        {
            uint64_t sum = 0;
            for(const auto &path: inpaths) {
                gzFile fp = gzopen(path.data(), "rb");
                kseq_t *ks = kseq_init(fp);
                while(kseq_read(ks) >= 0) {
                    sum += ks->seq.l;
                }
                kseq_destroy(ks);
                gzclose(fp);
            }
            assert(sum > hash_size || !std::fprintf(stderr, "sum: %" PRIu64". hash size: %zu\n", sum, hash_size));
        }
#endif
        LOG_INFO("Estimated cardinality: %zu\n", hash_size);
        if(tax_path.empty()) RUNTIME_ERROR("Tax path required. [See -T option.]");
        LOG_INFO("Parent map bulding from %s\n", tax_path.data());
        khash_t(p) *taxmap(build_parent_map(tax_path.data()));
        //LOG_INFO("I just feel like stopping this executable now for testing.\n");
        //goto fail;
        phase2_map.db_ = score_scheme::LEX == mode ? lca_map<score::Lex>(inpaths, taxmap, seq2taxpath.data(), sp, num_threads, canon, hash_size)
                                                   : lca_map<score::Entropy>(inpaths, taxmap, seq2taxpath.data(), sp, num_threads, canon, hash_size);
        phase2_map.write(dbpath.data(), write_fmt);
        //fail:
        kh_destroy(p, taxmap);
        return EXIT_SUCCESS;
    }
    LOG_INFO("Making minimized map\n");
    // Minimized layout: <phase1map> <out.path> <paths...>; the first two
    // positionals are the input map and the output, not genomes.
    if(paths_file.empty()) inpaths.assign(argv + optind + 2, argv + argc);
    Database<khash_t(64)> phase1_map{Database<khash_t(64)>(dbpath.data())};
    Database<khash_t(c)>  phase2_map{phase1_map};
    Spacer sp(k, wsz, phase1_map.s_);
    khash_t(p) *taxmap(tax_path.empty() ? nullptr: build_parent_map(tax_path.data()));
    phase2_map.db_ = minimized_map<score::Hash>(inpaths, phase1_map.db_, seq2taxpath.data(), taxmap, sp, num_threads, start_size, canon);
    std::string dbpath2 = argv[optind + 1];
    if(endswith(dbpath2, suf))     write_fmt = ZLIB;
    if(write_fmt && !endswith(dbpath2, ".gz"))
        dbpath2 += suf, LOG_INFO("Writing gzipped, but without a .gz suffix. Adding it.\n");
    // Write minimized map
    phase2_map.write(dbpath2.data(), write_fmt);
    if(taxmap) kh_destroy(p, taxmap);
    return EXIT_SUCCESS;
}


int phase1_main(int argc, char *argv[]) {
    int c, taxmap_preparsed(0), use_hll(0), mode(score_scheme::LEX), wsz(-1), k(31), num_threads(1), sketch_size(24);
    bool canon(true);
    std::ios_base::sync_with_stdio(false);
    std::string spacing;

    if(argc < 5) {
        usage:
        std::fprintf(stderr, "Usage: %s <flags> <seq2tax.path> <taxmap.path> <out.path> <paths>\nFlags:\n"
                     "-k: Set k.\n"
                     "-p: Number of threads. [1] (Set to -1 to use all threads.)\n-S: add a spacer of the format "
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
    if(num_threads < 0) num_threads = std::thread::hardware_concurrency();
    wsz = k;
    if(k > 32) LOG_EXIT("k must be <= 32 for this alphabet (64-bit encodings).\n");
    khash_t(p) *taxmap(taxmap_preparsed ? khash_load<khash_t(p)>(argv[optind + 1])
                                        : build_parent_map(argv[optind + 1]));
    spvec_t sv(parse_spacing(spacing.data(), k));
    Spacer sp(k, wsz, sv);
    std::vector<std::string> inpaths(argv + optind + 3, argv + argc);
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
    for(const char *const *p(argv); *p; ++p) if((std::strcmp(*p, "-h") && std::strcmp(*p, "--help")) == 0) goto usage;
    if(argc == 1) {
        usage:
        std::fprintf(stderr, "Produces a histogram of how many kmers have been assigned to all taxids\n"
                             "Usage: bonsai %s <database.db> [outfile (omit to emit to stdout)]\n", *argv);
        std::exit(EXIT_FAILURE);
    }
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
    SORT(std::begin(structs), std::end(structs), [] (const elcount &a, const elcount &b) {return a.second < b.second;});
    std::fputs("Name\tCount\n", ofp);
    for(const auto &i: structs) std::fprintf(ofp, "%u\t%u\n", i.first, i.second);
    if(ofp != stdout) std::fclose(ofp);
     return EXIT_SUCCESS;
 }

int err_main(int argc, char *argv[]) {
    std::fprintf(stderr, "[bonsai:%s] No valid subcommand provided. Options: prebuild/p1/phase, build/p2/phase2, classify, metatree\n", BONSAI_VERSION);
    return EXIT_FAILURE;
}

int zomg_main(int argc, char *argv[]) {
    // This is testing code which will eventually prepare the reformed taxonomy (without and with)
    // new nodes for downstream work and classificaiton.
    std::vector<std::string> paths {"1", "2", "3"};
    TaxonomyReformation tr("Foobar", paths, nullptr);
    return EXIT_FAILURE;
}

int metatree_usage(const char *arg) {
    std::fprintf(stderr, "Usage: %s <db.path> <taxmap> <nameidmap> <out_taxmap> <out_taxkey>\n"
                         "\n"
                         "-F: Parse file paths from file instead of further arguments at command-line.\n"
                         "-d: Do not perform inversion (assume it's already been done.)\n"
                         "-f: Store binary dumps in folder <arg>.\n"
                         "-L: Set an lca to restrict analysis to one portion of the subtree.\n"
                         "    Repeating this option multiple times accepts genomes which are descendents of any taxid provided by this option.\n"
                         "-p: nthreads [1] (set to -1 to use all threads.)\n"
                         "-n: nthreads [1]\n"
                         "-r: Reform the taxonomy so that every genome gets a node (TODO.md item 2), write\n"
                         "    <arg>.names.gz, <arg>.o2n.gz, and <arg>.n2o.gz, and emit placements against the reformed ids.\n"
                 , arg);
    std::exit(EXIT_FAILURE);
    return EXIT_FAILURE;
}

bool accepted_pass(const khash_t(p) *taxmap, const std::vector<tax_t> &accepted, tax_t id) {
    if(accepted.empty()) return true;
    for(const auto el: accepted) if(lca(taxmap, el, id) == el) return true;
    return false;
}

int metatree_main(int argc, char *argv[]) {
    if(argc < 5) metatree_usage(*argv);
    int c, num_threads(1), k(31), nelem(0);
    size_t heap_size = 1 << 15;
    std::vector<tax_t> accept_lcas;
    FILE *ofp(stdout);
    std::string paths_file, folder, spacing, reform_prefix;
    std::ios_base::sync_with_stdio(false);
    while((c = getopt(argc, argv, "L:p:w:k:s:f:F:r:n:o:h?")) >= 0) {
        switch(c) {
            case '?': case 'h':     return metatree_usage(*argv);
            case 'f': folder      = optarg;                       break;
            case 'k': k           = std::atoi(optarg);            break;
            case 'F': paths_file  = optarg;                       break;
            case 'o': ofp         = std::fopen(optarg, "w");      break;
            case 'p': num_threads = std::atoi(optarg);            break;
            case 'n': nelem       = std::strtoull(optarg, 0, 10); break;
            case 'L': accept_lcas.push_back(std::atoi(optarg));   break;
            case 'r': reform_prefix = optarg;                     break;
        }
    }
    if(num_threads <= 0) num_threads = std::thread::hardware_concurrency();
    if(k > 32) LOG_EXIT("k must be <= 32 for this alphabet (64-bit encodings).\n");
    Spacer sp(k, k, nullptr);
    omp_set_num_threads(num_threads);
    khash_t(name) *name_hash(build_name_hash(argv[optind + 2]));
    LOG_DEBUG("Parsed name hash.\n");
    LOG_DEBUG("Got inpaths. Now building parent map\n");
    khash_t(p) *taxmap(build_parent_map(argv[optind + 1]));
    LOG_DEBUG("Got taxmap. Now getting maxtax\n");
    tax_t max_tax(get_max_val(taxmap));
    LOG_DEBUG("Now getting tax depths\n");
    auto tax_depths(get_tax_depths(taxmap, argv[optind + 1]));
    LOG_DEBUG("Got tax depths\n");
    std::unordered_set<tax_t> used_taxes;
    std::vector<std::string> inpaths(paths_file.size() ? get_paths(paths_file.data())
                                                       : std::vector<std::string>(argv + optind + 5, argv + argc));
    std::unordered_set<std::string> save;
    for(size_t i(0); i < inpaths.size(); ++i) {
        const auto &path(inpaths[i]);
#if !NDEBUG
        if((i % 500) == 0) LOG_DEBUG("At index %zu/%zu, save size is %zu\n", i, inpaths.size(), save.size());
#endif
        tax_t id;
        if((id = get_taxid(path.data(), name_hash)) != UINT32_C(-1)) {
            if(accepted_pass(taxmap, accept_lcas, id)) {
                save.insert(path), used_taxes.insert(id);
            }
        }
    }
    {
        inpaths = std::vector<std::string>(save.begin(), save.end());
        std::unordered_set<std::string> tmp;
        std::swap(tmp, save);
    }
    std::vector<tax_t> raw_taxes(used_taxes.begin(), used_taxes.end());
    for(auto tax: raw_taxes) {
        tax_t id;
        while((id = get_parent(taxmap, tax)))
            used_taxes.insert(id), tax = id;
    }
    LOG_DEBUG("Finished filtering genomes\n");

    if(inpaths.size() == 0)
        throw std::runtime_error("No input paths. I need to process genomes to tell you about them.");
    LOG_INFO("Processing %zu genomes\n", inpaths.size());

// Core
    std::unordered_map<tax_t, std::forward_list<std::string>> tx2desc_map;
    std::vector<tax_t> taxes;
    std::unique_ptr<TaxonomyReformation> tr;
    if(!reform_prefix.empty()) {
        // Build a concise taxonomy (a node per genome) from the genomes,
        // their taxid annotations, and the NCBI taxonomy, then place k-mers
        // against the reformed ids.
        tr.reset(new TaxonomyReformation(argv[optind + 2], inpaths, taxmap));
        tr->write_name_map((reform_prefix + ".names.gz").data());
        tr->write_old_to_new((reform_prefix + ".o2n.gz").data());
        tr->write_new_to_old((reform_prefix + ".n2o.gz").data());
        taxmap = const_cast<khash_t(p) *>(tr->pmap());
        max_tax = get_max_val(taxmap);
        for(khiter_t ki(0); ki != kh_end(taxmap); ++ki)
            if(kh_exist(taxmap, ki) && kh_key(taxmap, ki) != 1) taxes.push_back(kh_key(taxmap, ki));
        for(const auto &pair: tr->newid_paths()) { // genome leaf -> path; walk ancestors
            for(tax_t t(pair.first), p(t); t; t = p) {
                auto it(tx2desc_map.find(t));
                if(it == tx2desc_map.end()) it = tx2desc_map.emplace(t, std::forward_list<std::string>{}).first;
                it->second.push_front(pair.second);
                const khint_t ki(kh_get(p, taxmap, t));
                p = ki == kh_end(taxmap) ? 0 : kh_val(taxmap, ki);
            }
        }
    } else {
        taxes = vector_set_filter(get_sorted_taxes(taxmap, argv[optind + 1]), used_taxes);
        tx2desc_map = tax2desc_genome_map(tax2genome_map(name_hash, inpaths), taxmap, taxes, tax_depths);
    }
#if !NDEBUG
    for(const auto tax: taxes) assert(kh_get(p, taxmap, tax) != kh_end(taxmap));
    ks::string ks;
    const int stderrfn(fileno(stderr));
    ks.resize(1 << 6);
    typename decltype(tx2desc_map)::iterator it;
    for(const auto tax: taxes) {
        ks.putsn("Tax ", 4);
        ks.putuw_(tax);
        ks.puts(" has for descendent genomes: ");
        if((it = tx2desc_map.find(tax)) == tx2desc_map.end()) ks.putsn_("N/A\n", 4);
        else {
            for(const auto el: it->second) {
                ks.putsn(el.data(), el.size()); ks.putc(',');
            }
            ks.back() = '\n';
        }
        if(ks.size() & (1 << 16)) ks.write(stderrfn), ks.clear();
    }
    ks.write(stderrfn);
    ks.clear();
#endif
    // (canonicalize, max_heapsz, to_add) -- positional drift here once bound
    // heap_size to canonicalize and nelem to max_heapsz_, crashing fill_heap.
    FMEmitter fme(taxmap, tx2desc_map, true, heap_size, nelem);
    // Group taxes by parent and emit per subtree; no cppitertools needed.
    pdqsort(taxes.begin(), taxes.end(), [taxmap](const tax_t a, const tax_t b) {
        return get_parent(taxmap, a) < get_parent(taxmap, b);
    });
    for(auto bit(taxes.begin()); bit != taxes.end();) {
        auto eit(bit);
        while(eit != taxes.end() && get_parent(taxmap, *eit) == get_parent(taxmap, *bit)) ++eit;
        fme.process_subtree(get_parent(taxmap, *bit), bit, eit, sp, num_threads, nullptr);
        bit = eit;
    }

    fme.run_collapse(max_tax, ofp);
    if(ofp != stdout) std::fclose(ofp);
    return EXIT_SUCCESS;
}


// First-token-of-header taxid lookup: reads the first fasta header and
// probes the name hash with its first whitespace-delimited token (a
// "|kraken:taxid|" suffix, if present, is stripped). Returns tax_t(-1) when
// the token is absent from the name hash.
// Create dir and any missing parents (mkdir -p semantics).
static void make_dirs(const std::string &dir) {
    std::string cur;
    cur.reserve(dir.size());
    for(const auto c: dir) {
        cur.push_back(c);
        if(c == '/' && cur.size() > 1) ::mkdir(cur.data(), 0755);
    }
    ::mkdir(dir.data(), 0755);
}

static tax_t export_header_taxid(const std::string &path, const khash_t(name) *name_hash, std::string &seqid_out) {
    const std::string seqid(first_header_token(path.data()));
    if(seqid.empty()) LOG_EXIT("Empty, unreadable, or headerless fasta file %s\n", path.data());
    seqid_out = seqid;
    const khint_t ki(kh_get(name, name_hash, seqid_out.data()));
    return ki == kh_end(name_hash) ? tax_t(-1): kh_val(name_hash, ki);
}

static std::unordered_map<tax_t, std::string> load_taxon_strings(const std::string &path) {
    std::unordered_map<tax_t, std::string> ret;
    std::FILE *fp(std::fopen(path.data(), "r"));
    if(fp == nullptr) LOG_EXIT("Could not open %s\n", path.data());
    char *line(nullptr);
    size_t cap(0);
    while(getline(&line, &cap, fp) >= 0) {
        char *tab(std::strchr(line, '\t'));
        if(tab == nullptr) continue;
        *tab = 0;
        char *end(tab + 1);
        while(*end && *end != '\n' && *end != '\r') ++end;
        *end = 0;
        ret[std::atoi(line)] = tab + 1;
    }
    std::free(line);
    std::fclose(fp);
    return ret;
}

int export_usage(const char *arg) {
    std::fprintf(stderr, "Usage: %s <flags> -T <taxmap> -M <names> -o <outdir> [genome.fa ...]\n"
                         "Emits nodes.dmp/names.dmp/seqid2taxid.map and kraken library directories\n"
                         "that kraken2-build and krakenuniq-build ingest directly (see\n"
                         "include/bonsai/kraken_export.h). With -B > 0, greedily grafts up to B\n"
                         "synthetic nodes (candidate lanes: observed k-mer classes + BBMinHash\n"
                         "tree-layered sketches) before exporting.\n"
                         "\n"
                         "-T: Taxonomy map (child->parent; NCBI nodes.dmp-style rows accepted).\n"
                         "-M: Name map: \"<seqid>\\t<taxid>\" per line (krakenuniq seqid2taxid.map format).\n"
                         "-F: Read genome paths from file instead of the command line.\n"
                         "-o: Output directory (required).\n"
                         "-k: k-mer size for the grafting fill. [%u]\n"
                         "-B: Graft budget (0 = export the tree as-is). [0]\n"
                         "-P: Synthetic node name prefix. [grafted_]\n"
                         "-N: Optional \"<taxid>\\t<name>\" file for taxonomy node names.\n"
                         "-R: Optional \"<taxid>\\t<rank>\" file for taxonomy node ranks.\n"
                         "\n"
                         "PlacementTree/KmerClasses use u64 genome masks: at most 64 genomes per export.\n",
                 arg, 25);
    std::exit(EXIT_FAILURE);
    return EXIT_FAILURE;
}

int export_main(int argc, char *argv[]) {
    int c, k(25), budget(0);
    std::string tax_path, name_path, outdir, paths_file, graft_prefix("grafted_"), names_file, ranks_file;
    std::ios_base::sync_with_stdio(false);
    while((c = getopt(argc, argv, "T:M:F:o:k:B:P:N:R:h?")) >= 0) {
        switch(c) {
            case 'h': case '?': return export_usage(*argv);
            case 'T': tax_path = optarg;       break;
            case 'M': name_path = optarg;      break;
            case 'F': paths_file = optarg;     break;
            case 'o': outdir = optarg;         break;
            case 'k': k = std::atoi(optarg);   break;
            case 'B': budget = std::atoi(optarg); break;
            case 'P': graft_prefix = optarg;   break;
            case 'N': names_file = optarg;     break;
            case 'R': ranks_file = optarg;     break;
        }
    }
    if(tax_path.empty() || name_path.empty() || outdir.empty()) return export_usage(*argv);
    if(k > 32) LOG_EXIT("k must be <= 32 for this alphabet (64-bit encodings).\n");
    std::vector<std::string> paths(paths_file.size() ? get_paths(paths_file.data())
                                                     : std::vector<std::string>(argv + optind, argv + argc));
    if(paths.empty()) return export_usage(*argv);
    khash_t(p) *taxmap(build_parent_map(tax_path.data()));
    khash_t(name) *name_hash(build_name_hash(name_path.data()));
    std::vector<std::string> seqids, kept;
    std::vector<tax_t> genome_taxa;
    for(const auto &path: paths) {
        std::string seqid;
        const tax_t tid(export_header_taxid(path, name_hash, seqid));
        if(tid == tax_t(-1)) {
            LOG_WARNING("Genome %s (seqid %s) is not in the name map; skipping.\n", path.data(), seqid.data());
            continue;
        }
        if(kh_get(p, taxmap, tid) == kh_end(taxmap)) {
            LOG_WARNING("Taxid %u for genome %s is not in the taxonomy; skipping.\n", tid, path.data());
            continue;
        }
        seqids.push_back(seqid);
        kept.push_back(path);
        genome_taxa.push_back(tid);
    }
    destroy_name_hash(name_hash);
    if(kept.empty()) LOG_EXIT("No genomes survived name/taxonomy filtering.\n");
    if(kept.size() > 64) LOG_EXIT("%zu genomes provided, but the placement tree is limited to 64 (u64 genome masks).\n", kept.size());
    paths.swap(kept);
    LOG_INFO("Exporting %zu genomes against %zu taxonomy nodes.\n", paths.size(), kh_size(taxmap));
    PlacementTree tree(taxmap, genome_taxa);
    if(budget > 0) {
        Spacer sp(k, k);
        KmerClasses kc;
        const double t0(now_s());
        kc.fill<score::Lex>(paths, sp);
        const double fill_t(now_s() - t0);
        const u64 delta0(tree.total_delta(kc));
        std::vector<u64> cands;
        for(const auto &pair: kc.count) if(pc64(pair.first) >= 2) cands.push_back(pair.first);
        const size_t nclasses(cands.size());
        const auto sks(genome_sketchers<score::Lex>(paths, sp));
        for(const auto &comp: sketch_cluster_candidates_tree(sks, tree, 0.008)) {
            u64 m(0);
            for(const auto g: comp) m |= u64(1) << g;
            cands.push_back(m);
        }
        LOG_INFO("Grafting: %zu classes, %zu sketch candidates, index fill %.3f s, initial delta %llu.\n",
                 nclasses, cands.size() - nclasses, fill_t, (unsigned long long)delta0);
        const GraftingResult gr(graft_greedy(tree, kc, cands, u32(budget)));
        const HammingDecomp hd(hamming_decomp(tree, kc));
        LOG_INFO("Grafted %zu nodes: delta %llu -> %llu (fn %llu, fp %llu).\n",
                 gr.added_masks.size(), (unsigned long long)delta0,
                 (unsigned long long)gr.deltas.back(), (unsigned long long)hd.fn, (unsigned long long)hd.fp);
    }
    std::unordered_map<tax_t, std::string> taxon_names, taxon_ranks;
    if(!names_file.empty()) taxon_names = load_taxon_strings(names_file);
    if(!ranks_file.empty()) taxon_ranks = load_taxon_strings(ranks_file);
    make_dirs(outdir);
    const KrakenExportResult result(KrakenExporter::write(tree, seqids, paths, outdir.data(), graft_prefix.data(),
                                                          taxon_names.empty() ? nullptr: &taxon_names,
                                                          taxon_ranks.empty() ? nullptr: &taxon_ranks));
    kh_destroy(p, taxmap);
    return result.seqid2taxid.empty() ? EXIT_FAILURE: EXIT_SUCCESS;
}

int main(int argc, char *argv[]) {
    const static std::unordered_map<std::string, bns::MainFnPtr> md {
        {"phase1",   phase1_main},
        {"prebuild", phase1_main},
        {"p1",       phase1_main},
        {"phase2",   phase2_main},
        {"build",    phase2_main},
        {"p2",       phase2_main},
        {"lca",      phase1_main},
        {"hist",     hist_main},
        {"metatree", metatree_main},
        {"export",   export_main},
        {"classify", classify_main}
    };
    if(std::find_if(argv, argv + argc, [&](char *s) {return std::strcmp("-v", s) == 0 || std::strcmp("--version", s) == 0;}) != argv + argc) {
        std::fprintf(stdout, "bonsai|%s\n", BONSAI_VERSION);
        return EXIT_SUCCESS;
    }
    return (argc > 1 && md.find(argv[1]) != md.end() ? md.find(argv[1])->second
                                                     : err_main)(argc - 1, argv + 1);
}
