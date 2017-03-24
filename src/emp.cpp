#include "lib/feature_min.h"
#include "lib/util.h"
#include "lib/database.h"
#include "lib/classifier.h"
#include "lib/tree_climber.h"
#include "lib/bitmap.h"
#include "lib/tx.h"
#include "lib/khpp.h"
#include <functional>

using namespace emp;

using namespace std::literals;

int classify_main(int argc, char *argv[]) {
    int co, num_threads(16), emit_kraken(1), emit_fastq(0), emit_all(0), chunk_size(1 << 20), per_set(32);
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
            case 'o': ofp = std::fopen(optarg, "w"); break;
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
    if(ofp != stdout) std::fclose(ofp);
    kh_destroy(p, taxmap);
    LOG_INFO("Successfully completed classify!\n");
    return EXIT_SUCCESS;
}

std::vector<std::string> get_paths(const char *path) {
    gzFile fp(gzopen(path, "rb"));
    char buf[1024], *line;
    std::vector<std::string> ret;
    while((line = gzgets(fp, buf, sizeof buf))) {
        ret.emplace_back(line);
        ret[ret.size() - 1].pop_back();
    }
    gzclose(fp);
    return ret;
}

int phase2_main(int argc, char *argv[]) {
    int c, mode(score_scheme::LEX), wsz(-1), num_threads(-1), k(31);
    std::size_t start_size(1<<16);
    std::string spacing, tax_path, seq2taxpath, paths_file;
    // TODO: update documentation for tax_path and seq2taxpath options.
    if(argc < 4) {
        usage:
        std::fprintf(stderr, "Usage: %s <flags> [tax_path if lex else <phase1map.path>] <out.path> <paths>\nFlags:\n"
                     "-k: Set k.\n"
                     "-p: Number of threads\n"
                     "-t: Build for taxonomic minimizing\n-f: Build for feature minimizing\n"
                     "-F: Load paths from file provided instead further arguments on the command-line.\n"
                     , *argv);
        std::exit(EXIT_FAILURE);
    }
    while((c = getopt(argc, argv, "w:M:S:p:k:T:F:tfHh?")) >= 0) {
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
            case 'F': paths_file = optarg; break;
        }
    }
    if(wsz < 0 || wsz < k) LOG_EXIT("Window size must be set and >= k for phase2.\n");
    std::vector<std::string> inpaths(paths_file.size() ? get_paths(paths_file.data())
                                                       : std::vector<std::string>(argv + optind + 2, argv + argc));
    LOG_DEBUG("Got paths\n");
    if(score_scheme::LEX == mode) {
        if(seq2taxpath.empty()) LOG_EXIT("seq2taxpath required for lexicographic mode for final database generation.");
        Spacer sp(k, wsz, spvec_t(k - 1, 0));
        Database<khash_t(c)>  phase2_map(sp);
#if 0
        std::size_t hash_size(use_hll ? estimate_cardinality<lex_score>(inpaths, k, k, sp.s_, nullptr, num_threads, 24): 1 << 16);
#else
        // Force using hll so that we can use __sync_bool_compare_and_swap to parallelize.
        std::size_t hash_size(estimate_cardinality<lex_score>(inpaths, k, k, sp.s_, nullptr, num_threads, 24));
        LOG_DEBUG("Estimated cardinality: %zu\n", hash_size);
#endif
        LOG_DEBUG("Parent map bulding from %s\n", argv[optind]);
        khash_t(p) *taxmap(build_parent_map(argv[optind]));
        phase2_map.db_ = lca_map<lex_score>(inpaths, taxmap, seq2taxpath.data(), sp, num_threads, hash_size);
        phase2_map.write(argv[optind + 1]);
        kh_destroy(p, taxmap);
        return EXIT_SUCCESS;
    }
    Database<khash_t(64)> phase1_map(Database<khash_t(64)>(argv[optind]));
    Spacer sp(k, wsz, phase1_map.s_);
    Database<khash_t(c)>  phase2_map(phase1_map);
    khash_t(p) *taxmap(tax_path.empty() ? nullptr: build_parent_map(tax_path.data()));
    phase2_map.db_ = minimized_map<hash_score>(inpaths, phase1_map.db_, sp, num_threads, start_size);
    // Write minimized map
    phase2_map.write(argv[optind + 1]);
    if(taxmap) kh_destroy(p, taxmap);
    return EXIT_SUCCESS;
}


int hll_main(int argc, char *argv[]) {
    int c, wsz(-1), k(31), num_threads(-1), sketch_size(24);
    std::string spacing;
    if(argc < 2) {
        usage: LOG_EXIT("Usage: %s <opts> <paths>\nFlags:"
                        "-k:\tkmer length (Default: 31. Max: 31)\n"
                        "-w:\twindow size (Default: -1)  Must be -1 (ignored) or >= kmer length.\n"
                        "-s:\tspacing (default: none). format: <value>x<times>,<value>x<times>,...\n"
                        "   \tOmitting x<times> indicates 1 occurrence of spacing <value>\n"
                        "-S:\tsketch size (default: 24). (Allocates 2 << [param] bytes of memory per HyperLogLog.\n"
                        "-p:\tnumber of threads.\n"
                        , argv[0]);
    }
    while((c = getopt(argc, argv, "w:s:S:p:k:tfh?")) >= 0) {
        switch(c) {
            case 'h': case '?': goto usage;
            case 'k': k = atoi(optarg); break;
            case 'p': num_threads = atoi(optarg); break;
            case 's': spacing = optarg; break;
            case 'S': sketch_size = atoi(optarg); break;
            case 'w': wsz = atoi(optarg); break;
        }
    }
    if(wsz < k) wsz = k;
    std::vector<std::string> inpaths(argv + optind, argv + argc);
    spvec_t sv(spacing.empty() ? spvec_t(k - 1, 0): parse_spacing(spacing.data(), k));
    std::size_t est(estimate_cardinality<lex_score>(inpaths, k, wsz, sv, nullptr, num_threads, sketch_size));
    std::fprintf(stderr, "Estimated number of unique exact matches: %zu\n", est);
    return EXIT_SUCCESS;
}

int phase1_main(int argc, char *argv[]) {
    int c, taxmap_preparsed(0), use_hll(0), mode(score_scheme::LEX), wsz(-1), k(31), num_threads(-1), sketch_size(24);
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
    while((c = getopt(argc, argv, "s:S:p:k:tfTHh?")) >= 0) {
        switch(c) {
            case 'h': case '?': goto usage;
            case 'k': k = atoi(optarg); break;
            case 'p': num_threads = atoi(optarg); break;
            case 's': spacing = optarg; break;
            case 'S': sketch_size = atoi(optarg); break;
            case 'T': taxmap_preparsed = 1; break;
            case 'H': use_hll = 1; break;
            case 't': mode = score_scheme::TAX_DEPTH; break;
            case 'f': mode = score_scheme::FEATURE_COUNT; break;
            //case 'w': wsz = atoi(optarg); break;
        }
    }
    if(wsz < 0) wsz = k;
    khash_t(p) *taxmap(taxmap_preparsed ? khash_load<khash_t(p)>(argv[optind + 1])
                                        : build_parent_map(argv[optind + 1]));
    spvec_t sv(parse_spacing(spacing.data(), k));
    Spacer sp(k, wsz, sv);
    std::vector<std::string> inpaths(argv + optind + 3, argv + argc);
    std::size_t hash_size(use_hll ? estimate_cardinality<lex_score>(inpaths, k, k, sv, nullptr, num_threads, sketch_size): 1 << 16);
    if(use_hll) std::fprintf(stderr, "Estimated number of elements: %zu\n", hash_size);

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

void metatree_usage(char *arg) {
    std::fprintf(stderr, "Usage: %s <db.path> <taxmap> <nameidmap> <out_taxmap> <out_taxkey>\n"
                         "\n"
                         "-F: Parse file paths from file instead of further arguments at command-line.\n"
                         "-d: Do not perform inversion (assume it's already been done.)\n"
                         "-f: Store binary dumps in folder <arg>.\n"
                 , arg);
    std::exit(EXIT_FAILURE);
}

//typedef std::tuple<std::uint64_t, int, std::vector<std::uint64_t>> indvec_t;

template struct kh::khpp_t<std::vector<std::uint64_t> *, std::uint64_t, ptr_wang_hash_struct<std::vector<std::uint64_t> *>>;
using pkh_t = kh::khpp_t<std::vector<std::uint64_t> *, std::uint64_t, ptr_wang_hash_struct<std::vector<std::uint64_t> *>>;

struct potential_node_t {
    const tax_t                 p_;
    std::vector<std::uint64_t> *bits_;
    std::uint64_t               score_;
    bool operator<(potential_node_t &other) {
        return std::tie(score_, p_, bits_) < std::tie(other.score_, other.p_, other.bits_);
    }
    bool operator>(potential_node_t &other) {
        return std::tie(score_, p_, bits_) > std::tie(other.score_, other.p_, other.bits_);
    }
    bool operator==(potential_node_t &other) {
        return std::tie(score_, p_, bits_) == std::tie(other.score_, other.p_, other.bits_);
    }
};

struct tree_glob_t {
    using tax_path_map_t = std::unordered_map<std::uint32_t, std::forward_list<std::string>>;
    const tax_t                                    parent_;
    const std::string                              parent_path_;
    khash_t(all)                                  *acceptable_;
    std::set<tax_t>                                taxes_;
    count::Counter<std::vector<std::uint64_t>>     counts_;


    static constexpr const char *KMER_SUFFIX = ".kmers.bin";


    std::string make_parent(const std::string &fld, const tax_t parent) {
        return std::string(fld.empty() ? "": fld + '/') + std::to_string(parent) + KMER_SUFFIX;
    }
    tree_glob_t(khash_t(p) *tax, tax_t parent, const std::string &fld, const Spacer &sp, tax_path_map_t &tpm,
                const std::unordered_map<tax_t, std::vector<tax_t>> &invert, int num_threads=16):
        parent_(parent), parent_path_(make_parent(fld, parent_)), acceptable_(tree::load_binary_kmerset(parent_path_.data()))
    {
        get_taxes(invert);
        tax_path_map_t tmp;
        for(auto tax: taxes_) tmp.emplace(tax, tpm[tax]);
        counts_ = std::move(bitmap_t(kgset_t(tmp, sp, num_threads, acceptable_)).to_counter());
        khash_destroy(acceptable_);
        acceptable_ = nullptr;
    }
    void add(const tax_t tax) {
        //if(get_parent(tax_, tax) != parent_) LOG_EXIT("Unexpected node whose parent (%u) is not as expected (%u)\n", get_parent(tax_, tax), parent_);
        taxes_.insert(tax);
    }
    template<typename Container>
    void add_children(Container &&c) {
        for(auto tax: c)    add(tax); // Add all the taxes.
    }
    template<typename Container>
    void add_children(Container &c) {
        for(auto tax: c)    add(tax); // Add all the taxes.
    }
    template<typename It>
    void add_children(It first, It end) {
        while(first != end) add(*first);
    }
    void get_taxes(const std::unordered_map<tax_t, std::vector<tax_t>> &invert) {
        add_children(get_all_descendents(invert, parent_));
    }
    ~tree_glob_t() {if(acceptable_) khash_destroy(acceptable_);}
};

struct tree_adjudicator_t {
    std::vector<tree_glob_t>   subtrees_;
    std::vector<adjmap_t>      fwds_;
    std::vector<adjmap_t>      revs_;
    std::set<potential_node_t> nodes_;
    const std::uint64_t        original_tax_count_; // Needed for scoring
    const std::uint64_t        max_el_;             // Needed to know which nodes to add.
    int                        recalculate_scores_; // If yes, recalculate scores 
    tree_adjudicator_t(std::uint64_t orig, int recalculate=1):
        original_tax_count_(orig), max_el_(roundup64(orig)), recalculate_scores_(recalculate)
    {
    }
    template<typename... U>
    void process_subtree(U&&... Args) {
        subtrees_.emplace_back(std::forward<U>(Args)...);
        auto &last = *(subtrees_.end() - 1);
        fwds_.emplace_back(last.counts_, false);
        revs_.emplace_back(last.counts_, true);
        // Add potential_node_t's from each subtree using default scoring
    }
    void process() {
        for(auto &tree: subtrees_) process_subtree(tree);
    }
};



int metatree_main(int argc, char *argv[]) {
    if(argc < 5) metatree_usage(*argv);
    int c, dry_run(0), num_threads(-1);
    std::string paths_file, folder, spacing;
    while((c = getopt(argc, argv, "w:k:s:f:F:h?d")) >= 0) {
        switch(c) {
            case '?': case 'h': metatree_usage(*argv);
            case 'f': folder = optarg; break;
            case 'F': paths_file = optarg; break;
            case 'd': dry_run = 1; break;
            case 'p': num_threads = atoi(optarg); break;
        }
    }
    khash_t(name) *name_hash(build_name_hash(argv[optind + 2]));
    khash_t(p) *tmp_taxmap(build_parent_map(argv[optind + 1]));
    std::vector<std::string> inpaths(paths_file.size() ? get_paths(paths_file.data())
                                                       : std::vector<std::string>(argv + optind + 5, argv + argc));
    if(inpaths.empty()) LOG_EXIT("Need input files from command line or file. See usage.\n");
    khash_t(p) *taxmap(tree::pruned_taxmap(inpaths, tmp_taxmap, name_hash));
    std::unordered_map<tax_t, std::vector<tax_t>> inverted_map(invert_parent_map(taxmap));
    kh_destroy(p, tmp_taxmap);
    std::size_t num_nodes(kh_size(taxmap));
    tree::SortedNodeGuide guide(taxmap);
    Database<khash_t(c)> db(argv[optind]);
    Spacer sp(db.k_, db.w_, db.s_);
    for(auto &i: sp.s_) --i;
    std::vector<tax_t> nodes(std::move(guide.get_nodes()));
    for(auto tax: nodes) {
        std::fprintf(stderr, "node %u has parent %u and depth %u\n", tax, get_parent(taxmap, tax), node_depth(taxmap, tax));
    }
    std::vector<std::string> to_fetch;
#if 1
    if(!dry_run) to_fetch = std::move(tree::invert_lca_map(db, folder.data()));
    else {
        to_fetch.reserve(kh_size(taxmap));
        char buf[256];
        for(khiter_t ki(0); ki < kh_end(taxmap); ++ki) {
            if(!kh_exist(taxmap, ki)) continue;
            std::sprintf(buf, "%s%u.kmers.bin", folder.data(), kh_val(taxmap, ki));
            if(access(buf, F_OK) != -1) to_fetch.emplace_back(buf); // Only add to list if file exists.
        }
    }
#else
    to_fetch = std::move(tree::par_invert(db, folder.data()));
#endif
    if(to_fetch.empty()) LOG_EXIT("No binary files to grab from.\n");
    LOG_DEBUG("Fetched! Making tx2g\n");
    std::unordered_map<tax_t, strlist> tx2g(tax2genome_map(name_hash, inpaths));
    LOG_DEBUG("Made! Printing\n");
    for(auto &kv: tx2g)
        for(auto &path: kv.second)
            std::fprintf(stderr, "Taxid %u has path %s and first line %s\n", kv.first, path.data(), get_firstline(path.data()).data());
    std::size_t index(0);
    std::vector<tax_t> parents;
    std::vector<std::vector<tax_t>> descendents;
    for(auto tax_iter(std::begin(nodes)), end_iter(tax_iter);tax_iter + 1 < std::end(nodes);tax_iter = end_iter + 1) {
        const tax_t parent_tax(get_parent(taxmap, *tax_iter));
        end_iter = tax_iter;
        while(end_iter != std::end(nodes) && get_parent(taxmap, *++end_iter) == parent_tax);
        assert(get_parent(taxmap, *(end_iter - 1)) == parent_tax);
        assert(get_parent(taxmap, *end_iter) != parent_tax);
        std::set<tax_t> taxes(tax_iter, end_iter);
        std::unordered_map<tax_t, std::forward_list<std::string>> range_map;
        std::vector<std::forward_list<std::string>> list_vec;
        std::vector<tax_t> parent_descendents(get_all_descendents(inverted_map, parent_tax));
        for(auto tax: parent_descendents) {
            auto m(tx2g.find(tax));
            if(m == tx2g.end()) continue;
            range_map.emplace(tax, m->second);
        }
        for(auto &pair: range_map) list_vec.push_back(pair.second); // Copying strings. Still not bad because there are only thousands of strings total.
        auto n([](std::unordered_map<tax_t, std::forward_list<std::string>> &m){
            std::size_t ret(0);for(auto &pair: m) for(auto &str: pair.second) ++ret; return ret;
        }(range_map));
        std::fprintf(stderr, "%zu genomes who have %u as their parent\n", n, parent_tax);
        if(n < 3) continue;
        parents.push_back(parent_tax);
        std::string parent_path(folder);
        if(!parent_path.empty()) parent_path += '/';
        parent_path += std::to_string(parent_tax) + ".kmers.bin";
        LOG_INFO("Get acceptable\n");
        khash_t(all) *acceptable(tree::load_binary_kmerset(parent_path.data()));
        bitmap_t bitmap;
        {
            kgset_t kgs(range_map, sp, num_threads, acceptable);
            bitmap = kgs;
            descendents.emplace_back(std::move(kgs.get_taxes()));
        }
        count::Counter<std::vector<std::uint64_t>> counts(bitmap.to_counter());
        adjmap_t fwd_adj(counts);
        adjmap_t rev_adj(counts, true);
    }
    destroy_name_hash(name_hash);
    kh_destroy(p, taxmap);
    LOG_DEBUG("I ran to completion... How???\n");
    return EXIT_SUCCESS;
}

int hist_main(int argc, char *argv[]) {
    Database<khash_t(c)> db(argv[1]);
    khash_t(c) *map(db.db_);
    std::FILE *ofp(stdout);
    count::Counter<std::uint32_t> counter;
    if(argc > 2) ofp = std::fopen(argv[2], "w");
    for(khiter_t ki(0); ki != kh_end(map); ++ki) if(kh_exist(map, ki)) counter.add(kh_val(map, ki));
    auto &cmap(counter.get_map());
    using elcount = std::pair<tax_t, std::uint32_t>;
    std::vector<elcount> structs;
    for(auto& i: cmap) structs.emplace_back(i.first, i.second);
    std::sort(std::begin(structs), std::end(structs), [] (elcount &a, elcount &b) {
        return a.second < b.second;
    });
    std::fputs("Name\tCount\n", ofp);
    for(auto &i: structs) std::fprintf(ofp, "%u\t%u\n", i.first, i.second);
    if(ofp != stdout) std::fclose(ofp);
    return EXIT_SUCCESS;
}

static std::vector<std::pair<std::string, int (*)(int, char **)>> mains {
    {"phase1", phase1_main},
    {"p1",     phase1_main},
    {"phase2", phase2_main},
    {"p2",     phase2_main},
    {"lca", phase1_main},
    {"hll", hll_main},
    {"hist", hist_main},
    {"metatree", metatree_main},
    {"classify", classify_main}
};
int main(int argc, char *argv[]) {

    if(argc > 1) for(auto &i: mains) if(i.first == argv[1]) return i.second(argc - 1, argv + 1);
    std::fprintf(stderr, "No valid subcommand provided. Options: phase1, phase2, classify, hll, metatree\n");
    return EXIT_FAILURE;
}
