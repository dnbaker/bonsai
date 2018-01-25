#include "metamain.h"
#include "util.h"
#include "flextree.h"
#include "bitmap.h"
#include "cppitertools/itertools.hpp"
#include <algorithm>
#include <omp.h>

namespace emp {

int metatree_usage(const char *arg) {
    std::fprintf(stderr, "Usage: %s <db.path> <taxmap> <nameidmap> <out_taxmap> <out_taxkey>\n"
                         "\n"
                         "-F: Parse file paths from file instead of further arguments at command-line.\n"
                         "-d: Do not perform inversion (assume it's already been done.)\n"
                         "-f: Store binary dumps in folder <arg>.\n"
                         "-L: Set an lca to restrict analysis to one portion of the subtree.\n"
                         "    Repeating this option multiple times accepts genomes which are descendents of any taxid provided by this option.\n"
                         "-p: nthreads [1]\n"
                         "-n: nthreads [1]\n"
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
    std::string paths_file, folder, spacing;
    std::ios_base::sync_with_stdio(false);
    while((c = getopt(argc, argv, "L:p:w:k:s:f:F:n:h?")) >= 0) {
        switch(c) {
            case '?': case 'h':     return metatree_usage(*argv);
            case 'f': folder      = optarg;                       break;
            case 'k': k           = std::atoi(optarg);            break;
            case 'F': paths_file  = optarg;                       break;
            case 'o': ofp         = std::fopen(optarg, "w");      break;
            case 'p': num_threads = std::atoi(optarg);            break;
            case 'n': nelem       = std::strtoull(optarg, 0, 10); break;
            case 'L': accept_lcas.push_back(std::atoi(optarg));   break;
        }
    }
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
    std::vector<tax_t> taxes(get_sorted_taxes(taxmap, argv[optind + 1]));
    taxes = vector_set_filter(taxes, used_taxes);
    std::cerr << "Got sorted taxes\n";
    auto tx2desc_map(tax2desc_genome_map(tax2genome_map(name_hash, inpaths), taxmap, taxes, tax_depths));
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
    // TODO: Add new taxonomy creation
    FMEmitter fme(taxmap, tx2desc_map, heap_size, nelem);
    std::vector<tax_t> tmptaxes;
    for(auto &&pair: iter::groupby(taxes, [tm=taxmap](const tax_t a){return get_parent(tm, a);})) {
        // Copying just because I don't trust the lifetime management of iter::groupby.
        for(const auto tax: pair.second) tmptaxes.push_back(tax);
        std::cerr << "Going through taxes with parent = " << static_cast<tax_t>(pair.first) << '\n';
        fme.process_subtree(pair.first, begin(tmptaxes), end(tmptaxes), sp, num_threads, nullptr);
        tmptaxes.clear();
    }
    fme.run_collapse(max_tax, ofp);
    if(ofp != stdout) std::fclose(ofp);
    return EXIT_SUCCESS;
}


}
