#include "util.h"
#include "getopt.h"

using namespace emp;

void print_vec(const std::vector<std::string> &vec) {
    for(const auto &el: vec) std::cerr << el << ", ";
    std::cerr << '\n';
}

int main(int argc, char *argv[]) {
    int c;
    std::vector<std::string>  names;
    std::vector<tax_t>       taxids;
    while((c = getopt(argc, argv, "n:t:h?")) >= 0) {
        switch(c) {
            case 'n': names.emplace_back(optarg); break;
            case 't': taxids.push_back(std::atoi(optarg)); break;
            case 'h': case '?': goto usage;
        }
            
    }
    if(argc < 3) {
        usage:
            std::fprintf(stderr, "Functionality: calculates and emits the taxid for the lca of all provided ids.\n"
                                 "Usage: %s <opts> tax_nodes.txt [optional: <nameidmap.txt>]\nFlags:\n-n:\tAdd a name from nameidmap to list.\n:-t:\t"
                                 "Add a taxid <int>.\n", argv[0]);
            std::exit(EXIT_FAILURE);
    }
    print_vec(names);
    khash_t(p) *tax(build_parent_map(argv[optind]));
    if(tax == nullptr) LOG_EXIT("Could not open taxmap. (See warning logs.)\n");
    khash_t(name) *name_hash(argv[optind + 1] ? build_name_hash(argv[optind + 1]) : nullptr);
    // std::cerr << "Name hash size: " << kh_size(name_hash) << '\n';
    if(name_hash) {
        khiter_t ki;
        const char *p;
        for(const auto &name: names) {
            std::cerr << "Name: " << name << '\n';
            if((ki = kh_get(name, name_hash, name.data())) != kh_end(name_hash)) {
                assert(ki != kh_end(name_hash));
                taxids.push_back(kh_val(name_hash, ki));
            }
            else {
                if((p = std::strchr(name.data(), '.'))) {
                    std::string trname(name.data(), p);
                    if((ki = kh_get(name, name_hash, trname.data())) != kh_end(name_hash)) {
                        assert(ki < kh_size(name_hash));
                        taxids.push_back(kh_val(name_hash, ki));
                    }
                    continue;
                }
                std::fprintf(stderr, "Warning: name '%s' not found. Ignoring.\n", name.data());
            }
        }
    } else {
        std::fprintf(stderr, "No name hash provided. (path: %s)\n", argv[optind + 1]);
    }
#if 0
    tax_t current_lca(taxids[0]);
    for(const auto id: taxids) {
        auto tmptax(lca(tax, id, current_lca));
        LOG_INFO("lca of %u and new %u is %u\n", current_lca, id, tmptax);
        current_lca = tmptax;
    }
    std::fprintf(stderr, "lca: %u\n", current_lca);
#else
    std::fprintf(stderr, "lca: %u\n", lca(tax, taxids));
#endif

    khash_destroy(tax);
    destroy_name_hash(name_hash);
}
