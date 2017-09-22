#include "lib/util.h"
#include "getopt.h"

using namespace emp;

int main(int argc, char *argv[]) {
    int c;
    std::vector<std::string>  names;
    std::vector<tax_t>       taxids;
    while((c = getopt(argc, argv, "n:t:h?")) >= 0) {
        switch(c) {
            case 'n': names.emplace_back(argv[optind]); break;
            case 't': taxids.push_back(std::atoi(argv[optind])); break;
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
    const bool has_nameidmap(optind < argc);
    if(names.size() && !has_nameidmap) {
        LOG_EXIT("Error: nameidmap not provided but names <-n options> provided.\n");
    }
    khash_t(p) *tax(build_parent_map(argv[optind]));
    if(tax == nullptr) LOG_EXIT("Could not open taxmap. (See warning logs.)\n");
    khash_t(name) *name_hash(argv[optind + 1] ? build_name_hash(argv[optind + 1]): nullptr);
    if(name_hash) {
        khiter_t ki;
        destroy_name_hash(name_hash);
        for(const auto &name: names) {
            if((ki = kh_get(name, name_hash, name.data())) == kh_end(name_hash)) taxids.push_back(kh_val(name_hash, ki));
            else std::fprintf(stderr, "Warning: name '%s' not found. Ignoring.\n", name.data());
        }
    }
    std::fprintf(stderr, "lca: %u\n", lca(tax, taxids));
    khash_destroy(tax);
}
