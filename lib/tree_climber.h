#ifndef TREE_CLIMBER_H__
#define TREE_CLIMBER_H__

#include "lib/tx.h"
#include "lib/database.h"
#include <shared_mutex>

#if __cplusplus < 201402L || __GNUC__ < 6
#define shared_mutex shared_timed_mutex
#endif

namespace emp {

namespace tree {



khash_t(p) *pruned_taxmap(std::vector<std::string> &paths, khash_t(p) *taxmap, khash_t(name) *name_hash);

class SortedNodeGuide {

    std::vector<tax_t> nodes_;
    std::vector<std::uint32_t> offsets_;
    std::vector<tax_t> parent_map_;

public:

    SortedNodeGuide(std::vector<tax_t> &nodes, std::vector<std::uint32_t> & offsets):
        nodes_(std::move(nodes)), offsets_(std::move(offsets))
    {}

    SortedNodeGuide(khash_t(p) *taxmap) {
        std::unordered_map<tax_t, std::uint32_t> tmp;
        nodes_.reserve(kh_size(taxmap));
        for(khiter_t ki(0); ki != kh_end(taxmap); ++ki) {
            if(kh_exist(taxmap, ki)) {
                nodes_.push_back(kh_key(taxmap, ki));
                tmp.emplace(kh_key(taxmap, ki), nodes_.size() - 1); // Convert this back into proper parent map.
            }
        }
        std::sort(std::begin(nodes_), std::end(nodes_), [taxmap] (const tax_t a, const tax_t b) {
            tax_t aa(node_depth(taxmap, a)), bb(node_depth(taxmap, b));
            if(aa != bb) return aa > bb;
            if((aa = get_parent(taxmap, a)) != // Set and compare lexicographically by parents.
               (bb = get_parent(taxmap, b)))
                return aa < bb;
            return a < b;
        });
        tax_t u, last(std::numeric_limits<tax_t>::max());
        for(std::size_t i(0), e(nodes_.size()); i < e; ++i)
            if((u = node_depth(taxmap, nodes_[i])) != last)
                offsets_.push_back(i), last = u;
    }

    const std::vector<tax_t> &get_nodes() { return nodes_;}
    const std::vector<std::uint32_t> &get_offsets() { return offsets_;}
    int lt(tax_t a, tax_t b) {
        return -1;
    }
};

/*
 * In order to do this collapsing/traversal, we need to:
 *   1. Determine the order of the nodes which we'll need to visit.
 *   2. Determine which of them are in sets (and need to be processed together).
*/
std::vector<std::string> invert_lca_map(Database<khash_t(c)> &db, const char *folder, int prebuilt=0);

std::vector<std::uint64_t> load_binary_kmers(const char *path);
khash_t(all) *load_binary_kmerset(const char *path);
std::vector<std::string> par_invert(Database<khash_t(c)> &db, const char *folder, int num_threads=16, std::size_t chunk_size=1<<16);

} // namespace tree

} // namespace emp

#endif // #ifndef TREE_CLIMBER_H__
