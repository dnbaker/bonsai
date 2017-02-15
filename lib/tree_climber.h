#ifndef TREE_CLIMBER_H__
#define TREE_CLIMBER_H__

#include "lib/tx.h"
#include "lib/database.h"

namespace emp {

namespace tree {

INLINE std::uint32_t get_parent(khash_t(p) *taxmap, std::uint32_t key) noexcept {
    // Returns maximum value if not found.
    khiter_t ki;
    return ((ki = kh_get(p, taxmap, key)) != kh_end(taxmap)) ? kh_val(taxmap, ki)
                                                             : std::numeric_limits<std::uint32_t>::max();
}


khash_t(p) *pruned_taxmap(std::vector<std::string> &paths, khash_t(p) *taxmap, khash_t(name) *name_hash);

class SortedNodeGuide {

    std::vector<std::uint32_t> nodes_, offsets_;
    std::vector<std::uint32_t> parent_map_;

public:

    SortedNodeGuide(std::vector<std::uint32_t> &nodes, std::vector<std::uint32_t> & offsets):
        nodes_(std::move(nodes)), offsets_(std::move(offsets))
    {}

    SortedNodeGuide(khash_t(p) *taxmap) {
        std::unordered_map<std::uint32_t, std::uint32_t> tmp;
        nodes_.reserve(kh_size(taxmap));
        for(khiter_t ki(0); ki != kh_end(taxmap); ++ki) {
            if(kh_exist(taxmap, ki)) {
                nodes_.push_back(kh_key(taxmap, ki));
                tmp.emplace(kh_key(taxmap, ki), nodes_.size() - 1); // Convert this back into proper parent map.
            }
        }
        std::sort(std::begin(nodes_), std::end(nodes_), [taxmap] (const std::uint32_t a, const std::uint32_t b) {
            std::uint32_t aa(node_depth(taxmap, a)), bb(node_depth(taxmap, b));
            if(aa != bb) return aa > bb;
            if((aa = get_parent(taxmap, a)) != // Set and compare lexicographically by parents.
               (bb = get_parent(taxmap, b)))
                return aa < bb;
            return a < b;
        });
        std::uint32_t u, last(std::numeric_limits<std::uint32_t>::max());
        for(std::size_t i(0), e(nodes_.size()); i < e; ++i)
            if((u = node_depth(taxmap, nodes_[i])) != last)
                offsets_.push_back(i), last = u;
    }

    const std::vector<std::uint32_t> &get_nodes() { return nodes_;}
    const std::vector<std::uint32_t> &get_offsets() { return offsets_;}
    int lt(std::uint32_t a, std::uint32_t b) {
        return -1;
    }
};

/*
 * In order to do this collapsing/traversal, we need to:
 *   1. Determine the order of the nodes which we'll need to visit.
 *   2. Determine which of them are in sets (and need to be processed together).
*/
std::vector<std::string> invert_lca_map(Database<khash_t(c)> &db, const char *folder);

std::vector<std::uint64_t> load_binary_kmers(const char *path);

} // namespace tree

} // namespace emp

#endif // #ifndef TREE_CLIMBER_H__
