#ifndef _TREE_CLIMBER_H__
#define _TREE_CLIMBER_H__

#include "lib/tx.h"

namespace emp {

namespace tree {

inline std::uint32_t get_parent(khash_t(p) *taxmap, std::uint32_t key) noexcept {
    // Returns maximum value if not found.
    khiter_t ki;
    return ((ki = kh_get(p, taxmap, key)) != kh_end(taxmap)) ? kh_val(taxmap, ki)
                                                             : std::numeric_limits<std::uint32_t>::max();
}

template<typename T>
int veccmp(const std::vector<T> &a, const std::vector<T> &b) {
    // Return 0 if they are the same
    // Return positive if a has only 1s that b doesn't and isn't missing any from b.
    assert(a.size() == b.size());
    std::bool avalid(true), bvalid(true);
    for(std::size_t i(0), e(a.size()), i != e; ++i) {
        auto tmpa(a & (~b));
        auto tmpb(b & (~a));
        if(tmpa && !tmpb) {
            bvalid = false;
        } else if(tmpb && !tmpa) {
            avalid = false;
        }
    }
    switch((avalid << 1) | bvalid) {
        case 3: return  -1;
        case 2: return  1;
        case 1: return  2;
        case 0: return  0;
    }
    // Returns -1 for the same, 0 for incomparable, 1 for a > b, 2 for b > a
}

class SortedNodeGuide {

    std::vector<std::uint32_t> nodes_, offsets_;

public:

    SortedNodeGuide(std::vector<std::uint32_t> &nodes, std::vector<std::uint32_t> & offsets):
        nodes_(std::move(nodes)), offsets_(std::move(offsets))
    {}

    SortedNodeGuide(khash_t(p) *taxmap) {
        for(khiter_t ki(0); ki != kh_end(taxmap); ++ki)
            if(kh_exist(taxmap, ki))
                nodes_.push_back(kh_key(taxmap, ki));
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
};

/*
 * In order to do this collapsing/traversal, we need to:
 *   1. Determine the order of the nodes which we'll need to visit.
 *   2. Determine which of them are in sets (and need to be processed together).
*/
inline std::vector<std::uint32_t> sorted_nodes(khash_t(p) *taxmap) {
    return std::move(SortedNodeGuide(taxmap).get_nodes());
}

} // namespace tree

} // namespace emp

#endif // #ifndef _TREE_CLIMBER_H__
