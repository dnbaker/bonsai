#ifndef _TREE_CLIMBER_H__
#define _TREE_CLIMBER_H__

#include "lib/tx.h"

namespace emp {

namespace tree {

class SortedNodeGuide {
    std::vector<std::uint32_t> nodes_, offsets_;
public:
    SortedNodeGuide(std::vector<std::uint32_t> &nodes, std::vector<std::uint32_t> & offsets):
        nodes_(std::move(nodes)), offsets_(std::move(offsets))
    {}
};

namespace {
    uint32_t get_parent(khash_t(p) *taxmap, uint32_t key) {
        khiter_t ki;
        return ((ki = kh_get(p, taxmap, key)) != kh_end(taxmap)) ? kh_val(taxmap, ki): std::numeric_limits<uint32_t>::max();
    }
}

/*
 * In order to do this collapsing/traversal, we need to:
 *   1. Determine the order of the nodes which we'll need to visit.
 *   2. Determine which of them are in sets (and need to be processed together).
*/
SortedNodeGuide sorted_nodes(khash_t(p) *taxmap) {
    struct Node {
        const std::uint32_t id_, depth_;
        Node(std::uint32_t id, std::uint32_t depth): id_(id), depth_(depth) {}
    };
    std::vector<Node> nodes;
    nodes.reserve(taxmap->n_occupied);
    for(khiter_t ki(0); ki != kh_end(taxmap); ++ki)
        if(kh_exist(taxmap, ki))
            nodes.emplace_back(kh_key(taxmap, ki),
                               node_depth(taxmap, kh_key(taxmap, ki)));
    std::vector<std::uint32_t> ids;
    for(khiter_t ki(0); ki != kh_end(taxmap); ++ki) if(kh_exist(taxmap, ki)) ids.push_back(kh_key(taxmap, ki));
    std::sort(std::begin(ids), std::end(ids), [taxmap] (const std::uint32_t a, const std::uint32_t b) {
        const uint32_t nda(node_depth(taxmap, a)), ndb(node_depth(taxmap, b));
        if(nda != ndb) return nda > ndb;
        return get_parent(taxmap, a) < get_parent(taxmap, b);
    });
#if 0
    std::sort(std::begin(nodes), std::end(nodes), [taxmap] (const Node &a, const Node &b) {
        if(a.depth_ != b.depth_) {
            return a.depth_ > b.depth_;
        }
    });
#endif
    std::vector<std::uint32_t> offsets;
    unsigned u, last(std::numeric_limits<uint32_t>::max());
    for(size_t i(0); i < ids.size(); ++i) {
        if((u = node_depth(taxmap, ids[i])) != last) {
            offsets.push_back(i);
            last = u;
        }
    }
    return SortedNodeGuide(ids, offsets);
}
} // namespace tree

} // namespace emp

#endif // #ifndef _TREE_CLIMBER_H__
