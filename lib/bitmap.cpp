#include "bitmap.h"

namespace emp {

std::uint64_t adjlist_node_addn(std::vector<std::uint64_t> &bitstring,
                                adjmap_t &am, count::Counter<std::vector<std::uint64_t>> &counts) {
    const auto m(am.find(bitstring));
    const auto node(counts.find(bitstring));
    if(unlikely(m == am.end()))        return UINT64_C(-1);
    if(unlikely(node == counts.end())) return UINT64_C(-1);
    std::uint64_t ret(node->second);
    for(auto i: m->second) ret += counts.find(*i)->second;
    return ret;
}

std::uint64_t bitdiff_node_addn(std::vector<std::uint64_t> &bitstring,
                                count::Counter<std::vector<std::uint64_t>> &counts) {
    return counts.get_nelem() - popcnt::vec_popcnt(bitstring);
}

std::uint64_t score_node_addn(std::vector<std::uint64_t> &bitstring,
                              adjmap_t &am, count::Counter<std::vector<std::uint64_t>> &counts) {
            // Consider splitting this into two functions for clarity (bit diff and adjacency list traversal).
    const std::uint64_t node_addn(adjlist_node_addn(bitstring, am, counts));
    return node_addn == BF ? BF: node_addn * bitdiff_node_addn(bitstring, counts);
}


count::Counter<std::vector<std::uint64_t>>
bitmap_t::to_counter() {
    count::Counter<std::vector<std::uint64_t>> ret;
    for(auto &pair: core_) ret.add(pair.second);
    ret.set_nelem(set_.size());
    return ret;
}

} // namespace emp
