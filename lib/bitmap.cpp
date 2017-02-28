#include "bitmap.h"
#include <cassert>

namespace emp {

std::uint64_t adjlist_node_addn(std::vector<std::uint64_t> &bitstring,
                                adjmap_t &am, count::Counter<std::vector<std::uint64_t>> &counts, size_t nelem) {
    assert(bitstring.size() << 6 >= nelem);
    const auto m(am.find(bitstring));
    const auto node(counts.find(bitstring));
    if(unlikely(m == am.end()) || node == counts.end()) return UINT64_C(-1);
    std::uint64_t ret(node->second * (nelem - popcnt::vec_popcnt(node->first)));
    for(const auto i: m->second) ret += counts.find(*i)->second * popcnt::vec_popcnt(*i);
    return ret;
}


count::Counter<std::vector<std::uint64_t>>
bitmap_t::to_counter() {
    count::Counter<std::vector<std::uint64_t>> ret;
    for(auto &pair: core_) ret.add(pair.second);
    ret.set_nelem(set_.size());
    return ret;
}

} // namespace emp
