#include "bitmap.h"
#include <cassert>

namespace emp {

u64 score_node_addn(const bitvec_t &bitstring,
                              const adjmap_t &am, const count::Counter<bitvec_t> &counts, std::size_t nelem) {
    assert(bitstring.size() << 6 >= nelem);
    const auto m(am.find(bitstring));
    const auto node(counts.find(bitstring));
    if(unlikely(m == am.end()) || node == counts.end()) return UINT64_C(-1);
    u64 ret(node->second * (nelem - popcnt::vec_popcnt(node->first)));
    for(const auto i: m->second) ret += counts.find(*i)->second * popcnt::vec_popcnt(*i);
    return ret;
}


count::Counter<bitvec_t>
bitmap_t::to_counter() {
    count::Counter<bitvec_t> ret;
    for(auto &pair: core_) ret.add(pair.second);
    ret.set_nelem(core_.size());
    return ret;
}

} // namespace emp
