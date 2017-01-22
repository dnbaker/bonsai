#include "bitmap.h"

namespace emp {

std::uint64_t score_node_addn(std::vector<std::uint64_t> &bitstring,
                              adjmap_t &am, count::Counter<std::vector<std::uint64_t>> &counts) {
    std::uint64_t ret;
    const auto m(am.find(bitstring));
    const auto node(counts.find(bitstring));
    if(m == am.end())        goto fail;
    if(node == counts.end()) goto fail;
    ret = node->second;
    for(auto i: m->second) ret += counts.find(*i)->second;
            // Consider splitting this into two functions for clarity (bit diff and adjacency list traversal).
    return ret * (counts.get_nelem() - popcnt::vec_popcnt(bitstring));

    fail: return UINT64_C(-1);
}

adjmap_t adj_list(count::Counter<std::vector<std::uint64_t>> &counts) {
    const auto map(counts.get_map());
    adjmap_t ret;
    for(auto i(map.cbegin()), ie(map.cend()); i != ie; ++i) {
        auto j(i);
        while(++j != map.cend()) {
            adjmap_t::iterator m;
            switch(veccmp(i->first, j->first)) {
                case 0: break; // They are identical: this shouldn't happen because of the way we iterate.
                case 1: 
                    if((m = ret.find(i->first)) == ret.end())
                            ret[i->first].push_back(const_cast<std::vector<std::uint64_t> *>(&j->first));
                        //m = ret.emplace(i->first, std::vector<std::vector<std::uint64_t> *>{&j->first});
                    else m->second.push_back(const_cast<std::vector<std::uint64_t> *>(&j->first));
                    break;
                case 2:
#if 0
                    if((m = ret.find(j->first)) == ret.end())
                        ret.emplace(j->first, std::vector<std::vector<std::uint64_t> *>{&i->first});
                    else m->second.push_back(const_cast<std::vector<std::uint64_t> *>(&i->first));
#endif
                    break;
                case 3: break; // Do nothing: neither is a strict parent of the other.
            }
        }
    }
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
