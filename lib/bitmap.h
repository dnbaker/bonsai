#ifndef _BITMAP_H__
#define _BITMAP_H__

#include "lib/tx.h"

namespace emp {

// Vector comparison function
template<typename T>
int veccmp(const std::vector<T> &a, const std::vector<T> &b);


template<typename T>
class AdjacencyList {
    std::unordered_map<T, std::vector<T*>> map_;
    std::size_t m_, nelem_;

public:

    auto find(T &elem) const {return map_.find(elem);}
    auto end()         const {return map_.end();}

    // Could also use a tree map depending on the size of the graph/number of keys.
    AdjacencyList(count::Counter<std::vector<std::uint64_t>> &counts):
        m_(0), nelem_(counts.get_nelem()) {
        std::set<T*> ptrs;
        const auto map(counts.get_map());
        for(auto i(map.cbegin()), ie(map.cend()); i != ie; ++i) {
            ++m_;
            auto j(i);
            while(++j != map.cend()) {
                typename std::unordered_map<T, std::vector<T*>>::iterator m;
                switch(veccmp(i->first, j->first)) {
                    case 1:
                        if((m = map_.find(i->first)) == map_.end())
                                map_[i->first].push_back(const_cast<std::vector<std::uint64_t> *>(&j->first));
                        else m->second.push_back(const_cast<std::vector<std::uint64_t> *>(&j->first));
                        break;
                    case 2:
                        if((m = map_.find(j->first)) == map_.end())
                                map_[j->first].push_back(const_cast<std::vector<std::uint64_t> *>(&i->first));
                        else m->second.push_back(const_cast<std::vector<std::uint64_t> *>(&i->first));
                        break;
                    case 0:        // They are identical: this shouldn't happen because of the way we iterate.
                    case 3: break; // Do nothing: neither is a strict parent of the other.
                }
            }
        }
    }
};

using adjmap_t = AdjacencyList<std::vector<std::uint64_t>>;


/*
 * We need to somehow build a tree of bitstrings:
 * I want to try two data structures.
 * First, I want to make a map from each bitstring to all of its descendents.
 */
template<typename T>
int veccmp(const std::vector<T> &a, const std::vector<T> &b) {
    // Return 0 if they are the same
    // Return positive if a has only 1s that b doesn't and isn't missing any from b.
    bool avalid(true), bvalid(true);
    for(std::size_t i(0), e(a.size()); i != e; ++i) {
        auto tmpa(a[i] & (~b[i]));
        auto tmpb(b[i] & (~a[i]));
        if(tmpa && !tmpb) {
            bvalid = false;
        } else if(tmpb && !tmpa) {
            avalid = false;
        }
    }
    switch((avalid << 1) | bvalid) {
        default: case 3: return -1;
        case 2: return  1;
        case 1: return  2;
        case 0: return  0;
    }
    return 0; // This never happens.
    // Returns -1 for the same, 0 for incomparable, 1 for a > b, 2 for b > a
}

class bitmap_t {
    std::unordered_map<std::uint64_t, std::vector<std::uint64_t>> core_;
    kgset_t &set_;

    public:

    std::unordered_map<std::uint64_t, std::vector<std::uint64_t>> fill(kgset_t &set) {
        std::unordered_map<std::uint64_t, std::vector<std::uint64_t>> tmp;
        const unsigned len((set.size() + 63) >> 6);
        khash_t(all) *h;
        auto &vec(set.get_core());

        for(std::size_t i(0); i < set.size(); ++i) {
            h = vec[i];
            for(khiter_t ki(0); ki != kh_end(h); ++ki) {
                if(kh_exist(h, ki)) {
                    auto m(tmp.find(kh_key(h, ki)));
                    if(m == tmp.end()) m = tmp.emplace(kh_key(h, ki),
                                                       std::move(std::vector<std::uint64_t>(len))).first;
                    m->second[i >> 6] |= 1u << (i & 63u);
                }
            }
        }
        return tmp;
    }

    bitmap_t(kgset_t &set): set_(set) {
        auto tmp(fill(set));
        unsigned bitsum;
#if !NDEBUG
        size_t n_passed(0);
        size_t total(tmp.size());
#endif
        for(auto &i: tmp)
            if((bitsum = popcnt::vec_popcnt(i.second)) != 1 &&
                bitsum != set.size())
#if !NDEBUG
                    core_.emplace(i.first, i.second), ++n_passed;
#else
                    core_.emplace(i.first, i.second);
#endif
        LOG_DEBUG("Keeping %zu of %zu kmers for bit patterns which are not exactly compressed by the taxonomy heuristic.\n",
                  n_passed, total);
    }

    count::Counter<std::vector<std::uint64_t>> to_counter();
};



std::uint64_t score_node_addn(std::vector<std::uint64_t> &bitstring,
                              adjmap_t &am, count::Counter<std::vector<std::uint64_t>> &counts);

}

#endif // #ifndef _BITMAP_H__
