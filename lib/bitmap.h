#ifndef _BITMAP_H__
#define _BITMAP_H__

#include "lib/tx.h"

namespace emp {

// Vector comparison function
template<typename T>
int veccmp(const std::vector<T> &a, const std::vector<T> &b);


template<typename T>
class AdjacencyList {
    std::unordered_map<const T*, std::vector<const T*>> map_;
    std::size_t m_, nelem_;
    bool is_reverse_;

public:

    auto find(const T *elem) const {return map_.find(elem);}
    auto find(const T &elem) const {return find(&elem);}
    auto end()         const {return map_.end();}
    AdjacencyList(): m_(0), nelem_(0), is_reverse_(0) {}
#if 0
    AdjacencyList(const AdjacencyList<T> &other) = default;
    AdjacencyList(AdjacencyList<T> &&other)      = default;
#endif

    AdjacencyList(count::Counter<T> &counts, bool reverse=false):
        m_(0), nelem_(counts.get_nelem()), is_reverse_(reverse) {
        // O(n^2)
        std::set<T*> ptrs;
        auto &map(counts.get_map());
        if(reverse == false) {
            for(auto i(map.cbegin()), ie(map.cend()); i != ie; ++i) {
                ++m_;
                auto j(i);
                while(++j != map.end()) {
                    typename std::unordered_map<T*, std::vector<T*>>::iterator m;
                    switch(veccmp(i->first, j->first)) {
                        case 1:
                            map_[&i->first].push_back(&j->first);
                            break;
                        case 2:
                            map_[&j->first].push_back(&i->first);
                            break;
                    }
                }
            }
        } else {
            for(auto i(map.cbegin()), ie(map.cend()); i != ie; ++i) {
                ++m_;
                auto j(i);
                while(++j != map.end()) {
                    typename std::unordered_map<T*, std::vector<T*>>::iterator m;
                    switch(veccmp(i->first, j->first)) {
                        case 1:
                            map_[&j->first].push_back(&i->first);
                            break;
                        case 2:
                            map_[&i->first].push_back(&j->first);
                            break;
                    }
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
    static const uint8_t ret[4]{3, 2, 1, 0};
    for(std::size_t i(0), e(a.size()); i != e; ++i) {
        auto tmpa(a[i] & (~b[i]));
        auto tmpb(b[i] & (~a[i]));
        if(tmpa && !tmpb) {
            bvalid = false;
        } else if(tmpb && !tmpa) {
            avalid = false;
        }
    }
    return ret[(avalid << 1) | bvalid];
}

class bitmap_t {
    std::unordered_map<std::uint64_t, std::vector<std::uint64_t>> core_;

public:
    std::unordered_map<std::uint64_t, std::vector<std::uint64_t>> fill(const kgset_t &set) {
        std::unordered_map<std::uint64_t, std::vector<std::uint64_t>> tmp;
        const unsigned len((set.size() + 63) >> 6);
        khash_t(all) *h;
        const auto &vec(set.get_core());
        LOG_DEBUG("Size of set: %zu\n", set.size());

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

    bitmap_t(){}
    bitmap_t(const kgset_t &set) {
        const auto tmp(fill(set));
        unsigned bitsum;
#if !NDEBUG
        std::size_t n_passed(0), total(tmp.size());
#endif
        for(auto &i: tmp)
            if((bitsum = popcnt::vec_popcnt(i.second)) != 1 &&
                bitsum != set.size())
#if !NDEBUG
                    core_.emplace(i.first, i.second), ++n_passed;
            else    LOG_DEBUG("bitsum is %u while the set size is %zu\n", bitsum, set.size());
#else
                    core_.emplace(i.first, i.second);
#endif
        LOG_DEBUG("Keeping %zu of %zu kmers for bit patterns which are not exactly compressed by the taxonomy heuristic.\n",
                  n_passed, total);
    }
    bitmap_t(bitmap_t &&other)            = default;
    bitmap_t &operator=(bitmap_t &&other) = default;

    count::Counter<std::vector<std::uint64_t>> to_counter();
};



std::uint64_t score_node_addn(const std::vector<std::uint64_t> &bitstring,
                              const adjmap_t &am, const count::Counter<std::vector<std::uint64_t>> &counts, std::size_t nelem);

}

#endif // #ifndef _BITMAP_H__
