#pragma once
#include "bitcmp.h"
#include "kgset.h"

namespace emp {


template<typename T>
class AdjacencyList {
    std::unordered_map<const T*, lazy::vector<const T*>, wang_hash_struct> map_;
    size_t m_, nelem_;
    bool is_reverse_;

public:

    auto find(const T *elem) const {return map_.find(elem);}
    auto find(const T &elem) const {return find(&elem);}
    auto end()         const {return map_.end();}
    AdjacencyList(): m_(0), nelem_(0), is_reverse_(0) {}
    enum orientation: size_t {
        FORWARD = 0,
        REVERSE = 1
    };
#if 0
    AdjacencyList(const AdjacencyList<T> &other) = default;
    AdjacencyList(AdjacencyList<T> &&other)      = default;
#endif

    AdjacencyList(const count::Counter<T> &counts, bool reverse=false):
        m_(0), nelem_(counts.get_nelem()), is_reverse_(reverse) {
        // O(n^2)
        std::set<T*> ptrs;
        const auto &map(counts.get_map());
        LOG_INFO("Trying to do stuff.\n");
        if(reverse == false) {
            for(auto i(map.cbegin()), ie(map.cend()); i != ie; ++i) {
                ++m_;
                auto j(i);
                while(++j != map.end()) {
                    LOG_DEBUG("Calling veccmp in forward\n");
                    switch(veccmp(i->first, j->first)) {
                        case 1:
                            map_[&i->first].push_back(&j->first); // i is a strict parent of j.
                            break;
                        case 2:
                            map_[&j->first].push_back(&i->first); // j is a strict parent of i.
                            break;
                    }
                    LOG_DEBUG("Called veccmp in forward\n");
                }
            }
        } else {
            for(auto i(map.cbegin()), ie(map.cend()); i != ie; ++i) {
                ++m_;
                auto j(i);
                while(++j != map.end()) {
                    LOG_DEBUG("Calling veccmp in reverse\n");
                    switch(veccmp(i->first, j->first)) {
                        case 1:
                            map_[&j->first].push_back(&i->first); // j is a strict parent of i.
                            break;
                        case 2:
                            map_[&i->first].push_back(&j->first); // i is a strict parent of j.
                            break;
                    }
                    LOG_DEBUG("Called veccmp in reverse\n");
                }
            }
        }
        LOG_DEBUG("About to shrink_to_fit.\n");
        for(auto &el: map_) {
            el.second.shrink_to_fit();
        }
        LOG_DEBUG("Finished constructor.");
    }
};

using adjmap_t = AdjacencyList<bitvec_t>;

class bitmap_t {
    std::unordered_map<u64, bitvec_t> core_;

public:
    std::unordered_map<u64, bitvec_t> fill(const kgset_t &set) {
        std::unordered_map<u64, bitvec_t> tmp;
        const unsigned len((set.size() + 63) >> 6);
        khash_t(all) *h;
        const auto &vec(set.get_core());
        LOG_DEBUG("Size of set: %zu\n", set.size());

        for(size_t i(0); i < set.size(); ++i) {
            h = vec[i];
            for(khiter_t ki(0); ki != kh_end(h); ++ki) {
                if(kh_exist(h, ki)) {
                    auto m(tmp.find(kh_key(h, ki)));
                    if(m == tmp.end()) m = tmp.emplace(kh_key(h, ki),
                                                       bitvec_t(len, lazy::LAZY_VEC_INIT)).first;
                    m->second[i >> 6] |= 1u << (i & 63u);
                }
            }
        }
        return tmp;
    }

    auto &get_map() {return core_;}
    auto &cget_map() const {return static_cast<const decltype(core_)&>(core_);}

    bitmap_t(){}
    bitmap_t(const kgset_t &set) {
        const auto tmp(fill(set));
        unsigned bitsum;
#if !NDEBUG
        size_t n_passed(0), total(tmp.size());
        std::unordered_map<unsigned, size_t> counts;
        std::unordered_set<std::string> stringset;
#endif
        // Only keeps kmers from kgset if they don't have 1 or set.size() bits set.
        for(auto &i: tmp) {
            bitsum = pop::vec_popcnt(i.second);
            if(bitsum != 1u && bitsum != set.size()) {
                core_.emplace(i.first, i.second);
#if !NDEBUG
                ++n_passed, stringset.insert(bitvec2str(i.second));
#endif
            }
#if !NDEBUG
            ++counts[bitsum];
#endif
        }
        LOG_DEBUG("Keeping %zu of %zu kmers whose bit patterns are not exactly compressed by the taxonomy heuristic.\n",
                  n_passed, total);
#if !NDEBUG
        for(const auto &pair: counts)
            LOG_DEBUG("Count %u appeared %u times\n", pair.first, pair.second);
        for(const auto &str: stringset) LOG_DEBUG("'%s' found\n", str.data());
#endif
        
    }
    bitmap_t(bitmap_t &&other)            = default;
    bitmap_t &operator=(bitmap_t &&other) = default;

    count::Counter<bitvec_t> to_counter();
};



u64 score_node_addn(const bitvec_t &bitstring,
                    const adjmap_t &am, const count::Counter<bitvec_t> &counts, size_t nelem);

} // namespace emp
