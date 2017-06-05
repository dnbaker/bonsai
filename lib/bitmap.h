#ifndef _BITMAP_H__
#define _BITMAP_H__

#include "lib/tx.h"


namespace emp {

inline constexpr int log2_64(uint64_t value)
{
    // https://stackoverflow.com/questions/11376288/fast-computing-of-log2-for-64-bit-integers
    const int tab64[64] = {
        63,  0, 58,  1, 59, 47, 53,  2,
        60, 39, 48, 27, 54, 33, 42,  3,
        61, 51, 37, 40, 49, 18, 28, 20,
        55, 30, 34, 11, 43, 14, 22,  4,
        62, 57, 46, 52, 38, 26, 32, 41,
        50, 36, 17, 19, 29, 10, 13, 21,
        56, 45, 25, 31, 35, 16,  9, 12,
        44, 24, 15,  8, 23,  7,  6,  5
    };
    value |= value >> 1;
    value |= value >> 2;
    value |= value >> 4;
    value |= value >> 8;
    value |= value >> 16;
    value |= value >> 32;
    return tab64[((uint64_t)((value - (value >> 1))*0x07EDD5E59A4E28C2)) >> 58];
}

// Vector comparison function
template<typename T>
int veccmp(const std::vector<T> &a, const std::vector<T> &b);
template<typename T, typename size_type=std::size_t>
struct ptr_hash {
    static constexpr size_t SHIFT = log2_64(alignof(T));
    size_type operator()(const T *a) const {
        return reinterpret_cast<size_type>(a) >> SHIFT;
    }
};


template<typename T>
class AdjacencyList {
    std::unordered_map<const T*, std::vector<const T*>, ptr_hash<T>> map_;
    std::size_t m_, nelem_;
    bool is_reverse_;

public:

    auto find(const T *elem) const {return map_.find(elem);}
    auto find(const T &elem) const {return find(&elem);}
    auto end()         const {return map_.end();}
    AdjacencyList(): m_(0), nelem_(0), is_reverse_(0) {}
    enum orientation: std::size_t {
        FORWARD = 0,
        REVERSE = 1
    };
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
                    switch(veccmp(i->first, j->first)) {
                        case 1:
                            map_[&i->first].push_back(&j->first); // i is a strict parent of j.
                            break;
                        case 2:
                            map_[&j->first].push_back(&i->first); // j is a strict parent of i.
                            break;
                    }
                }
            }
        } else {
            for(auto i(map.cbegin()), ie(map.cend()); i != ie; ++i) {
                ++m_;
                auto j(i);
                while(++j != map.end()) {
                    switch(veccmp(i->first, j->first)) {
                        case 1:
                            map_[&j->first].push_back(&i->first); // j is a strict parent of i.
                            break;
                        case 2:
                            map_[&i->first].push_back(&j->first); // i is a strict parent of j.
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
    for(std::size_t i(0), e(a.size()); i != e; ++i) {
        // First term/second bit: a does not have any bits set b does not.
        // Second term/first bit: b does not have any bits set a does not.
        // If both are 1, then they are equal on this word, so we do nothing.
        // If both are 0, then neither is a strict parent, return 3.
        switch(((!(a[i] & (~b[i]))) << 1) | !(b[i] & (~a[i]))) {
            case 2: avalid = false; break;
            case 1: bvalid = false; break;
            case 0: return 3;
        }
    }
    static const uint8_t ret[4]{3, 2, 1, 0};
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
            else if(bitsum > 1 && bitsum != set.size()) LOG_DEBUG("bitsum is %u while the set size is %zu\n", bitsum, set.size());
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
