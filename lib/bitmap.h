#ifndef _BITMAP_H__
#define _BITMAP_H__

#include "lib/tx.h"

namespace emp {

class bitmap_t {
    std::unordered_map<uint64_t, std::vector<std::uint64_t>> core_;
    kgset_t &set_;

    public:

    std::unordered_map<std::uint64_t, std::vector<std::uint64_t>> fill(kgset_t &set) {
        std::unordered_map<std::uint64_t, std::vector<std::uint64_t>> tmp;
        const unsigned len((set.size() + 63) >> 6);
        khash_t(all) *h;
        auto &vec(set.get_core());

        for(size_t i(0); i < set.size(); ++i) {
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
#if !NDEBUG
        size_t n_passed(0), total(tmp.size());
#endif
        unsigned bitsum;
        for(auto &i: tmp) {
            bitsum = popcnt::vec_popcnt(i.second);
            if(bitsum != 1 && bitsum != set.size()) {
#if !NDEBUG
                ++n_passed;
                //LOG_DEBUG("Number passed: %zu. bitsum: %u\n", n_passed, bitsum);
#endif
                core_.emplace(i.first, i.second);
            }
        }
        LOG_DEBUG("Keeping %zu of %zu kmers for bit patterns which are not exactly compressed by the taxonomy heuristic.\n",
                  n_passed, total);
    }
    count::Counter<std::vector<std::uint64_t>> to_counter() {
        count::Counter<std::vector<std::uint64_t>> ret;
        for(auto &pair: core_) ret.add(pair.second);
        ret.set_nelem(set_.size());
        return ret;
    }
};



}

#endif // #ifndef _BITMAP_H__
