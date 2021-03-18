#ifndef BONSAI_SETSKETCHINDEX_H__
#define BONSAI_SETSKETCHINDEX_H__
#include "sketch/setsketch.h"

namespace bns {

template<typename FT=double, typename KeyT=uint64_t, typename IdT=uint32_t>
struct SetSketchIndex {
private:
    size_t m_;
    using SSType = sketch::CSetSketch<FT>;
    using HashMap = ska::flat_hash_map<KeyT, std::vector<IdT>>;
    using HashV = std::vector<HashMap>;
    std::vector<HashV> packed_maps_;
    std::vector<uint64_t> regs_per_reg_;
    size_t total_ids_ = 0;
public:
    size_t m() const {return m_;}
    size_t size() const {return total_ids_;}
    SetSketchIndex(size_t m, bool densified=false): m_(m) {
        uint64_t rpr = 1;
        const size_t nrpr = densified ? m: size_t(ilog2(sketch::integral::roundup(m)));
        regs_per_reg_.reserve(nrpr);
        packed_maps_.reserve(nrpr);
        for(;rpr <= m_;) {
            regs_per_reg_.push_back(rpr);
            packed_maps_.emplace_back(HashV(m_ / rpr));
            if(densified) {
                ++rpr;
            } else {
                rpr <<= 1;
            }
        }
    }
    template<typename Sketch>
    void update(const Sketch &item) {
        if(item.size() != m_) throw std::invalid_argument(std::string("Item has wrong size: ") + std::to_string(item.size()) + ", expected" + std::to_string(m_));
        const auto my_id = total_ids_++;
        const size_t n_subtable_lists = regs_per_reg_.size();
        using ResT = typename std::decay_t<decltype(item[0])>;
        for(size_t i = 0; i < n_subtable_lists; ++i) {
            auto &subtab = packed_maps_[i];
            const size_t nsubs = subtab.size();
            const size_t nelem = regs_per_reg_[i];
            for(size_t j = 0; j < nsubs; ++j) {
                KeyT myhash = XXH3_64bits(&item[nelem * j], nelem * sizeof(ResT));
                subtab[j][myhash].push_back(my_id);
            }
        }
    }
    template<typename Sketch>
    std::pair<std::vector<IdT>, std::vector<uint32_t>>
    query_candidates(const Sketch &item, size_t maxcand) const {
        /*
        *  Returns ids matching input minhash sketches, in order from most specific/least sensitive
        *  to least specific/most sensitive
        *  Can be then used, along with sketches, to select nearest neighbors
        */
        using ResT = typename std::decay_t<decltype(item[0])>;
        ska::flat_hash_set<IdT> rset; rset.reserve(maxcand);
        std::vector<IdT> passing_ids;
        std::vector<uint32_t> items_per_row;
        for(std::ptrdiff_t i = regs_per_reg_.size();--i >= 0;) {
            auto &m = packed_maps_[i];
            const size_t nelem = regs_per_reg_[i];
            const size_t nsubs = m.size();
            const size_t items_before = passing_ids.size();
            for(size_t j = 0; j < nsubs; ++j) {
                KeyT myhash = XXH3_64bits(&item[nelem * j], nelem * sizeof(ResT));
                auto it = m[j].find(myhash);
                if(it == m[j].end()) continue;
                for(const auto id: it->second) {
                    auto rit2 = rset.find(id);
                    if(rit2 == rset.end()) {
                        rset.insert(id);
                        passing_ids.push_back(id);
                    }
                }
            }
            items_per_row.push_back(passing_ids.size() - items_before);
            if(rset.size() >= maxcand) {
                std::fprintf(stderr, "Candidate set of size %zu, > maxc %zu\n", rset.size(), maxcand);
                break;
            }
        }
        return std::make_pair(passing_ids, items_per_row);
    }
};

} // namespace bns

#endif /* BONSAI_SETSKETCHINDEX_H__ */
