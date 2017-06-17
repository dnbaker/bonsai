#ifndef __FLEX_TREE_H
#define __FLEX_TREE_H

#include "lib/bits.h"
#include "lib/counter.h"

namespace emp {

struct fnode_t;
using popcnt::vec_bitdiff;
using popcnt::vec_popcnt;
using NodeType = std::pair<const bitvec_t, fnode_t>;

struct fnode_t {
    std::uint64_t                 n_;  // Number of kmers at this point in tree.
    const NodeType             *laa_;  // lowest added ancestor
    std::vector<NodeType *> subsets_;
    const std::uint32_t          bc_;

    fnode_t(std::uint32_t bc, const std::uint64_t n):
        n_{n}, laa_{nullptr}, bc_{bc} {}

    bool added() const {return laa_ && &laa_->second == this;}
};


INLINE std::uint64_t get_score(const NodeType &node) {
    if(node.second.added()) return 0;
    std::uint64_t ret(node.second.n_);
    for(auto s: node.second.subsets_) {
        if(s->second.added()) continue;
        if(s->second.laa_ == nullptr) ret += s->second.n_;
        else if(s->second.laa_ != s &&
                vec_bitdiff(s->first, node.first) <
                    vec_bitdiff(s->first, s->second.laa_->first)) {
            ret += s->second.n_;
        }
    }
    return node.second.laa_ ? vec_bitdiff(node.second.laa_->first, node.first) * ret
                            : node.second.bc_ - vec_popcnt(node.first) * ret;
}

class FlexMap {

    struct node_lt {
        bool operator()(const NodeType *a, const NodeType *b) const {
            return get_score(*a) > get_score(*b);
        }
    };

    std::unordered_map<bitvec_t, fnode_t> map_;
    std::uint64_t                                           n_;
    std::uint32_t                                           bitcount_;

public:
    void prepare_data() {
        build_adjlist();
        fill_heap();
    }
    void build_adjlist() {
        for(auto i(map_.begin()), ie(map_.end()); i != ie; ++i) {
            auto j(i);
            while(++j != map_.end()) {
                switch(veccmp(i->first, j->first)) {
                    case 1:
                        i->second.subsets_.emplace_back(&*j);
#if !NDEBUG
                    {
                        for(auto ii(i->first.cbegin()), ji(j->first.cbegin()), ei(i->first.cend()); ii != ei; ++ii, ++ji) {
                            assert((*ii & (~*ji)));  // Assert that ii has bits set ji does not.
                            assert(!(*ji & (~*ii))); // Assert that ji has no bits set which are not set in i.
                                                     // If both tests pass, then we did this correctly.
                        }
                    }
#endif
                        break;
                    case 2:
                        j->second.subsets_.emplace_back(&*i);
                        break;
                }
            }
        }
    }
    void add_to_heap(std::set<NodeType *, node_lt> &heap) const {
        for(auto &pair: map_) {
            heap.insert(&pair);
        }
    }
public:
    FlexMap(std::uint32_t bitcount): n_{0}, bitcount_{bitcount} {
    }

    INLINE void add(bitvec_t &&elem) {
#if __GNUC__ >= 7
        if(auto match = map_.find(elem); match == map_.end())
#else
        auto match(map_.find(elem));
        if(match == map_.end())
#endif
            map_.emplace(std::move(elem),
                         std::forward<fnode_t>(fnode_t(bitcount_, UINT64_C(1))));
        else ++match->second.n_;
        ++n_;
    }
};

class FMEmitter {
    std::vector<FlexMap>      subtrees_;
    std::set<NodeType *, node_lt> heap_;
    khash_t(p)                    *tax_;
    void run_collapse(std::FILE* fp=stdout, std::size_t nelem=0) {
        assert(heap_.size());
        if(nelem == 0) {
            auto bccpy(kh_size(tax_));
            kroundup32(bccpy);
            ++bccpy;
            kroundup32(bccpy); //
            nelem = bccpy - kh_size(tax_);
            LOG_DEBUG("elements to add defaulting to %zu more [nodes addable before reaching 2 * nearest power of two (%zu).]\n", nelem, static_cast<std::size_t>(bccpy));
        }
        std::unordered_set<NodeType *> to_reinsert;
        ks::KString ks;
        for(std::size_t i(0); i < nelem; ++i) {
            const auto bptr(*heap_.begin());
            if((bptr->second.added()) {
                LOG_WARNING("Cannot add more nodes. [Best candidate is impossible.] Breaking from loop.");
                break;
            } else bptr->second.laa_ = bptr;
            for(auto other: bptr->second.subsets_) {
                to_reinsert.insert(other);
            }
            //Emit results
            format_emitted_node(ks, bptr);
            if(ks.size() > 1 << 16) {
                std::fwrite(ks.data(), 1, ks.size(), fp);
                ks.clear();
            }
            
            // Make a list of all pointers to remove and reinsert to the map.
            heap_.erase(heap_.begin());
            for(const auto el: to_reinsert) heap_.erase(el);
            for(const auto el: to_reinsert) heap_.insert(el);
            to_reinsert.clear();
            heap_.insert(bptr);
        }
        std::fwrite(ks.data(), 1, ks.size(), fp);
        ks.clear();
    }
    void format_emitted_node(ks::KString &ks, NodeType *node) {
        
    }
    FMEmitter(khash_t(p) *tax): tax_{tax} {}
};

} // namespace emp

#endif // __FLEX_TREE_H
