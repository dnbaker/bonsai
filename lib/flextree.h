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
    std::vector<NodeType *> parents_;
    const std::uint32_t          pc_;  // cached popcount of bitvector
    const std::uint32_t          bc_;  // bitcount in family
    const std::uint32_t          si_;

    fnode_t(const bitvec_t &bits, std::uint32_t bc, std::uint32_t subtree_index, const std::uint64_t n=0):
        n_{n}, laa_{nullptr}, pc_{static_cast<std::uint32_t>(bits.size())},
        bc_{bc},  si_{subtree_index} {}

    bool added() const {return laa_ && &laa_->second == this;}
};

INLINE std::uint64_t get_score(const NodeType &node) {
    if(node.second.added()) return 0;
    std::uint64_t ret(node.second.n_);
    for(auto s: node.second.subsets_)
        if(s->second.added() == false)
            if(s->second.laa_ == nullptr ||
               (s->second.laa_ != s && s->second.laa_->second.pc_ > node.second.pc_))
                ret += s->second.n_;
    return ((node.second.laa_ ? node.second.laa_->second.pc_
                              : node.second.bc_) - node.second.pc_) * ret;
}

struct node_lt {
    bool operator()(const NodeType *a, const NodeType *b) const {
        return get_score(*a) > get_score(*b);
    }
};

class FlexMap {

    std::unordered_map<bitvec_t, fnode_t> map_;
    std::uint64_t                           n_;
    std::uint32_t                    bitcount_;
    const std::uint32_t                    id_;

public:
    void prepare_data() {
        build_adjlist();
    }
    void build_adjlist() {
        for(auto i(map_.begin()), ie(map_.end()); i != ie; ++i) {
            auto j(i);
            while(++j != map_.end()) {
                switch(veccmp(i->first, j->first)) {
                    case 1:
                        i->second.subsets_.emplace_back(&*j);
                        j->second.parents_.emplace_back(&*i);
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
                        i->second.parents_.emplace_back(&*j);
                        break;
                }
            }
        }
    }
    void add_to_heap(std::set<const NodeType *, node_lt> &heap) const {
        for(auto &pair: map_) {
            heap.insert(&pair);
        }
    }
public:
    FlexMap(std::uint32_t id, std::uint32_t bitcount): n_{0}, bitcount_{bitcount}, id_{id} {
    }

    INLINE void add(bitvec_t &&elem) {
#if __GNUC__ >= 7
        if(auto match = map_.find(elem); match == map_.end())
#else
        auto match(map_.find(elem));
        if(match == map_.end())
#endif
            map_.emplace(std::move(elem),
                         std::forward<fnode_t>(fnode_t(elem, bitcount_, id_, UINT64_C(1))));
        else ++match->second.n_;
        ++n_;
    }
    template<typename T>
    void process(const T &container, const khash_t(name) *name_hash) {
        const size_t nelem(size(container));
        std::map<std::string, tax_t> pathmap;
        std::unordered_map<tax_t, unsigned> taxmap;
        for(const auto &el: container) {
            const tax_t tax(get_taxid(el.data(), name_hash));
            pathmap.emplace(el, tax);
        }
    }
};

class FMEmitter {
    std::vector<FlexMap>      subtrees_;
    std::set<NodeType *, node_lt> heap_;
    std::unordered_set<tax_t>    added_;
    khash_t(p)              *const tax_;
    // Need taxid to paths map
    // 
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
            if(bptr->second.added()) {
                LOG_WARNING("Cannot add more nodes. [Best candidate is impossible.] Breaking from loop.\n");
                break;
            } else bptr->second.laa_ = bptr;
#if NDEBUG
           static_assert(false, "raise NotImplementedError('Take nodes from parents of affected nodes out of the tree and put them back in.')");
#endif
            for(const auto parent: bptr->second.parents_) {
                to_reinsert.insert(parent);
            }
            for(auto other: bptr->second.subsets_) {
                to_reinsert.insert(other);
                if(other->second.laa_ == nullptr ||
                   other->second.laa_->second.pc_ > bptr->second.pc_)
                    other->second.laa_ = bptr;
                for(const auto parent: other->second.parents_) {
                    to_reinsert.insert(parent);
                }
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
            assert_sorted<decltype(heap_), node_lt>(heap_);
        }
        std::fwrite(ks.data(), 1, ks.size(), fp);
        ks.clear();
    }
    FlexMap &emplace_subtree(unsigned ngenomes) {
#if __cplusplus < 201700LL
        subtrees_.emplace_back(ngenomes, subtrees_.size());
        return subtrees_.back();
#else
        return subtrees_.emplace_back(ngenomes, subtrees_.size());
#endif
    }
    void format_emitted_node(ks::KString &ks, const NodeType *node) const {
        
    }
    // Also need a map of taxid to tax level.
    FMEmitter(khash_t(p) *tax): tax_{tax} {}
};

} // namespace emp

#endif // __FLEX_TREE_H
