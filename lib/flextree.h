#ifndef __FLEX_TREE_H
#define __FLEX_TREE_H

#include "lib/util.h"
#include "lib/bits.h"
#include "lib/counter.h"

namespace emp {

struct fnode_t;
using popcnt::vec_bitdiff;
using popcnt::vec_popcnt;
using NodeType = std::pair<const bitvec_t, fnode_t>;
using namespace std::literals;

struct fnode_t {
    std::uint64_t                 n_;  // Number of kmers at this point in tree.
    const NodeType             *laa_;  // lowest added ancestor
    std::vector<NodeType *> subsets_;
    std::vector<NodeType *> parents_;
    const std::uint32_t          pc_;  // cached popcount of bitvector
    const std::uint32_t          bc_;  // bitcount in family
    const std::uint32_t          si_;

    fnode_t(const bitvec_t &bits, std::uint32_t bc, std::uint32_t subtree_index, const std::uint64_t n=0):
        n_{n}, laa_{nullptr}, pc_{static_cast<std::uint32_t>(popcnt::vec_popcnt(bits))},
        bc_{bc}, si_{subtree_index} {}

    bool added()      const {return laa_ && &laa_->second == this;}
    std::string str() const {
        return "[fnode]{n:"s + std::to_string(n_) + (added() ? ",added:true,": ",added:false,bc:")
                             + std::to_string(bc_) + ",pc:" + std::to_string(pc_)
                             + ",subtree id:" + std::to_string(si_) + '}';
    }
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
    std::vector<tax_t>                    tax_;
    std::uint64_t                           n_;
    std::uint32_t                    bitcount_;
    const std::uint32_t                    id_;
    const tax_t                        parent_;

public:
    std::string str() const {
        ks::KString ks;
        ks.sprintf("FlexMap:{id:%u;parent:%u;bitcount:%zu;n: %" PRIu64 ";map_(%zu);taxes(%zu): [", id_, parent_, bitcount_, n_, map_.size(), tax_.size());
        for(const auto tax: tax_) ks.putuw_(tax), ks.putc_(',');
        ks.back() = ']';
        ks.putc('}');
        std::string ret(ks.data());
        return ret;
    }
    auto parent()           const {return parent_;}
    const auto &get_taxes() const {return tax_;}
    auto num_descendents()  const {return bitcount_;}
    void prepare_data() {
        build_adjlist();
    }
    void build_adjlist() {
        LOG_DEBUG("Building adjacency list with map size %zu\n", map_.size());
        for(auto i(map_.begin()), ie(map_.end()); i != ie; ++i) {
            auto j(i);
            while(++j != map_.end()) {
                assert(j->first.size() == i->first.size());
                switch(veccmp(i->first, j->first)) {
                    case BitCmp::FIRST_PARENT:
                        i->second.subsets_.emplace_back(&*j);
                        j->second.parents_.emplace_back(&*i);
                        break;
                    case BitCmp::SECOND_PARENT:
                        j->second.subsets_.emplace_back(&*i);
                        i->second.parents_.emplace_back(&*j);
                        break;
                    case BitCmp::EQUAL: LOG_EXIT("We should never have two identical bitmaps in the same hashmap.\n"); // Do nothing
                    case BitCmp::INCOMPARABLE: break;
                }
            }
        }
    }
    template<typename T, typename LT>
    void add_to_heap(std::set<T, LT> &heap) const {
#if 0
        for(typename std::set<T, LT>::iterator it(map_.begin()), eit(map_.end()); it != eit; ++it) {
            heap.insert(&*it);
        }
#endif
        if(map_.size() == 0) {
            LOG_WARNING("Heap of string '%s' is empty....\n", this->str());
        } else { // Redundant but explains intention.
            for(auto &pair: map_) {
                std::pair<const bitvec_t, fnode_t> &ref((std::pair<const bitvec_t, fnode_t> &)pair);
                LOG_DEBUG("Adding pair where key's first entry is is %" PRIu64 " and the node is %s\n", ref.first[0], ref.second.str().data());
                //heap.insert(const_cast<std::pair<const std::vector<long long unsigned int>, emp::fnode_t>*>(&pair));
                heap.insert(&ref);
            }
        }
    }
public:
    FlexMap(const tax_t parent, const std::uint32_t ntaxes, std::uint32_t id):
        n_{0}, bitcount_{ntaxes}, id_{id}, parent_{parent} {}

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
    void fill(const std::unordered_map<tax_t, strlist> &list, const Spacer &sp, int num_threads=-1,
              khash_t(all) *acc=nullptr) {
        if(tax_.size()) tax_.clear();
        for(const auto &pair: list) {
            LOG_DEBUG("Filling from flexmap #%u with tax %u and is %s\n", id_, pair.first, pair.second.empty() ? "empty": "not empty");
            tax_.push_back(pair.first);
        }
#if !NDEBUG
        bitmap_t tmp_bitmap(kgset_t(list, sp, num_threads, acc));
        auto &map(tmp_bitmap.get_map());
        LOG_DEBUG("Map size: %zu\n", map.size());
        for(auto &&pair: map) add(std::move(pair.second));
#else
        for(auto &&pair: bitmap_t(kgset_t(list, sp, num_threads, acc)).get_map()) add(std::move(pair.second));
#endif
    }
};

class FMEmitter {
    std::vector<FlexMap>                 subtrees_;
    std::set<NodeType *, node_lt>            heap_;
    std::unordered_set<tax_t>               added_;
    khash_t(p)                         *const tax_;
    const std::unordered_map<tax_t, strlist> &tpm_;


    /*
     * Emits an additional node to the tree where its paren
    */
    void format_emitted_node(ks::KString &ks, const NodeType *node, const std::uint64_t score, const tax_t taxid) const {
        std::uint64_t val;
        const auto &fm(subtrees_[node->second.si_]);
        ks.putuw_(taxid);
        ks.putc_('\t');
        ks.putl_(score);
        ks.putc_('\t');
        ks.putuw_(fm.parent());
        const auto &taxes(fm.get_taxes());
        for(size_t i = 0, e = node->first.size(); i < e; ++i) {
            if((val = node->first[i]) == 0) continue;
            for(unsigned char j = 0; j < '@'; ++j) { // @ is 64
#if !NDEBUG
                const size_t index((i << 6) + j);
                try {
                    if(node->first[i] & (1ul << j)) ks.putuw_(taxes.at(index)), ks.putc_(',');
                } catch (std::out_of_range &ex) {
                    LOG_EXIT("length of taxes: %zu. index: %zu.\n", taxes.size(), index);
                    throw;
                }
#else
                if(node->first[i] & (1ul << j)) ks.putuw_(taxes[(i << 6) + j]), ks.putc_(',');
#endif
            }
        }
        ks.back() = '\n';
        ks.terminate();
        // Maybe summary stats?
    }

    bool emplace_subtree(const tax_t parent, const std::uint32_t ntaxes) {
        if(ntaxes < 2) {
            LOG_DEBUG("Skipping subtree of one element. (parent: %u)\n", parent);
            return false;
        }
        subtrees_.emplace_back(parent, ntaxes, subtrees_.size());
        return true;
    }
public:
    // Also need a map of taxid to tax level.
    FMEmitter(khash_t(p) *tax, const std::unordered_map<tax_t, strlist> &taxpathmap): tax_{tax}, tpm_{taxpathmap} {}
    void run_collapse(tax_t maxtax, std::FILE* fp=stdout, std::size_t nelem=0) {
        ks::KString ks;
        run_collapse(maxtax, ks, fp, nelem);
    }
    void run_collapse(tax_t maxtax, ks::KString &ks, std::FILE* fp, std::size_t nelem) {
        assert(heap_.size());
        if(nelem == 0) {
            auto bccpy(kh_size(tax_));
            roundup64(bccpy);
            bccpy <<= 1;
            nelem = bccpy - kh_size(tax_);
            LOG_DEBUG("elements to add defaulting to %zu more [nodes addable before reaching 2 * nearest power of two (%zu).]\n", nelem, static_cast<std::size_t>(bccpy));
        }
        std::unordered_set<NodeType *> to_reinsert;
        ks.clear();
        ks.puts("#Taxid (inserted)\tScore\tParent\tChildren [comma-separated]\n");
        const int fd(fileno(fp));
        while(added_.size() < nelem) {
#if !NDEBUG
            ::emp::assert_sorted_impl<decltype(heap_), node_lt>(heap_);
            ::std::cerr << "Size of heap: " << heap_.size() << '\n';
            for(const auto node: heap_)
                ::std::cerr << node->second.str() << '\n';
#endif
            if(heap_.size() == 0) throw std::runtime_error("Heap is empty as collapse begins.");
            const auto bptr(*heap_.begin());
            const auto addn_score(get_score(*bptr));
            if(bptr->second.added()) {
                LOG_WARNING("Cannot add more nodes. [Best candidate is already added, which means we've inserted all possible.] Breaking from loop.\n");
                break;
            }
            bptr->second.laa_ = bptr;
            to_reinsert.insert(std::begin(bptr->second.parents_),
                               std::end(bptr->second.parents_));
            for(auto other: bptr->second.subsets_) {
                to_reinsert.insert(other);
                auto &desc(other->second);
                if((desc.laa_ == nullptr ||
                    desc.laa_->second.pc_ > bptr->second.pc_) && !desc.added())
                    desc.laa_ = bptr;
                // descendent's lca is set to bptr if its lca has fewer bits set and it has not been itself added to the map.
                // This check might be possible to remove later, as a descendent should always have fewer bits set.
                for(const auto parent: other->second.parents_)
                    if(!parent->second.added())
                        to_reinsert.insert(parent);
            }
            format_emitted_node(ks, bptr, addn_score, maxtax++);

            //  if(ks.size() >= (1 << 16))
            if(ks.size() & (1u << 16)) ks.write(fd), ks.clear();

            // Make a list of all pointers to remove and reinsert to the map.
            heap_.erase(heap_.begin());
            for(const auto el: to_reinsert) heap_.erase(el);
            heap_.insert(to_reinsert.begin(), to_reinsert.end()), to_reinsert.clear();
            heap_.insert(bptr);
        }
        ks.write(fd), ks.clear();
    }
    template<typename T>
    void process_subtree(const tax_t parent, T bit, T eit, const Spacer &sp, int num_threads=-1, khash_t(all) *acc=nullptr) {
        std::unordered_map<tax_t, strlist> tmpmap;
        typename std::unordered_map<tax_t, strlist>::const_iterator m;
        while(bit != eit) {
            if((m = tpm_.find(*bit++)) == tpm_.end()) LOG_DEBUG("No paths found for tax %u. Continuing.\n", *(bit - 1));
            else tmpmap.emplace(m->first, m->second);
        }
        if(emplace_subtree(parent, tmpmap.size())) {
            auto &ref(subtrees_.back());
            ref.fill(tmpmap, sp, num_threads, acc);
            ref.build_adjlist();
            ref.add_to_heap(heap_);
        }
    }
};

} // namespace emp

#endif // __FLEX_TREE_H
