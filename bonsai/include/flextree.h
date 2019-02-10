#ifndef __FLEX_TREE_H
#define __FLEX_TREE_H

#include "util.h"
#include "counter.h"
#include "bitmap.h"
#include "linear/linear.h"

namespace bns {

struct fnode_t;
using pop::vec_bitdiff;
using pop::vec_popcnt;
using NodeType = std::pair<const bitvec_t, fnode_t>;
using namespace std::literals;

struct fnode_t {
    u64                n_;  // Number of kmers at this point in tree.
    i64         desc_pts_;  // Change in score based on related node additions.
    const u32         pc_;  // popcount
    const u32         bc_;  // bitcount for family (number of clades in subtree)
    const u32      si_:31;  // subtree index
    u32          added_:1;  // Whether node has been added or been subsumed
    fnode_t(const bitvec_t &bits, u32 bc, u32 subtree_index, const u64 n=0):
        n_{n}, desc_pts_{0}, pc_{static_cast<u32>(pop::vec_popcnt(bits))},
        bc_{bc}, si_{subtree_index}, added_(false) {}
    ks::string str() const {
        return ks::sprintf("fnode_t[n:%zu,popcount:%u,familysize:%u", n_, pc_, bc_);
    }
    void subsume(NodeType &other) {
        const auto tmp((bc_ - pc_) * other.second.n_);
        desc_pts_ += tmp;
        other.second.desc_pts_ -= tmp;
    }
};

INLINE u64 get_score(const NodeType &node) {
    return (node.second.bc_ - node.second.pc_) * node.second.n_ + node.second.desc_pts_;
}
INLINE u64 get_score(const NodeType *node) {return get_score(*node);}

struct node_lt {
    bool operator()(const NodeType *a, const NodeType *b) const {
        return get_score(*a) < get_score(*b);
    }
};

class FlexMap {

    std::unordered_map<bitvec_t, fnode_t> map_;
    std::vector<tax_t>                    tax_;
    u64                                     n_;
    u32                              bitcount_;
    const u32                              id_;
    const tax_t                        parent_;

public:
    std::string str() const {
        ks::string ks;
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
    template<typename T, typename LT>
    void fill_heap(std::set<T, LT> &heap, size_t max_heap_size) const {
        for(auto &pair: map_) {
            std::pair<const bitvec_t, fnode_t> &ref(const_cast<std::pair<const bitvec_t, fnode_t> &>(pair));
            if(ref.second.added_) continue;
            //LOG_DEBUG("Adding pair where key's first entry is is %" PRIu64 " and the node is %s\n", ref.first[0], ref.second.str().data());
            //heap.insert(const_cast<std::pair<const std::vector<long long unsigned int>, bns::fnode_t>*>(&pair));
            if(heap.size() < max_heap_size) {
                heap.insert(&ref);
            } else if(get_score(pair) < get_score(*heap.begin())) {
                heap.erase(heap.begin());
                heap.insert(&ref);
            }
        }
    }
public:
    FlexMap(const tax_t parent, const u32 ntaxes, u32 id):
        n_{0}, bitcount_{ntaxes}, id_{id}, parent_{parent} {}

    INLINE void add(bitvec_t &&elem) {
#if __GNUC__ >= 7 && __cplusplus >= 201703L
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
    void fill(const std::unordered_map<tax_t, std::forward_list<std::string>> &list, const Spacer &sp, bool canonicalize=true, int num_threads=-1,
              khash_t(all) *acc=nullptr) {
        if(tax_.size()) tax_.clear();
        for(const auto &pair: list) {
            LOG_DEBUG("Filling from flexmap #%u with tax %u and is %s\n", id_, pair.first, pair.second.empty() ? "empty": "not bnsty");
            tax_.push_back(pair.first);
        }
#if !NDEBUG
        bitmap_t tmp_bitmap(kgset_t(list, sp, canonicalize, num_threads, acc));
        auto &map(tmp_bitmap.get_map());
        LOG_DEBUG("Map size: %zu\n", map.size());
        for(auto &&pair: map) add(std::move(pair.second));
#else
        for(auto &&pair: bitmap_t(kgset_t(list, sp, canonicalize, num_threads, acc)).get_map()) add(std::move(pair.second));
#endif
    }
};

class FMEmitter {
    using HeapType = std::set<NodeType *, node_lt>;
    using TaxPathType = std::unordered_map<tax_t, std::forward_list<std::string>>;
    static const size_t BufferSize = 1 << 17;

    std::vector<FlexMap>   subtrees_;
    HeapType                   heap_;
    khash_t(p)           *const tax_;
    const TaxPathType          &tpm_;
    const size_t         max_heapsz_;
    unsigned            left_to_add_;
    const bool                canon_;


    /*
     * Emits an additional node to the tree where its paren
    */
    void format_emitted_node(ks::string &ks, const NodeType *node, const tax_t taxid) const {
        u64 val;
        const auto &fm(subtrees_[node->second.si_]);
        ks.putuw_(taxid);
        ks.putc_('\t');
        ks.putl_(get_score(*node));
        ks.putc_('\t');
        ks.putuw_(fm.parent());
        const auto &taxes(fm.get_taxes());
        for(size_t i = 0, e = node->first.size(); i < e; ++i) {
            if((val = node->first[i]) == 0) continue;
            for(unsigned j = 0; j < 64; ++j) { // @ is 64
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

    bool emplace_subtree(const tax_t parent, const u32 ntaxes) {
        if(ntaxes < 2) {
            LOG_DEBUG("Skipping subtree of one element. (parent: %u)\n", parent);
            return false;
        }
        subtrees_.emplace_back(parent, ntaxes, subtrees_.size());
        return true;
    }
public:
    // Also need a map of taxid to tax level.
    FMEmitter(khash_t(p) *tax, const std::unordered_map<tax_t, std::forward_list<std::string>> &taxpathmap, bool canonicalize=true,
              size_t max_heapsz=1<<8, size_t to_add=0):
        tax_{tax}, tpm_{taxpathmap}, max_heapsz_{max_heapsz},
        left_to_add_((to_add ? to_add
                             : roundup64(kh_size(tax_))) - kh_size(tax_)),
        canon_{canonicalize}
    {
    }
    void run_collapse(tax_t maxtax, std::FILE* fp=stdout) {
        ks::string ks;
        run_collapse(maxtax, ks, fp);
    }
    struct stlt {
        bool operator()(const NodeType *a, const NodeType *b) const {
            return a->second.pc_ > b->second.pc_;
        }
        // Places the highest bitcount items at the beginning.
    };
    void condense_subtree(HeapType &subtree) {
        using std::begin;
        using std::end;
        // Update scores
        // Mark
        using HeapIt = typename HeapType::iterator;
        HeapIt ait, bit, eait;
        using std::begin;
        using std::end;
        using std::find;
        //std::vector<linear::set<NodeType *>> subtree_subsets;
        linear::set<NodeType *> to_remove;
#if 0
        auto in = []<class T>(std::vector<T *> &vec, T *el) {
            return find(begin(vec), end(vec), el) != vec.end();
        };
        auto insert = [&]<class T>(std::vector<T *> &vec, bitvec_t *el) {
            if(!in(vec, el)) vec.emplace_back(el);
        };
#endif
        // 1. Make linear vector of pointers to objects.
        // 2. 
        for(ait = begin(subtree), eait = end(subtree); ait != eait; ++ait) {
            for(bit = ait, ++bit; bit != eait; ++bit) {
                if((*ait)->second.added_|| (*bit)->second.added_) continue;
                //if(&(*ait)->first) continue;
                switch(veccmp((*ait)->first, (*bit)->first)) {
                    case BitCmp::FIRST_PARENT:
                        (*ait)->second.subsume(**bit);
                        to_remove.insert(*ait);
                        to_remove.insert(*bit);
                        break;
                    case BitCmp::SECOND_PARENT:
                        (*bit)->second.subsume(**ait);
                        to_remove.insert(*ait);
                        to_remove.insert(*bit);
                        break;
                    case BitCmp::EQUAL:
                        LOG_DEBUG("Warning: bit patterns a and b should be not equal....\n");
                        [[fallthrough]];
                    case BitCmp::INCOMPARABLE: continue;//Do nothing
                }
            }
        }
        // I think I should instead think of some other way to organize them, but this works for a start.
        for(auto el: to_remove) subtree.erase(el);
    }
    void run_collapse(tax_t maxtax, ks::string &ks, std::FILE* fp) {
        const int fd(fileno(fp));
        ks.clear();
        ks.puts("#Taxid (inserted)\tScore\tParent\tChildren [comma-separated]\n");
        ks.write(fd), ks.clear();
#if COMPLICATED_BUT_WEIRD
        for(;;) {
            while(heap_.size() < max_heapsz_)  {
                // This could loop forever if all nodes are added.
                // I don't expect this to happen, but it's worth noting.
                for(auto &subtree: subtrees_) {
                    std::set<NodeType *, node_lt> tmptree;
                    subtree.fill_heap(tmptree, max_heapsz_);
                    condense_subtree(tmptree);
                    heap_.merge(tmptree);
                }
            }
            for(auto it(heap_.crbegin()), sit(heap_.crend()); it != sit;) {
                format_emitted_node(ks, *it++, ++maxtax);
                (*it)->second.added_ = 1;
                if(ks.size() & ~(BufferSize-1)) {
                    ks.write(fd), ks.clear();
                }
                if(--left_to_add_ == 0) return;
            }
            std::fputs("I could continue this loop to add more nodes, but I am exiting for the proof of concept.\n", stderr);
            break;
        }
#else
        while(left_to_add_) {
            for(auto &subtree: subtrees_) {
                subtree.fill_heap(heap_, max_heapsz_);
                // condense_subtree(heap_);
            }
            for(auto it(heap_.crbegin()), sit(heap_.crend()); it != sit;) {
                format_emitted_node(ks, *it++, ++maxtax);
                (*it)->second.added_ = 1;
                if(ks.size() & ~(BufferSize-1)) {
                    ks.write(fd), ks.clear();
                }
                if(--left_to_add_ == 0) break;
            }
        }
        ks.write(fd), ks.clear();
        return;
#endif
    }
    template<typename T>
    void process_subtree(const tax_t parent, T bit, T eit, const Spacer &sp, int num_threads=-1, khash_t(all) *acc=nullptr) {
        std::unordered_map<tax_t, std::forward_list<std::string>> tmpmap;
        typename std::unordered_map<tax_t, std::forward_list<std::string>>::const_iterator m;
        while(bit != eit) {
            if((m = tpm_.find(*bit++)) == tpm_.end()) LOG_DEBUG("No paths found for tax %u. Continuing.\n", *(bit - 1));
            else tmpmap.emplace(m->first, m->second);
        }
        if(emplace_subtree(parent, tmpmap.size())) {
            auto &ref(subtrees_.back());
            ref.fill(tmpmap, sp, canon_, num_threads, acc);
        }
    }
};

} // namespace bns

#endif // __FLEX_TREE_H
