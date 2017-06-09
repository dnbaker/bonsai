#ifndef __FLEX_TREE_H
#define __FLEX_TREE_H

#include "lib/bits.h"
#include "lib/counter.h"

namespace emp {

struct fnode_t;
using NodeType = std::pair<const std::vector<std::uint64_t>, fnode_t>;
using popcnt::vec_bitdiff;
using popcnt::vec_popcnt;

template<typename T> class TD;

struct fnode_t {
    std::uint64_t                        n_;  // Number of kmers at this point in tree.
    const NodeType                    *lua_;  // lowest unadded ancestor
    std::vector<NodeType *>        subsets_;
    const std::uint32_t                 bc_;

    fnode_t(std::uint32_t bc, const std::uint64_t n):
        n_{n}, lua_{nullptr}, bc_{bc} {}

    bool added() const {return lua_ && &lua_->second == this;}
};


INLINE std::uint64_t get_score(const NodeType &node) {
    std::uint64_t ret(node.second.n_);
    for(auto s: node.second.subsets_) {
        if(s->second.added()) continue;
        if(s->second.lua_ == nullptr) ret += s->second.n_;
        else if(s->second.lua_ != s &&
                vec_bitdiff(s->first, node.first) <
                    vec_bitdiff(s->first, s->second.lua_->first)) {
            ret += s->second.n_;
            s->second.lua_ = &node;
        }
    }
    return node.second.lua_ ? vec_bitdiff(node.second.lua_->first, node.first) * ret
                            : node.second.bc_ - vec_popcnt(node.first) * ret;
}


class FlexMap {

    struct node_lt {
        bool operator()(const NodeType *a, const NodeType *b) const {
            return get_score(*a) > get_score(*b);
        }
    };

    std::unordered_map<std::vector<std::uint64_t>, fnode_t> map_;
    std::set<NodeType *, node_lt>                           heap_;
    std::uint64_t                                           n_;
    std::uint32_t                                           bitcount_;

public:
    void add(std::vector<std::uint64_t> &&elem) {
#if __GNUC__ >= 7
        if(auto match = map_.find(elem); match == map_.end())
#else
        auto match(map_.find(elem));
        if(match == map_.end())
#endif
            map_.emplace(std::move(elem),
                         bitcount_, UINT64_C(1));
        else ++match->second.n_;
        ++n_;
    }
    void fill_heap() {
        for(auto &pair: map_) {
            heap_.insert(&pair);
        }
    }
    FlexMap(std::uint32_t bitcount): n_{0}, bitcount_{bitcount} {
    }
};

} // namespace emp

#endif // __FLEX_TREE_H
