#ifndef __FLEX_TREE_H
#define __FLEX_TREE_H

#include "lib/bits.h"
#include "lib/counter.h"

namespace emp {

struct fnode_t;
using node_type = std::pair<std::vector<std::uint64_t>, fnode_t>;

struct fnode_t {
    std::uint64_t                  n_;  // Number of kmers at this point in tree.
    node_type                   *lua_;  // lowest unadded ancestor
    std::vector<node_type *> subsets_;

    // Constructor
    fnode_t(const std::uint64_t n=0):
        n_{n}, lua_{nullptr} {}


};

INLINE std::uint64_t get_score(const node_type &node) {
    std::uint64_t ret(node.second.n_);
    for(auto s: node.second.subsets_) {
        if(s->second.lua_ == nullptr) ret += s->second.n_;
        else if(popcnt::vec_bitdiff(s->first, node.first) <
                popcnt::vec_bitdiff(s->first, s->second.lua_->first)) {
            ret += s->second.n_;
            s->second.lua_ = const_cast<node_type*>(&node);
        }
    }
    return ret;
}

struct node_lt {
    bool operator()(const node_type *a, const node_type *b) const {
        return get_score(*a) > get_score(*b);
    }
};


class FlexMap {
    std::unordered_map<std::vector<std::uint64_t>, fnode_t>             map_;
    std::set<std::pair<std::vector<std::uint64_t>, fnode_t> *, node_lt> heap_;
};

} // namespace emp

#endif // __FLEX_TREE_H
