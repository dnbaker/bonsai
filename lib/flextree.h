#ifndef __FLEX_TREE_H
#define __FLEX_TREE_H

#include "lib/bits.h"
#include "lib/counter.h"

namespace emp {

struct fnode_t {
    const std::vector<std::uint64_t> *bits_; // Does not own.
    const std::uint64_t               n_;    // Number of kmers at this point in tree.
    fnode_t                          *lua_;  // lowest unadded ancestor

    // Constructor
    fnode_t(const std::vector<std::uint64_t> &vec, const std::uint64_t n):
        bits_{&vec}, n_{n}, lua_{nullptr} {}

};

} // namespace emp

namespace std {

  template <>
  struct hash<emp::fnode_t>
  {
    const struct hash<vector<uint64_t>> vechash_;
    uint64_t operator()(const emp::fnode_t& node) const
    {
        return vechash_(*node.bits_);
    }
  };

}

namespace emp {


class fadjlist_t {

    // This structure might not work with insertion and reinsertion.
    // Will rethink and correct.

    std::unordered_map<const fnode_t *, std::vector<const fnode_t *>> map_;
    
    // Score generation
    std::uint64_t score(const fnode_t &node) const {
        const std::uint64_t bd(node.lua_ ? popcnt::vec_bitdiff(*node.bits_,
                                                  *node.lua_->bits_)
                                         : popcnt::vec_popcnt(*node.bits_));
        auto m(map_.find(&node));
        std::uint64_t n(node.n_);
        if(m != map_.end())
            for(const auto i: m->second)
                n += i->n_;
        return n * bd;
    }

    fadjlist_t() {
        throw std::runtime_error("NotImplementedError: Make constructor pls.");
    }
};

} // namespace emp

#endif // __FLEX_TREE_H
