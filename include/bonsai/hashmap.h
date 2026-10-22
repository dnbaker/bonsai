#ifndef BONSAI_HASHMAP_H__
#define BONSAI_HASHMAP_H__
// Drop-in replacements for std::unordered_map / std::unordered_set, per the
// PI 2026-10-21 directive to remove the std unordered containers from the
// first-party tree. ankerl::unordered_dense ships in the sketch submodule
// (hll/include/ankerl, MIT; already a build dependency via sketch/bmh.h) and
// beats std::unordered_map on both speed and memory (open addressing,
// wyhash for strings and identity-mix for integers).
//
// umap/uset: flat tables - fastest, but references/pointers to elements are
// invalidated by rehashing, exactly like the std containers do NOT do (the
// std unordered containers are node-based). Use smap where element
// addresses must stay stable (e.g. per-taxon tallies whose references are
// cached by threads); segmented_map never moves elements.
#include "ankerl/unordered_dense.h"

namespace bns {

template<typename K, typename V, typename H = ankerl::unordered_dense::hash<K>, typename E = std::equal_to<K>>
using umap = ankerl::unordered_dense::map<K, V, H, E>;

template<typename K, typename V, typename H = ankerl::unordered_dense::hash<K>, typename E = std::equal_to<K>>
using smap = ankerl::unordered_dense::segmented_map<K, V, H, E>;

template<typename K, typename H = ankerl::unordered_dense::hash<K>, typename E = std::equal_to<K>>
using uset = ankerl::unordered_dense::set<K, H, E>;

} // namespace bns

#endif /* BONSAI_HASHMAP_H__ */
