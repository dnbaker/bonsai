/*
Grafting optimization for the synthetic-node placement technique
(doc/grafting-objective.md): choose <= B grafts (new synthetic nodes) to
minimize the total Hamming distance between the tree's effective k-mer
bitmaps and the true exhaustive bitmaps.

Operative model: the tree's objective evaluation (placement.h implied/err)
places each k-mer class at the deepest node containing it, so grafts can
steal classes from other grafts — the best-fit model. Gains are therefore
not necessarily submodular (see the counterexample test in
test/grafting.cpp); greedy is the primary algorithm and takes only
positive-gain steps, so H_total is non-increasing per graft.

Exactness note: gain evaluation is incremental but exact — a graft only
moves classes whose true bitmap intersects the graft set (stealing removes
only graft-set leaves from shrunk ancestors, so masks and implied nodes of
all other classes are unchanged), and those classes are re-evaluated on a
scratch copy of the tree. No per-genome k-mer bitmaps are materialized
beyond the class counts already in KmerClasses.
*/

#ifndef __BNS_GRAFTING_H
#define __BNS_GRAFTING_H

#include "placement.h"
#include "sketch_select.h"

#include <algorithm>
#include <numeric>
#include <optional>
#include <vector>

namespace bns {

// Class index for incremental gain evaluation: multi-genome classes with
// per-genome inverted lists.
struct GraftIndex {
    std::vector<u64> masks;                       // class bitmaps
    std::vector<u64> counts;                      // distinct k-mers per class
    std::vector<std::vector<u32>> by_genome;      // genome -> class indices
    u64 initial_delta;
};

inline GraftIndex build_graft_index(const PlacementTree &tree, const KmerClasses &kc) {
    GraftIndex gi;
    gi.by_genome.resize(kc.ngenomes);
    for(const auto &pair: kc.count) {
        if(pc64(pair.first) < 2) continue;
        const u32 id(u32(gi.masks.size()));
        gi.masks.push_back(pair.first);
        gi.counts.push_back(pair.second);
        for(u64 m(pair.first); m; m &= m - 1)
            gi.by_genome[__builtin_ctzll(m)].push_back(id);
    }
    gi.initial_delta = tree.total_delta(kc);
    return gi;
}

/*
Exact gain (delta reduction) of adding graft A to `tree`, without mutating
it. Only classes whose bitmap intersects A can change implied masks:
add_synthetic removes only A-leaves from existing masks, so a class
disjoint from A keeps every containing node and hence its implied mask.
Returns 0 for infeasible A (singleton or already a node).

before_cache, when provided, memoizes tree.err(mask) per class id; it is
only valid while `tree` is unchanged, so callers must reset it (fill with
u64(-1)) after mutating the tree.
*/
inline long double gain_of_add(const PlacementTree &tree, const GraftIndex &gi, u64 A,
                               std::vector<u64> *before_cache = nullptr) {
    if(pc64(A) < 2) return 0.0L;
    // affected classes via per-genome inverted lists
    std::vector<u32> affected;
    std::vector<u32> seen(gi.masks.size(), u32(-1));
    for(u64 m(A); m; m &= m - 1) {
        for(const u32 id: gi.by_genome[__builtin_ctzll(m)]) {
            if(seen[id] == u32(-1)) {seen[id] = 1; affected.push_back(id);}
        }
    }
    PlacementTree tmp(tree);
    if(tmp.add_synthetic(A) < 0) return 0.0L;
    long double gain(0.0L);
    for(const u32 id: affected) {
        const u64 S(gi.masks[id]);
        u64 before;
        if(before_cache) {
            u64 &cached((*before_cache)[id]);
            if(cached == u64(-1)) cached = tree.err(S);
            before = cached;
        } else {
            before = tree.err(S);
        }
        const u64 after(tmp.err(S));
        gain += (long double)(gi.counts[id]) * ((long double)(before) - (long double)(after));
    }
    return gain;
}

// Apply the best candidate greedily for up to budget steps; accepts a
// candidate only on strictly positive exact gain (H never increases).
struct GraftingResult {
    std::vector<u64> added_masks;
    std::vector<u64> deltas;      // total_delta after each accepted graft
    std::vector<u64> fn, fp;      // Hamming decomposition after each step
    u64 initial_delta;
};

struct HammingDecomp { u64 fn, fp; };

inline HammingDecomp hamming_decomp(const PlacementTree &tree, const KmerClasses &kc) {
    HammingDecomp ret{0, 0};
    for(const auto &pair: kc.count) {
        if(pc64(pair.first) < 2) continue;
        const u64 implied(tree.implied(pair.first));
        ret.fn += pair.second * pc64(pair.first & ~implied);
        ret.fp += pair.second * pc64(implied & ~pair.first);
    }
    return ret;
}

inline GraftingResult graft_greedy(PlacementTree &tree, const KmerClasses &kc,
                                   const std::vector<u64> &candidates, u32 budget) {
    GraftingResult ret;
    GraftIndex gi(build_graft_index(tree, kc));
    ret.initial_delta = gi.initial_delta;
    ret.deltas.push_back(gi.initial_delta);
    const HammingDecomp d0(hamming_decomp(tree, kc));
    ret.fn.push_back(d0.fn); ret.fp.push_back(d0.fp);
    std::vector<char> used(candidates.size(), 0);
    // tree.err(mask) per class, valid for one greedy round at a time.
    std::vector<u64> before_cache(gi.masks.size(), u64(-1));
    while(ret.added_masks.size() < budget) {
        long double best_gain(0.0L);
        u32 best(u32(-1));
        for(u32 i(0); i < candidates.size(); ++i) {
            if(used[i]) continue;
            const long double g(gain_of_add(tree, gi, candidates[i], &before_cache));
            if(g > best_gain) {best_gain = g; best = i;}
        }
        if(best == u32(-1)) break;
        used[best] = 1;
        const u64 A(candidates[best]);
        if(tree.add_synthetic(A) < 0) continue; // raced by an earlier graft; skip
        std::fill(before_cache.begin(), before_cache.end(), u64(-1));
        ret.added_masks.push_back(A);
        const u64 d(tree.total_delta(kc));
        ret.deltas.push_back(d);
        const HammingDecomp dd(hamming_decomp(tree, kc));
        ret.fn.push_back(dd.fn); ret.fp.push_back(dd.fp);
    }
    return ret;
}

/*
Feasible grafts (keep the family laminar, doc/grafting-objective.md):
- union grafts: every non-empty proper union of children of an internal node;
- split grafts: clade(v) minus one child clade, for every internal v.
Both are produced in mask form. Existing node masks, singletons, and the
full leaf set are excluded. With all_subsets (only for tiny leaf counts),
every non-empty proper subset that is not an existing node mask is added so
the exact reference explores the complete space.
*/
inline std::vector<u64> enumerate_grafts(const PlacementTree &tree, bool all_subsets=false) {
    std::vector<u64> node_masks;
    for(u32 id(0); id < tree.nnodes(); ++id)
        if(tree.is_alive(id) && !tree.is_leaf_node(id)) node_masks.push_back(tree.mask(id));
    auto exists = [&](u64 A) {
        return std::find(node_masks.begin(), node_masks.end(), A) != node_masks.end();
    };
    std::vector<u64> ret;
    auto push = [&](u64 A) {
        if(pc64(A) < 2) return;
        if(std::find(ret.begin(), ret.end(), A) != ret.end()) return;
        if(exists(A)) return;
        ret.push_back(A);
    };
    for(u32 id(0); id < tree.nnodes(); ++id) {
        if(!tree.is_alive(id) || tree.is_leaf_node(id)) continue;
        const auto &ch(tree.children_of(id));
        if(ch.empty()) continue;
        const u64 full(tree.mask(id));
        // union grafts: all non-empty proper subsets of children
        const size_t k(ch.size());
        if(k <= 20) {
            for(u64 s(1); s < (u64(1) << k) - 1; ++s) {
                u64 A(0);
                for(size_t i(0); i < k; ++i) if(s & (u64(1) << i)) A |= tree.mask(ch[i]);
                push(A);
            }
        }
        // split grafts: clade minus one child clade
        for(size_t i(0); i < k; ++i) push(full & ~tree.mask(ch[i]));
    }
    if(all_subsets) {
        const u32 n(tree.ngenomes());
        const u64 all(tree.genome_mask());
        for(u64 s(3); s < all; ++s) if((s & ~all) == 0) push(s);
    }
    return ret;
}

/*
Exact reference: best set of <= budget grafts under the operative
(best-fit) dynamics, by branch-and-bound over candidate subsets in index
order (no lookahead bound: gains are not submodular, so a single-gain bound
is not admissible). Suitable for small candidate sets and budgets.
*/
struct ExactGraftResult {
    std::vector<u64> grafts;
    u64 final_delta;
};

inline ExactGraftResult graft_exact(const PlacementTree &tree0, const KmerClasses &kc,
                                    const std::vector<u64> &candidates, u32 budget) {
    ExactGraftResult best; best.final_delta = tree0.total_delta(kc);
    std::vector<u64> chosen;
    std::function<void(const PlacementTree &, u32, u32, u64)> dfs = [&](const PlacementTree &tree, u32 start, u32 depth, u64 delta) {
        if(delta < best.final_delta) {
            best.final_delta = delta;
            best.grafts = chosen;
        }
        if(depth == budget) return;
        for(u32 i(start); i < candidates.size(); ++i) {
            PlacementTree tmp(tree);
            if(tmp.add_synthetic(candidates[i]) < 0) continue;
            chosen.push_back(candidates[i]);
            dfs(tmp, i + 1, depth + 1, tmp.total_delta(kc));
            chosen.pop_back();
        }
    };
    dfs(tree0, 0, 0, best.final_delta);
    return best;
}

// Sketch-side ranking signal for a graft A: estimated shared distinct
// k-mers between the graft's clade sketch and its parent-position sketch
// (register operations only). Used to order candidates before exact
// verification; never used as an objective value.
inline double sketch_shared_estimate(const PlacementTree &tree,
                                       const std::vector<GenomeSketch> &leaf_sketches,
                                       u64 A) {
    const u32 n(tree.ngenomes());
    GenomeSketch clade(leaf_sketches[__builtin_ctzll(A)].getp(), GENOME_SKETCH_B);
    GenomeSketch parent(leaf_sketches[__builtin_ctzll(A)].getp(), GENOME_SKETCH_B);
    const u32 lca_node(tree.lca(A));
    const u64 parent_mask(lca_node == tree.root() ? tree.genome_mask(): tree.mask(tree.parent_of(lca_node)));
    for(u32 g(0); g < n; ++g) {
        if(A & (u64(1) << g)) clade += leaf_sketches[g];
        if(parent_mask & (u64(1) << g)) parent += leaf_sketches[g];
    }
    auto fc(clade.finalize()), fp(parent.finalize());
    return fc.intersection_size(fp);
}

} // namespace bns

#endif // __BNS_GRAFTING_H
