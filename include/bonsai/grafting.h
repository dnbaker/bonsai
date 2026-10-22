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
#include <chrono>
#include <numeric>
#include <optional>
#include <vector>

#ifdef _OPENMP
#include <omp.h>
#endif

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
u64(-1)) after mutating the tree. seen_scratch/generation stamp the
affected-class dedup set so a full O(nclasses) fill is not paid per
candidate; pass nullptr/0 to allocate one internally.
*/
inline long double gain_of_add(const PlacementTree &tree, const GraftIndex &gi, u64 A,
                               std::vector<u64> *before_cache = nullptr,
                               std::vector<u32> *seen_scratch = nullptr, u32 generation = 0,
                               std::vector<u32> *affected_scratch = nullptr) {
    if(pc64(A) < 2) return 0.0L;
    // affected classes via per-genome inverted lists
    std::vector<u32> internal_affected;
    std::vector<u32> &affected(affected_scratch ? *affected_scratch: internal_affected);
    affected.clear();
    std::vector<u32> internal_seen;
    std::vector<u32> &seen(seen_scratch ? *seen_scratch: internal_seen);
    if(seen_scratch == nullptr) seen.assign(gi.masks.size(), u32(-1));
    for(u64 m(A); m; m &= m - 1) {
        for(const u32 id: gi.by_genome[__builtin_ctzll(m)]) {
            const bool fresh(seen_scratch ? seen[id] != generation: seen[id] == u32(-1));
            if(fresh) {
                if(seen_scratch) seen[id] = generation;
                else seen[id] = 1;
                affected.push_back(id);
            }
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

/*
Admissible upper bound on gain_of_add for the same tree: the gain is
sum(count * (before - after)) and after >= 0, so the gain cannot exceed
sum(count * before) over the same affected set. Costs only the inverted-
list scan (no tree copy, no add_synthetic), which lets the greedy skip
exact evaluation of candidates that cannot beat the current best. Fills the
same before_cache; valid only while `tree` is unchanged.
*/
inline long double gain_bound(const PlacementTree &tree, const GraftIndex &gi, u64 A,
                              std::vector<u64> &before_cache,
                              std::vector<u32> &seen, u32 generation) {
    if(pc64(A) < 2) return 0.0L;
    long double bound(0.0L);
    for(u64 m(A); m; m &= m - 1) {
        for(const u32 id: gi.by_genome[__builtin_ctzll(m)]) {
            if(seen[id] == generation) continue;
            seen[id] = generation;
            u64 &cached(before_cache[id]);
            if(cached == u64(-1)) cached = tree.err(gi.masks[id]);
            bound += (long double)(gi.counts[id]) * (long double)(cached);
        }
    }
    return bound;
}

// Apply the best candidate greedily for up to budget steps; accepts a
// candidate only on strictly positive exact gain (H never increases).
struct GraftingResult {
    std::vector<u64> added_masks;
    std::vector<u64> deltas;      // total_delta after each accepted graft
    std::vector<u64> fn, fp;      // Hamming decomposition after each step
    u64 initial_delta;
    // Measurement only (no algorithmic role): per-round candidate-scan
    // statistics, so the cost structure of the greedy is visible at scale.
    u64 scanned = 0;              // candidates considered (not previously used)
    u64 bound_skipped = 0;        // rejected by the admissible gain_bound
    u64 exact_evals = 0;          // paid tree-copy + add_synthetic evaluations
    std::vector<double> round_seconds;
    std::vector<double> est_gains;   // sketch-only greedy: the estimate behind each accepted graft
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
                                   const std::vector<u64> &candidates, u32 budget, double min_gain = 0.) {
    GraftingResult ret;
    GraftIndex gi(build_graft_index(tree, kc));
    ret.initial_delta = gi.initial_delta;
    ret.deltas.push_back(gi.initial_delta);
    const HammingDecomp d0(hamming_decomp(tree, kc));
    ret.fn.push_back(d0.fn); ret.fp.push_back(d0.fp);
    std::vector<char> used(candidates.size(), 0);
    // tree.err(mask) per class, valid for one greedy round at a time.
    std::vector<u64> before_cache(gi.masks.size(), u64(-1));
    // The per-round candidate scan is embarrassingly parallel: every
    // candidate is evaluated against the same unmutated tree (the single
    // accepted graft is applied after the join), so threads only need
    // per-thread gain slots and affected-set stamps; before_cache is filled
    // eagerly at round start and read-only during the scan. The admissible
    // gain_bound filters candidates that cannot beat the running best
    // before paying for a tree copy + add_synthetic + after-err walks.
    u32 nth(1);
#ifdef _OPENMP
    nth = u32(std::min(omp_get_max_threads(), std::max(1, int(candidates.size()) / 64)));
#endif
    std::vector<std::vector<u32>> tseen(nth);
    for(auto &v: tseen) v.assign(gi.masks.size(), 0);
    std::vector<std::vector<u32>> taffected(nth);
    std::vector<long double> tgain(nth, 0.0L);
    std::vector<u32> tbest(nth, u32(-1));
    std::vector<u32> tgen(nth, 0);
    while(ret.added_masks.size() < budget) {
        // Eagerly memoize every class's current err for this round.
        const auto round_t0(std::chrono::steady_clock::now());
        for(u32 id(0); id < gi.masks.size(); ++id)
            before_cache[id] = tree.err(gi.masks[id]);
        std::fill(tgain.begin(), tgain.end(), 0.0L);
        std::fill(tbest.begin(), tbest.end(), u32(-1));
#ifdef _OPENMP
        #pragma omp parallel for schedule(dynamic, 8) num_threads(int(nth))
#endif
        for(u32 i = 0; i < candidates.size(); ++i) {
            if(used[i]) continue;
#ifdef _OPENMP
            const u32 tid = u32(omp_get_thread_num());
#else
            const u32 tid(0);
#endif
            std::vector<u32> &seen(tseen[tid]);
            u32 &gen(tgen[tid]);
            if(++gen == 0) { // wrapped: re-zero the stamp array
                std::fill(seen.begin(), seen.end(), 0);
                gen = 1;
            }
#ifdef _OPENMP
            #pragma omp atomic
            ++ret.scanned;
#else
            ++ret.scanned;
#endif
            if(gain_bound(tree, gi, candidates[i], before_cache, seen, gen) <= tgain[tid]) {
#ifdef _OPENMP
                #pragma omp atomic
                ++ret.bound_skipped;
#else
                ++ret.bound_skipped;
#endif
                continue;
            }
            if(++gen == 0) {
                std::fill(seen.begin(), seen.end(), 0);
                gen = 1;
            }
            const long double g(gain_of_add(tree, gi, candidates[i], &before_cache, &seen, gen, &taffected[tid]));
#ifdef _OPENMP
            #pragma omp atomic
            ++ret.exact_evals;
#else
            ++ret.exact_evals;
#endif
            if(g > tgain[tid]) {tgain[tid] = g; tbest[tid] = i;}
        }
        long double best_gain(min_gain);   // acceptance threshold (round-5 parity knob; 0 = strictly-positive exact gains)
        u32 best(u32(-1));
        for(u32 t(0); t < nth; ++t)
            if(tbest[t] != u32(-1) && tgain[t] > best_gain) {best_gain = tgain[t]; best = tbest[t];}
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
        ret.round_seconds.push_back(std::chrono::duration<double>(std::chrono::steady_clock::now() - round_t0).count());
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

At fanout > EXACT_UNION_MAX (20), enumerating all 2^k - 2 child unions is
infeasible and used to drop the union lane entirely - the measured quality
cliff of doc/bonsai-gtdb-benchmark.md (stage F: fanout 24/40 recover ~8%
of initial delta with splits alone). Instead enumerate the polynomial
lane: pairwise and triplet child unions, O(k^2 + k^3), which compose
hierarchically (an accepted graft becomes a child of the same node, so
later rounds can union it with more children). Targeted unions of larger
arity come from the observed-class lane (class masks are candidates
themselves) and from enumerate_seeded_grafts below.
*/
inline constexpr unsigned EXACT_UNION_MAX = 20;

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
        // union grafts: all non-empty proper subsets of children at low
        // fanout; pairwise + triplet subsets above the cap
        const size_t k(ch.size());
        if(k <= EXACT_UNION_MAX) {
            for(u64 s(1); s < (u64(1) << k) - 1; ++s) {
                u64 A(0);
                for(size_t i(0); i < k; ++i) if(s & (u64(1) << i)) A |= tree.mask(ch[i]);
                push(A);
            }
        } else {
            for(size_t i(0); i < k; ++i)
                for(size_t j(i + 1); j < k; ++j)
                    push(tree.mask(ch[i]) | tree.mask(ch[j]));
            for(size_t i(0); i < k; ++i)
                for(size_t j(i + 1); j < k; ++j)
                    for(size_t l(j + 1); l < k; ++l)
                        push(tree.mask(ch[i]) | tree.mask(ch[j]) | tree.mask(ch[l]));
        }
        // split grafts: clade minus one child clade
        for(size_t i(0); i < k; ++i) push(full & ~tree.mask(ch[i]));
    }
    if(all_subsets) {
        const u64 all(tree.genome_mask());
        for(u64 s(3); s < all; ++s) if((s & ~all) == 0) push(s);
    }
    return ret;
}

/*
Observed-class-seeded union hulls: for each seed mask S (typically an
observed k-mer class) and each internal node v, the graft is the union of
the children of v whose clades intersect S - the smallest laminar
child-union covering S's footprint at v. O(|S| + fanout) per (S, v), so it
names high-arity unions directly where the poly lane (pairs/triplets)
would need many rounds to compose. Duplicates and existing node masks are
dropped; |seeds| x nnodes grafts at most.
*/
inline std::vector<u64> enumerate_seeded_grafts(const PlacementTree &tree, const std::vector<u64> &seeds) {
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
    for(const u64 S: seeds) {
        if(!S) continue;
        for(u32 id(0); id < tree.nnodes(); ++id) {
            if(!tree.is_alive(id) || tree.is_leaf_node(id)) continue;
            const auto &ch(tree.children_of(id));
            if(ch.size() < 2) continue;
            u64 A(0);
            for(const u32 c: ch) {
                const u64 cm(tree.mask(c));
                if(S & cm) A |= cm;
            }
            push(A);
        }
    }
    return ret;
}

/*
Round 4 (PI directive 2026-10-21): sketch-only naming. Candidate generation
and gain estimation use per-genome bbminhash INTERSECTION CARDINALITIES
alone - no exact k-mer sets, no inverted-list index, at any point in the
path. This is the scaling law: huge databases afford sketches, not indexes.
Design note and failure modes: doc/sketch-only-naming.md.
*/

// BBitMinHasher keeps core_ protected; the sketch-only lane needs
// register-level construction (intersection composites) without re-adding
// keys through the hashing add().
struct RegisterSketch: GenomeSketch {
    using GenomeSketch::GenomeSketch;
    auto &mutable_core() {return this->core_;}
};

// Register-wise union sketch over a genome mask (clade sketch).
template<typename SketchT>
inline RegisterSketch clade_sketch(const std::vector<SketchT> &sks, u64 mask) {
    const u32 g0(u32(__builtin_ctzll(mask)));
    RegisterSketch ret(sks[g0].getp(), GENOME_SKETCH_B);
    auto &dst(ret.mutable_core());
    dst = sks[g0].core();
    for(u64 m(mask & (mask - 1)); m; m &= m - 1) {
        const auto &src(sks[__builtin_ctzll(m)].core());
        for(size_t i(0); i < dst.size(); ++i) dst[i] = std::min(dst[i], src[i]);
    }
    return ret;
}

// Pairwise intersection cardinality, nshared * m / tmp with the submodule's
// harmonic tmp BUT with registers empty in BOTH sketches excluded from the
// equality count: the vendored BBitMinHasher::intersection_size counts
// both-empty registers as matches (its own jaccard_index excludes them),
// which biases intersections up by both_empty * m / tmp - 1847 for two
// disjoint 10k-key sketches at p = 13 (measured, doc/sketch-only-naming.md).
inline double pair_intersection_count(const GenomeSketch &x, const GenomeSketch &y) {
    const auto &a(x.core());
    const auto &b(y.core());
    const size_t m(a.size());
    const uint64_t EMPTY(std::numeric_limits<uint64_t>::max());
    const long double numinv(1. / std::ldexp((long double)1., 64 - x.getp()));
    long double tmp(0);
    size_t nshared(0);
    for(size_t i(0); i < m; ++i) {
        const uint64_t va(a[i]), vb(b[i]);
        const uint64_t v(std::min(va, vb));
        tmp = std::fma(numinv, v != EMPTY ? (long double)v: (long double)(EMPTY >> x.getp()), tmp);
        nshared += (va != EMPTY || vb != EMPTY) && va == vb;
    }
    return nshared * double(m) / double(tmp);
}

// k-way intersection cardinality from one-permutation registers, by the
// k-th order statistic (max of the member minima; Li-Owen-Zhang "k-min"
// form, equal-cardinality approximation): with V the register-wise
// maximum and f = V / 2^(64-p) the hash fraction,
//     E[2^-f] ~= k! (n / (m ln 2))^k   =>   n_hat = m ln2 (mean(2^-f)/k!)^(1/k)
// A bucket whose maximum is EMPTY (any member has no key there)
// contributes 0 - intersection keys cannot exist in that bucket. For k =
// 1 this reduces to the classic min-hash cardinality estimator.
inline double kway_intersection_count(const std::vector<GenomeSketch> &sks, u64 mask, u64 extra = 0) {
    const u32 k(pc64(mask) + pc64(extra));
    if(!k) return 0.;
    const size_t m(sks[0].core().size());
    const u32 p(sks[0].getp());
    const uint64_t EMPTY(std::numeric_limits<uint64_t>::max());
    RegisterSketch ctx_sk(sks[0].getp(), GENOME_SKETCH_B);
    const uint64_t *cxp(nullptr);
    if(extra) {
        ctx_sk = clade_sketch(sks, extra);
        cxp = ctx_sk.core().data();
    }
    long double sum(0);
    for(size_t i(0); i < m; ++i) {
        uint64_t v(0);
        for(u64 mm(mask); mm; mm &= mm - 1) {
            const uint64_t r(sks[__builtin_ctzll(mm)].core()[i]);
            if(r == EMPTY) {v = EMPTY; break;}
            v = std::max(v, r);
        }
        if(v != EMPTY && cxp) {
            const uint64_t r(cxp[i]);
            if(r == EMPTY) v = EMPTY; else v = std::max(v, r);
        }
        if(v != EMPTY) sum += exp2(-std::ldexp((long double)v, -int(64 - p)));
    }
    if(sum <= 0.) return 0.;
    const long double mean(sum / static_cast<long double>(m));
    return double(std::ldexp((long double)m, 0) * 0.6931471805599453L *
                  std::exp((std::log(mean) - std::lgammal(static_cast<long double>(k) + 1.0L)) / k));
}

// Exact per-node covers: the node's full sorted distinct-k-mer set. Small
// pangenomes fit them (round 7): acceptance error becomes 0 by construction
// (the sketch estimator is skipped entirely for covered nodes). Covers are
// PRODUCTION-BUILD material (the builder streams every genome's k-mers
// anyway); the bench builds them from the raw genomes for evaluation.
struct ExactCovers {
    bns::umap<u32, std::shared_ptr<std::vector<uint64_t>>> full;   // the node's whole distinct set (union)
    bns::umap<u32, std::shared_ptr<std::vector<uint64_t>>> core;   // k-mers in EVERY genome of the node (intersection)
    const std::vector<uint64_t> *get(u32 node) const {
        const auto it = full.find(node);
        return it == full.end() ? nullptr: it->second.get();
    }
    size_t bytes() const {
        bns::uset<const void *> seen;
        size_t n(0);
        for(const auto &kv: full) if(seen.insert(kv.second.get()).second) n += kv.second->size() * sizeof(uint64_t);
        for(const auto &kv: core) if(seen.insert(kv.second.get()).second) n += kv.second->size() * sizeof(uint64_t);
        return n;
    }
};

// Exact differential: k-mers present in EVERY part clade but not in every
// parent-clade genome: |inter parts| - |inter parent|. Set operations over
// the sorted covers; parent intersection is a subset of every part
// intersection, so this is a difference of sizes (no key walk).
inline double exact_differential(const ExactCovers &covers, const PlacementTree &tree,
                                 const std::vector<u64> &part_masks, u64 parent_mask) {
    if(part_masks.size() < 2) return 0.;
    std::vector<const std::vector<uint64_t> *> sets;
    for(const u64 pm: part_masks) {
        u32 any(u32(-1));
        for(u64 m(pm); m; m &= m - 1)
            if(any == u32(-1)) {any = u32(__builtin_ctzll(m)); break;}
        u32 anc(u32(-1));
        for(u32 x(tree.leaf_node(any)); x != u32(-1); x = tree.parent_of(x))
            if((tree.mask(x) & pm) == pm) {anc = x; break;}
        const auto *c(covers.get(anc));
        if(!c) return -1.;   // uncovered node: caller falls back to sketch
        sets.push_back(c);
    }
    // |intersection| via successive set_intersection on the smallest first
    std::vector<uint64_t> acc;
    {
        std::vector<size_t> ord(sets.size());
        std::iota(ord.begin(), ord.end(), size_t(0));
        std::sort(ord.begin(), ord.end(), [&](size_t a, size_t b) {return sets[a]->size() < sets[b]->size();});
        acc = *sets[ord[0]];
        std::vector<uint64_t> tmp;
        for(size_t i(1); i < ord.size() && !acc.empty(); ++i) {
            tmp.clear();
            std::set_intersection(acc.begin(), acc.end(), sets[ord[i]]->begin(), sets[ord[i]]->end(), std::back_inserter(tmp));
            acc.swap(tmp);
        }
    }
    // parent core: k-mers in EVERY parent genome (present in all parent
    // genomes implies present in all part genomes, so core(parent) is a
    // subset of inter(parts) and the differential is a size difference)
    u32 pany(u32(-1));
    for(u64 m(parent_mask); m; m &= m - 1) if(pany == u32(-1)) {pany = u32(__builtin_ctzll(m)); break;}
    u32 panc(u32(-1));
    for(u32 x(tree.leaf_node(pany)); x != u32(-1); x = tree.parent_of(x))
        if((tree.mask(x) & parent_mask) == parent_mask) {panc = x; break;}
    const auto pcit(covers.core.find(panc));
    if(pcit == covers.core.end() || !pcit->second) return double(acc.size());
    const auto &pc(*pcit->second);
    std::vector<uint64_t> pin;
    std::set_intersection(acc.begin(), acc.end(), pc.begin(), pc.end(), std::back_inserter(pin));
    return double(acc.size() - pin.size());
}

// Fold a p-register sketch down by 2^fold buckets: the folded register is
// the bucket-wise min, which is EXACTLY the sketch that would have been
// built at p - fold (one-permutation bucket nesting). Replication across
// register budgets therefore needs only the top-p sketch.
inline void fold_core(const uint64_t *src, size_t n, unsigned fold, uint64_t *dst) {
    const size_t out(n >> fold);
    const uint64_t EMPTY(std::numeric_limits<uint64_t>::max());
    for(size_t i(0); i < out; ++i) {
        uint64_t v(EMPTY);
        for(size_t j(0); j < (size_t(1) << fold); ++j) v = std::min(v, src[(i << fold) + j]);
        dst[i] = v;
    }
}

// Differential count on PREBUILT clade sketches (possibly folded), parts
// as masks; ctx given as per-genome leaf sketches (also foldable). Same
// estimator as pair_exclusive_count, on supplied register arrays.
inline double differential_from_sketches(const uint64_t *a, const uint64_t *b,
                                         const std::vector<const uint64_t *> &cx,
                                         size_t m, unsigned p) {
    const uint64_t EMPTY(std::numeric_limits<uint64_t>::max());
    const long double numinv(1. / std::ldexp((long double)1., 64 - p));
    size_t nshared(0);
    long double tmp(0);
    for(size_t i(0); i < m; ++i) {
        const uint64_t va(a[i]), vb(b[i]);
        const uint64_t v(std::min(va, vb));
        tmp = std::fma(numinv, v != EMPTY ? (long double)v: (long double)(EMPTY >> p), tmp);
        if(va == vb && va != EMPTY) {
            bool diverges(cx.empty());
            for(const uint64_t *g: cx) if(g[i] != va) {diverges = true; break;}
            nshared += diverges;
        }
    }
    return nshared * double(m) / double(tmp);
}

// Direct pairwise DIFFERENTIAL count: distinct k-mers shared by clades A
// and B that are NOT core to the parent clade (at least one ctx genome
// diverges from the common minimum). A bucket qualifies when A and B agree
// on a non-empty value v and some ctx genome's register differs from v.
// Pure pair-private sharing (reg_C > v for all C) is the special case; the
// filter also catches {A,B,C} accessory islands that pair-private measures
// miss (measured blind spot on the E. coli pangenome: subset-shared islands
// dominate, pair-exclusive sees nothing). Noise scales with the
// differential count itself - unlike difference-of-large-intersections,
// which at pangenome scale drowns a 10^3-kmer island in +-10^3 noise from
// the 10^6-kmer core (doc/realdata-grafting.md). A ctx register below v
// means that genome has a smaller bucket key: it may still contain the
// shared key, so the count overestimates slightly - documented bias,
// acceptable for edge formation and gain ranking.
inline double pair_exclusive_count(const std::vector<GenomeSketch> &sks, u64 ma, u64 mb, u64 ctx) {
    if(!ma || !mb) return 0.;
    const RegisterSketch A(clade_sketch(sks, ma)), B(clade_sketch(sks, mb));
    const auto &a(A.core());
    const auto &b(B.core());
    const size_t m(a.size());
    const uint64_t EMPTY(std::numeric_limits<uint64_t>::max());
    std::vector<const uint64_t *> cx;
    for(u64 mm(ctx); mm; mm &= mm - 1) cx.push_back(sks[__builtin_ctzll(mm)].core().data());
    const long double numinv(1. / std::ldexp((long double)1., 64 - A.getp()));
    size_t nshared(0);
    long double tmp(0);
    for(size_t i(0); i < m; ++i) {
        const uint64_t va(a[i]), vb(b[i]);
        const uint64_t v(std::min(va, vb));
        tmp = std::fma(numinv, v != EMPTY ? (long double)v: (long double)(EMPTY >> A.getp()), tmp);
        if(va == vb && va != EMPTY) {
            bool diverges(cx.empty());
            for(const uint64_t *g: cx) if(g[i] != va) {diverges = true; break;}
            nshared += diverges;
        }
    }
    return nshared * double(m) / double(tmp);
}

// Sketch-only gain proxy for graft A against its current parent clade P,
// from the validated PAIR intersection estimator alone (the k-way one-
// permutation order-statistic estimators fail at these resolutions -
// measured in doc/sketch-only-naming.md: with m = 8192 registers and
// n ~ 8-60k keys per genome (n ~ m or denser), the k-th-min estimator
// returns ~2063 for a DISJOINT pair (true ~0): max-of-minima cannot tell
// "small intersection" from "no intersection", and per-bucket occupancy
// breaks the exponential-spacing model). A is given as its component's
// member child-clade masks: the shared core is the minimum pairwise
// intersection over member clades; the contamination is the maximum
// pairwise intersection of a member with the context clade P \ A. gain =
// core - contam, the k-mers the component shares with ITSELF beyond what
// the surrounding clade explains.
inline double sketch_exclusive_estimate(const std::vector<GenomeSketch> &sks, const std::vector<u64> &parts, u64 parent) {
    const u64 A(std::accumulate(parts.begin(), parts.end(), u64(0), [](u64 a, u64 b) {return a | b;}));
    const u64 ctx(parent & ~A);
    if(parts.size() < 2 || !ctx) return 0.;
    double gain(std::numeric_limits<double>::max());
    for(size_t i(0); i < parts.size(); ++i)
        for(size_t j(i + 1); j < parts.size(); ++j)
            gain = std::min(gain, pair_exclusive_count(sks, parts[i], parts[j], ctx));
    return std::max(0., gain);
}

// Round 5 (doc/grafting-two-tier.md): two-tier per-node sketches and
// adaptive tau.
enum class NodeTier: uint8_t {SPARSE, DENSE};   // EXACT covers are plumbed via providers; without one, exact-tier nodes use their raw sketch
struct SketchNodePolicy {
    unsigned exact_max = 512;      // clade n <= this: exact-cover tier in production
    unsigned dense_min = 25000;    // clade n >= this: estimate on a densified COPY (union rule: never pre-densify)
    double tau_k = 4.;             // adaptive tau: background + k * sigma
    double tau_floor = 32.;        // when the node's background is ~0
    bool exclusive_edges = false;  // real pangenomes: edge on the pair's PRIVATE sharing (ctx-cancelled) -
                                   // raw counts are core-genome dominated and no module can clear them
    bool pair_candidates = false;  // emit each qualifying pair as its own graft (no union-find):
                                   // pangenome differential signals chain single-linkage into whole-clade
                                   // components; greedily applied pairs are the robust unit
};
inline NodeTier node_tier(double n_est, const SketchNodePolicy &pol) {
    return n_est >= double(pol.dense_min) ? NodeTier::DENSE: NodeTier::SPARSE;
}

// Per-node sketch index: raw clade-union sketches (register-wise min,
// sentinel = uint64 max = +inf; NEVER pre-densified), fill-based clade
// cardinality estimates (linear counting - sketch-only; production should
// prefer the BNSKMC01 sidecar folded up the tree), tiers, and the adaptive
// tau per internal node (median pairwise child intersection + k sigma).
struct SketchIndex {
    const PlacementTree *tree = nullptr;
    const std::vector<GenomeSketch> *sks = nullptr;
    const ExactCovers *covers = nullptr;   // round 7: exact pangenome mode (acceptance error 0 where covered)
    SketchNodePolicy pol;
    std::vector<u64> clade_mask;
    std::vector<RegisterSketch> clade_sk;
    std::vector<double> n_est;
    std::vector<NodeTier> tier;
    std::vector<double> tau;      // per node; 0 for leaves
};

// Acceptance-scale differential: exact covers when available (error 0),
// sketch estimate otherwise. Single entry point so edges, taus and gains
// never mix units.
inline double node_differential(const SketchIndex &idx, const std::vector<u64> &parts, u64 parent) {
    if(idx.covers) {
        const double e(exact_differential(*idx.covers, *idx.tree, parts, parent));
        if(e >= 0.) return e;
    }
    return sketch_exclusive_estimate(*idx.sks, parts, parent);
}

inline SketchIndex build_sketch_index(const PlacementTree &tree, const std::vector<GenomeSketch> &sks,
                                      const SketchNodePolicy &pol = SketchNodePolicy{}) {
    const u32 n(u32(sks.size()));
    if(n == 0 || n != tree.ngenomes()) throw std::invalid_argument("sketch/genome count mismatch");
    SketchIndex ret;
    ret.tree = &tree;
    ret.sks = &sks;
    ret.pol = pol;
    const uint64_t EMPTY(std::numeric_limits<uint64_t>::max());
    ret.clade_mask.assign(tree.nnodes(), 0);
    for(u32 g(0); g < n; ++g) ret.clade_mask[tree.leaf_node(g)] = u64(1) << g;
    std::function<void(u32)> post = [&](u32 x) {
        if(tree.is_leaf_node(x)) return;
        for(const u32 c: tree.children_of(x)) {
            post(c);
            ret.clade_mask[x] |= ret.clade_mask[c];
        }
    };
    post(tree.root());
    ret.clade_sk.reserve(tree.nnodes());
    for(u32 x(0); x < tree.nnodes(); ++x)
        ret.clade_sk.push_back(tree.is_leaf_node(x) && !ret.clade_mask[x]
                               ? RegisterSketch(sks[0].getp(), GENOME_SKETCH_B)
                               : clade_sketch(sks, ret.clade_mask[x]));
    // fill-based clade cardinality (sketch-only tiering input)
    ret.n_est.assign(tree.nnodes(), 0.);
    const double m(double(ret.clade_sk[0].core().size()));
    for(u32 x(0); x < tree.nnodes(); ++x) {
        size_t empty(0);
        for(const auto v: ret.clade_sk[x].core()) empty += v == EMPTY;
        const double f(1. - double(empty) / m);
        ret.n_est[x] = f > 0. ? -m * std::log(1. - f): 0.;
    }
    ret.tier.resize(tree.nnodes());
    for(u32 x(0); x < tree.nnodes(); ++x) ret.tier[x] = node_tier(ret.n_est[x], pol);
    // adaptive tau: median pairwise child intersection + k sigma
    ret.tau.assign(tree.nnodes(), 0.);
    for(u32 v(0); v < tree.nnodes(); ++v) {
        if(!tree.is_alive(v) || tree.is_leaf_node(v)) continue;
        const auto &ch(tree.children_of(v));
        if(ch.size() < 2) continue;
        std::vector<double> counts;
        for(size_t i(0); i < ch.size(); ++i)
            for(size_t j(i + 1); j < ch.size(); ++j)
                counts.push_back(pol.exclusive_edges
                                 ? node_differential(ret, std::vector<u64>{ret.clade_mask[ch[i]], ret.clade_mask[ch[j]]}, ret.clade_mask[v])
                                 : pair_intersection_count(ret.clade_sk[ch[i]], ret.clade_sk[ch[j]]));
        std::sort(counts.begin(), counts.end());
        const double bg(counts.size() & 1 ? counts[counts.size() / 2]
                                          : 0.5 * (counts[counts.size() / 2 - 1] + counts[counts.size() / 2]));
        const double sigma(std::sqrt(std::max(bg, 1.)));
        // RAW mode: the median is background shared by every child - add a
        // noise band on top. DIFFERENTIAL mode: the filter already cancelled
        // the core; the median is SIGNAL, so the threshold is a pure
        // counting-noise band (k x sqrt(median)), floored.
        ret.tau[v] = pol.exclusive_edges ? std::max(pol.tau_floor, pol.tau_k * sigma)
                                         : (bg > 0.5 ? bg + pol.tau_k * sigma: pol.tau_floor);
    }
    return ret;
}

// Tier-aware pair count: DENSE-tier sketches are densified ONCE, on a
// copy, at estimation time (union rule); SPARSE stays plain (2-8x sharper
// at small n). Exact covers (production providers) would slot in here.
inline double tier_pair_count(const SketchIndex &idx, u32 a, u32 b) {
    if(idx.tier[a] == NodeTier::DENSE || idx.tier[b] == NodeTier::DENSE) {
        RegisterSketch x(idx.sks->at(0).getp(), GENOME_SKETCH_B), y(idx.sks->at(0).getp(), GENOME_SKETCH_B);
        x.mutable_core() = idx.clade_sk[a].core();
        y.mutable_core() = idx.clade_sk[b].core();
        x.densify(); y.densify();
        return pair_intersection_count(x, y);
    }
    return pair_intersection_count(idx.clade_sk[a], idx.clade_sk[b]);
}

struct SketchGraft {
    u64 mask;            // graft leaf set
    std::vector<u64> parts;   // member child-clade masks (gain decomposition)
    double est_gain;     // sketch-only exclusive-count estimate
    double tau;          // acceptance threshold at enumeration (per-node adaptive, or absolute)
    u32 node;            // node whose children seeded it (gain context)
};

/*
Sketch-only seeded hulls: at each internal node v, cluster v's children by
their clade-sketch intersection profile and emit each connected component
of size >= 2 as a graft mask with its sketch-only gain against v's clade.
O(sum_v k_v^2) pairwise intersections of O(2^p) register work. No k-mer
sets anywhere.

Edge design (measured, doc/sketch-only-naming.md): raw pairwise-
intersection threshold (the validated submodule estimator). Threshold
sensitivity is the lane's central limitation: a shared clade background
(a genus block in every species) chains single-linkage into one whole-
clade component (skipped, no graft needed - measured on the GTDB-shaped
set), while a module ABOVE that background only separates at a tau between
background and background+module, and a module BELOW the largest
background anywhere is invisible. Absolute tau cannot fit both regimes;
the bench reports a tau grid.
*/
inline std::vector<SketchGraft> enumerate_sketch_seeded_grafts(const SketchIndex &idx, double tau_abs = -1.) {
    const PlacementTree &tree(*idx.tree);
    std::vector<SketchGraft> ret;
    for(u32 v(0); v < tree.nnodes(); ++v) {
        if(!tree.is_alive(v) || tree.is_leaf_node(v)) continue;
        const auto &ch(tree.children_of(v));
        const size_t k(ch.size());
        if(k < 2) continue;
        const double tau(tau_abs >= 0. ? tau_abs: idx.tau[v]);   // per-node adaptive when no absolute override
        std::vector<std::array<u64, 2>> edges;
        for(size_t i(0); i < k; ++i)
            for(size_t j(i + 1); j < k; ++j) {
                const double w(idx.pol.exclusive_edges
                               ? sketch_exclusive_estimate(*idx.sks, std::vector<u64>{idx.clade_mask[ch[i]], idx.clade_mask[ch[j]]}, idx.clade_mask[v])
                               : tier_pair_count(idx, ch[i], ch[j]));
                if(w >= tau) edges.push_back({u64(i), u64(j)});
            }
        if(idx.pol.pair_candidates) {
            for(const auto &e: edges) {
                const std::vector<u64> parts{idx.clade_mask[ch[size_t(e[0])]], idx.clade_mask[ch[size_t(e[1])]]};
                const double g(sketch_exclusive_estimate(*idx.sks, parts, idx.clade_mask[v]));
                ret.push_back(SketchGraft{parts[0] | parts[1], parts, g, tau, v});
            }
            continue;
        }
        std::vector<u32> par(k), sz(k, 1);
        std::iota(par.begin(), par.end(), 0);
        std::function<u32(u32)> find = [&](u32 x) {return par[x] == x ? x: par[x] = find(par[x]);};
        auto unite = [&](u32 a, u32 b) {
            a = find(a); b = find(b);
            if(a == b) return;
            if(sz[a] < sz[b]) std::swap(a, b);
            par[b] = a; sz[a] += sz[b];
        };
        for(const auto &e: edges) unite(u32(e[0]), u32(e[1]));
        bns::umap<u32, std::vector<u64>> comps;
        for(size_t i(0); i < k; ++i) comps[find(u32(i))].push_back(idx.clade_mask[ch[i]]);
        for(const auto &kv: comps) {
            if(kv.second.size() < 2) continue;
            const u64 cm(std::accumulate(kv.second.begin(), kv.second.end(), u64(0), [](u64 a, u64 b) {return a | b;}));
            if(cm == idx.clade_mask[v]) continue;   // the whole clade: nothing to graft
            const double g(node_differential(idx, kv.second, idx.clade_mask[v]));
            ret.push_back(SketchGraft{cm, kv.second, g, tau, v});
        }
    }
    return ret;
}
inline std::vector<SketchGraft> enumerate_sketch_seeded_grafts(const PlacementTree &tree, const std::vector<GenomeSketch> &sks, double tau) {
    return enumerate_sketch_seeded_grafts(build_sketch_index(tree, sks), tau);
}

/*
Sketch-only greedy: selection and gains are the bbminhash intersection
estimates alone (recomputed each round against the CURRENT tree - parent
clades shift as grafts are applied); add_synthetic still validates
feasibility. The caller evaluates the objective (exact H is the oracle in
benches; at scale an HLL/union sketch of the final tree plays that role).
*/
inline GraftingResult graft_greedy_sketch(PlacementTree &tree, const std::vector<GenomeSketch> &sks,
                                          std::vector<SketchGraft> candidates, u32 budget,
                                          const ExactCovers *covers = nullptr) {
    GraftingResult ret;
    ret.initial_delta = 0;   // no index: the caller owns the objective
    while(ret.added_masks.size() < budget && !candidates.empty()) {
        // recompute the gain of every remaining candidate against the
        // current tree: parent clade = the mask of the deepest live node
        // still containing it (candidates were seeded at `node`; its clade
        // may have shrunk).
        for(auto &c: candidates) {
            u32 anc(c.node);
            while(anc != u32(-1) && (tree.mask(anc) & c.mask) != c.mask) anc = tree.parent_of(anc);
            const u64 parent(anc == u32(-1) ? tree.genome_mask(): tree.mask(anc));
            if(covers) {
                const double e(exact_differential(*covers, tree, c.parts, parent));
                c.est_gain = e >= 0. ? e: sketch_exclusive_estimate(sks, c.parts, parent);
            } else c.est_gain = sketch_exclusive_estimate(sks, c.parts, parent);
        }
        size_t best(size_t(-1));
        for(size_t i(0); i < candidates.size(); ++i)
            if(candidates[i].est_gain > 0. && candidates[i].est_gain >= candidates[i].tau
               && (best == size_t(-1) || candidates[i].est_gain > candidates[best].est_gain)) best = i;
        if(best == size_t(-1)) break;
        const u64 A(candidates[best].mask);
        if(tree.add_synthetic(A) < 0) {candidates.erase(candidates.begin() + best); continue;}
        ret.added_masks.push_back(A);
        ret.est_gains.push_back(candidates[best].est_gain);
        candidates.erase(candidates.begin() + best);
        // drop masks that are no longer laminar with the tree (their gain
        // recomputation parent would be wrong); add_synthetic rejects them
        // anyway, so let the next round's failures clean up.
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
