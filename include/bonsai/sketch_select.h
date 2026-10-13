#ifndef __BNS_SKETCH_SELECT_H
#define __BNS_SKETCH_SELECT_H

/*
Sketching-accelerated structure selection for the synthetic-node
placement technique (placement.h).

Two accelerations are validated against the exact greedy:

1. Candidate generation from MinHash sketches: per-genome
   BBMinHash sketches (BBitMinHasher, b = 16 retained bits per register)
   over canonical k-mers; genomes are clustered
   by pairwise Jaccard (b-bit estimator) at >= a threshold, and the
   connected-component masks replace the exact observed-bitmap
   candidates. On module-structured sets the components recover the
   module masks; on iid-random sets nothing clusters, so the greedy
   still adds nothing.

2. Hash-subsampled disagreement metric: singleton k-mer classes are
   kept with probability p and scaled by 1/p (Horvitz-Thompson), while
   structural (shared) k-mers are forced kept, so greedy gains are
   unchanged and only the reported mean carries sampling error.
   subsample_kmers(force_structural=false) additionally subsamples
   shared k-mer copies so target bitmaps shrink: the placement-distortion
   regime.
*/

#include "placement.h"
#include "sketch/bbmh.h"

#include <numeric>
#include <optional>
#include <memory>
#include <random>

namespace bns {

/*
Per-genome BBMinHash sketches (sketch/bbmh.h): BBitMinHasher with 64-bit
register words, b = 16 retained quotient bits per register ("16-bit
registers"). sketch_size counts registers and must be a power of 2 (it sets
p). The live register file is 2^p x 8 bytes; the finalized sketch packs
b x 2^p bits, so at b = 16 the stored sketch holds 2 bytes per register
and the default 8192 registers below match the stored size of the previous
RangeMinHash<2048> default (2048 x 64-bit minima = 16 KiB). b = 16 keeps
the chance-agreement floor (2^-b) and the cross-set quotient-collision bias
(n^2 x 2^-b / 2^p) negligible, which low-Jaccard candidate thresholds need;
smaller b floods both. finalize() densifies sparse sketches, so genomes
that fill fewer than 2^p buckets keep an unbiased estimate.
*/
using GenomeSketch = sketch::BBitMinHasher<uint64_t>;
static constexpr unsigned GENOME_SKETCH_B = 16;

template<typename ScoreType>
std::vector<GenomeSketch::final_type> genome_sketches(const std::vector<std::string> &paths, const Spacer &sp, size_t sketch_size=8192) {
    if(!sketch::is_pow2(sketch_size)) throw std::invalid_argument("sketch_size must be a power of 2 for BBitMinHasher");
    const unsigned p(sketch::ilog2(sketch_size));
    std::vector<GenomeSketch::final_type> ret;
    ret.reserve(paths.size());
    for(const auto &path: paths) {
        GenomeSketch bb(p, GENOME_SKETCH_B);
        Encoder<ScoreType> enc(0, 0, sp, nullptr, true);
        enc.for_each([&](u64 v) {bb.addh(v);}, path.data());
        ret.push_back(bb.finalize());
    }
    return ret;
}

// Connected components of the graph linking genome pairs whose sketches
// agree at >= jaccard_threshold; component masks (size >= 2) are the
// candidate placements. All-pairs union-find; the correctness reference
// for small instances (see sketch_cluster_candidates_tree for the
// taxonomic-layered path used at scale).
inline std::vector<u64> sketch_cluster_candidates(const std::vector<GenomeSketch::final_type> &sketches, double jaccard_threshold) {
    const size_t n(sketches.size());
    std::vector<u32> par(n), sz(n, 1);
    std::iota(par.begin(), par.end(), 0);
    std::function<u32(u32)> find = [&](u32 x) {return par[x] == x ? x: par[x] = find(par[x]);};
    auto unite = [&](u32 a, u32 b) {
        a = find(a); b = find(b);
        if(a == b) return;
        if(sz[a] < sz[b]) std::swap(a, b);
        par[b] = a; sz[a] += sz[b];
    };
    for(u32 i(0); i < n; ++i)
        for(u32 j(i + 1); j < n; ++j)
            if(sketches[i].jaccard_index(sketches[j]) >= jaccard_threshold) unite(i, j);
    std::unordered_map<u32, u64> comp_mask;
    for(u32 i(0); i < n; ++i) comp_mask[find(i)] |= u64(1) << i;
    std::vector<u64> ret;
    for(const auto &pair: comp_mask) if(pc64(pair.second) >= 2) ret.push_back(pair.second);
    return ret;
}

// Unfinalized per-genome sketches, for merge-based use: an internal tree
// node's sketch is the register-wise min merge (BBitMinHasher::operator+=)
// of its children's sketches, so clade sketches at any depth are built from
// O(2^p) register work per merge instead of re-encoding raw k-mers.
template<typename ScoreType>
std::vector<GenomeSketch> genome_sketchers(const std::vector<std::string> &paths, const Spacer &sp, size_t sketch_size=8192) {
    if(!sketch::is_pow2(sketch_size)) throw std::invalid_argument("sketch_size must be a power of 2 for BBitMinHasher");
    const unsigned p(sketch::ilog2(sketch_size));
    std::vector<GenomeSketch> ret;
    ret.reserve(paths.size());
    for(const auto &path: paths) {
        ret.emplace_back(p, GENOME_SKETCH_B);
        Encoder<ScoreType> enc(0, 0, sp, nullptr, true);
        enc.for_each([&](u64 v) {ret.back().addh(v);}, path.data());
    }
    return ret;
}

/*
Layered candidate generation over the taxonomic tree. Bottom-up, each
internal node's sketch is the register-wise merge of its children's, and
sibling clades whose merged-sketch Jaccard clears the threshold union their
genome sets (single linkage). Comparisons stay inside sibling groups:
O(sum_v k_v^2) clade-pair Jaccards of O(2^p) register operations each,
versus O(n^2) genome-pair Jaccards for the all-pairs reference path above.
A genome pair whose sketches agree may be missed when the agreement is
diluted across its sibling clades, so this is the scale path; the all-pairs
path remains the correctness reference. Returns connected components
(size >= 2) as sorted genome-id lists; no per-genome k-mer bitmaps are
materialized anywhere in the pipeline.
*/
inline std::vector<std::vector<u32>> sketch_cluster_candidates_tree(const std::vector<GenomeSketch> &leaf_sketches, const PlacementTree &tree, double jaccard_threshold) {
    const u32 n(u32(leaf_sketches.size()));
    if(n == 0 || n != tree.ngenomes()) throw std::invalid_argument("sketch/genome count mismatch");
    for(const auto &s: leaf_sketches) if(s.getp() != leaf_sketches[0].getp()) throw std::invalid_argument("mixed register counts");
    std::vector<u32> par(n), sz(n, 1);
    std::iota(par.begin(), par.end(), 0);
    std::function<u32(u32)> find = [&](u32 x) {return par[x] == x ? x: par[x] = find(par[x]);};
    auto unite = [&](u32 a, u32 b) {
        a = find(a); b = find(b);
        if(a == b) return;
        if(sz[a] < sz[b]) std::swap(a, b);
        par[b] = a; sz[a] += sz[b];
    };
    const u32 nn(tree.nnodes());
    std::vector<std::optional<GenomeSketch>> merged(nn);            // union sketch per node
    std::vector<std::unique_ptr<GenomeSketch::final_type>> fin(nn); // finalized once per node (move-only)
    std::vector<std::vector<u32>> members(nn);                      // descendant genome ids
    std::vector<bool> have_fin(nn, false), have_mrg(nn, false);
    // genome id -> leaf node is 1:1; invert it.
    std::vector<u32> leaf_of(n);
    for(u32 g(0); g < n; ++g) leaf_of[g] = tree.leaf_node(g);
    std::vector<u32> genome_of_leaf(nn, u32(-1));
    for(u32 g(0); g < n; ++g) genome_of_leaf[leaf_of[g]] = g;
    std::function<void(u32)> post = [&](u32 x) {
        if(tree.is_leaf_node(x)) {
            // leaf node: copy the genome's sketch
            const u32 g(genome_of_leaf[x]);
            if(g != u32(-1)) {
                merged[x] = leaf_sketches[g];
                members[x] = {g};
            }
            have_mrg[x] = true;
            return;
        }
        for(const auto c: tree.children_of(x)) {
            post(c);
            if(!have_fin[c]) {fin[c] = std::make_unique<GenomeSketch::final_type>(merged[c]->finalize()); have_fin[c] = true;}
            if(!have_mrg[x]) {
                merged[x] = std::move(*merged[c]);
                have_mrg[x] = true;
            } else {
                *merged[x] += *merged[c];
            }
            for(auto &g: members[c]) members[x].push_back(g);
        }
        fin[x] = std::make_unique<GenomeSketch::final_type>(merged[x]->finalize());
        have_fin[x] = true;
        // sibling-clade comparisons at this node
        const auto &ch(tree.children_of(x));
        for(size_t i(0); i < ch.size(); ++i) {
            for(size_t j(i + 1); j < ch.size(); ++j) {
                if((*fin[ch[i]]).jaccard_index(*fin[ch[j]]) >= jaccard_threshold) {
                    if(!members[ch[i]].empty() && !members[ch[j]].empty())
                        unite(members[ch[i]].front(), members[ch[j]].front());
                }
            }
        }
    };
    post(tree.root());
    std::vector<u32> order(n);
    std::iota(order.begin(), order.end(), 0);
    std::sort(order.begin(), order.end(), [&](u32 a, u32 b){return find(a) < find(b);});
    std::vector<std::vector<u32>> ret;
    for(size_t i(0); i < order.size();) {
        size_t j(i + 1);
        while(j < order.size() && find(order[j]) == find(order[i])) ++j;
        if(j - i >= 2) ret.push_back(std::vector<u32>(order.begin() + i, order.begin() + j));
        i = j;
    }
    return ret;
}

// Horvitz-Thompson singleton subsampling of an exact KmerClasses:
// structural (shared) classes are kept exactly; each singleton class is
// kept with probability p and scaled by 1/p. Greedy gains are identical
// to the exact classes; only mean estimates carry sampling error.
inline KmerClasses subsample_singletons(const KmerClasses &kc, double p, uint64_t seed) {
    KmerClasses ret;
    ret.ngenomes = kc.ngenomes;
    std::mt19937_64 rng(seed);
    std::uniform_real_distribution<double> uni(0., 1.);
    for(const auto &pair: kc.count) {
        const bool structural(pc64(pair.first) >= 2);
        double scale(1.);
        if(!structural) {
            if(uni(rng) >= p) continue;
            scale = 1. / p;
        }
        const u64 newcount(u64(std::lround(pair.second * scale)));
        if(!newcount) continue;
        ret.count[pair.first] = u32(newcount);
        ret.weight[pair.first] = u64(std::llround(kc.weight.at(pair.first) * scale));
        ret.total_distinct += newcount;
    }
    return ret;
}

/*
Occurrence-level hash-prefix subsampling from scratch. Each k-mer copy is
kept with probability p and scaled by 1/p, except that k-mers shared in
the full data keep all their copies when force_structural is set (exact
target bitmaps for the shared structure; singleton classes are estimated).
With force_structural=false, shared copies are subsampled too, so target
bitmaps shrink and placements can distort.
*/
template<typename ScoreType>
KmerClasses subsample_kmers(const std::vector<std::string> &paths, const Spacer &sp, double p, uint64_t seed, bool force_structural, double *kept_fraction=nullptr) {
    const u32 n(paths.size());
    std::unordered_map<u64, u32> occ;
    u64 total_copies(0);
    for(const auto &path: paths) {
        Encoder<ScoreType> enc(0, 0, sp, nullptr, true);
        enc.for_each([&](u64 v) {++occ[v]; ++total_copies;}, path.data());
    }
    // Sample at occurrence granularity; remember each kept k-mer's scaled
    // occurrence count. Its class is final once all genomes are scanned.
    std::unordered_map<u64, u64> bm;              // k-mer -> genome bits
    std::unordered_map<u64, long double> kocc;    // k-mer -> scaled occurrences
    std::unordered_map<u64, long double> kcnt;    // k-mer -> scaled distinct-count contribution
    u64 kept_copies(0);
    std::mt19937_64 rng(seed);
    std::uniform_real_distribution<double> uni(0., 1.);
    for(u32 i(0); i < n; ++i) {
        const u64 bit(u64(1) << i);
        Encoder<ScoreType> enc(0, 0, sp, nullptr, true);
        enc.for_each([&](u64 v) {
            double scale(1.);
            if(!(force_structural && occ[v] > 1)) {
                if(uni(rng) >= p) return;
                scale = 1. / p;
            }
            ++kept_copies;
            if(bm[v] == 0) kcnt[v] = scale;
            bm[v] |= bit;
            kocc[v] += scale;
        }, paths[i].data());
    }
    KmerClasses ret;
    ret.ngenomes = n;
    std::unordered_map<u64, long double> class_count, class_weight;
    for(const auto &pair: bm) {
        class_count[pair.second] += kcnt[pair.first];
        class_weight[pair.second] += kocc[pair.first];
    }
    for(const auto &pair: class_count) {
        const u64 c(u64(std::llround(pair.second)));
        ret.count[pair.first] = u32(c);
        ret.weight[pair.first] = u64(std::llround(class_weight[pair.first]));
        ret.total_distinct += c;
    }
    if(kept_fraction) *kept_fraction = total_copies ? double(kept_copies) / total_copies : 0.;
    return ret;
}

} // namespace bns

#endif // __BNS_SKETCH_SELECT_H
