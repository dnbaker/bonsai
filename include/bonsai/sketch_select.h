#ifndef __BNS_SKETCH_SELECT_H
#define __BNS_SKETCH_SELECT_H

/*
Sketching-accelerated structure selection for the synthetic-node
placement technique (placement.h).

Two accelerations are validated against the exact greedy:

1. Candidate generation from MinHash sketches: per-genome
   RangeMinHash sketches over canonical k-mers; genomes are clustered
   by pairwise Jaccard (bottom-k estimator) at >= a threshold, and the
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
#include "sketch/mh.h"

#include <numeric>
#include <random>

namespace bns {

// Per-genome RangeMinHash sketches over canonical k-mers.
template<typename ScoreType>
std::vector<sketch::RangeMinHash<uint64_t>::final_type> genome_sketches(const std::vector<std::string> &paths, const Spacer &sp, size_t sketch_size=2048) {
    std::vector<sketch::RangeMinHash<uint64_t>::final_type> ret;
    ret.reserve(paths.size());
    for(const auto &path: paths) {
        sketch::RangeMinHash<uint64_t> rmh(sketch_size);
        Encoder<ScoreType> enc(0, 0, sp, nullptr, true);
        enc.for_each([&](u64 v) {rmh.addh(v);}, path.data());
        ret.push_back(rmh.finalize());
    }
    return ret;
}

// Connected components of the graph linking genome pairs whose sketches
// agree at >= jaccard_threshold; component masks (size >= 2) are the
// candidate placements.
inline std::vector<u64> sketch_cluster_candidates(const std::vector<sketch::RangeMinHash<uint64_t>::final_type> &sketches, double jaccard_threshold) {
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
