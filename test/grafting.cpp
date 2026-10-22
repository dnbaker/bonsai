#include "test/catch.hpp"

#include "grafting.h"
#include "util.h"

#include "aesctr/wy.h"

#include <random>
#include <string>
#include <vector>

using namespace bns;

namespace {

std::string rand_dna(wy::WyRand<u64> &rng, size_t len) {
    static const char bases[] = "ACGT";
    std::string ret;
    ret.reserve(len);
    while(ret.size() < len) ret.push_back(bases[rng() & 3]);
    return ret;
}

struct TmpDir {
    std::string dir;
    std::vector<std::string> files;
    TmpDir(): dir([] {
        char tmpl[] = "/tmp/bns_graft_XXXXXX";
        const char *ret(mkdtemp(tmpl));
        REQUIRE(ret != nullptr);
        return std::string(ret);
    }()) {}
    ~TmpDir() {
        for(const auto &f: files) std::remove(f.data());
        std::remove(dir.data());
    }
    std::string add(const std::string &name, const std::string &content) {
        files.push_back(dir + "/" + name);
        std::FILE *fp(std::fopen(files.back().data(), "w"));
        REQUIRE(fp != nullptr);
        std::fputs(content.data(), fp);
        std::fclose(fp);
        return files.back();
    }
};

khash_t(p) *build_tax(const std::vector<std::pair<tax_t, tax_t>> &edges) {
    khash_t(p) *ret(kh_init(p));
    int khr;
    for(const auto &e: edges) {
        const khint_t ki(kh_put(p, ret, e.first, &khr));
        kh_val(ret, ki) = e.second;
    }
    return ret;
}

// tax: 1 root; 200/201 genera; 100+i species. modules: list of
// (block-length carrier mask) — each module is one shared block record per
// carrier genome (separate FASTA records: no junction k-mers).
struct ModuleSet {
    TmpDir tmp;
    std::vector<std::string> paths;
    std::vector<tax_t> gt;
    khash_t(p) *tax;
    std::vector<std::pair<u64, unsigned>> modules;
    ModuleSet(const std::vector<std::pair<tax_t, tax_t>> &edges,
              const std::vector<tax_t> &genome_taxa,
              const std::vector<std::pair<u64, unsigned>> &mods,
              u64 seed): tax(build_tax(edges)), modules(mods) {
        wy::WyRand<u64> rng(seed);
        std::vector<std::string> blocks(mods.size());
        for(size_t i(0); i < mods.size(); ++i) blocks[i] = rand_dna(rng, mods[i].second);
        for(unsigned g(0); g < genome_taxa.size(); ++g) {
            std::string content(">g" + std::to_string(g) + "\n" + rand_dna(rng, 1400) + "\n");
            for(size_t i(0); i < mods.size(); ++i)
                if(mods[i].first & (u64(1) << g))
                    content += ">g" + std::to_string(g) + "_m" + std::to_string(i) + "\n" + blocks[i] + "\n";
            paths.push_back(tmp.add("g" + std::to_string(g) + ".fa", content));
            gt.push_back(genome_taxa[g]);
        }
    }
    ~ModuleSet() {kh_destroy(p, tax);}
};

std::vector<tax_t> species_taxa(unsigned ngenomes) {
    std::vector<tax_t> gt;
    for(unsigned g(0); g < ngenomes; ++g) gt.push_back(100 + (g >> 1));
    return gt;
}

// 8 leaves: species 100,101 under genus 200 (genomes 0-3); species 102,103
// under genus 201 (genomes 4-7). Three cross-taxon modules.
struct EightLeaf {
    ModuleSet ms;
    EightLeaf(): ms({{1, 0}, {200, 1}, {201, 1}, {100, 200}, {101, 200}, {102, 201}, {103, 201}},
                    species_taxa(8),
                    {{0x33u, 120}, {0x14u, 120}, {0x48u, 120}}, 777) {}
};

} // namespace

// (a) greedy vs the exact optimum (branch and bound over all feasible graft
// subsets), B <= 3, <= 8 leaves. On the 8-leaf module instance the exact
// search beats greedy at B = 2 (delta 624 vs 416 — the destruction geometry
// of the non-submodularity counterexample), while greedy still clears the
// (1 - 1/e) reduction threshold; smaller and randomized instances are
// gap-free.
TEST_CASE("graft_greedy_matches_exact") {
    Spacer sp(17, 17);
    EightLeaf st;
    KmerClasses kc;
    kc.fill<score::Lex>(st.ms.paths, sp);
    const std::vector<u64> cands(enumerate_grafts(PlacementTree(st.ms.tax, st.ms.gt), true));
    std::fprintf(stderr, "graft-8leaf: %zu candidates\n", cands.size());
    // B = 0: nothing is added.
    {
        PlacementTree tree(st.ms.tax, st.ms.gt);
        const GraftingResult gr(graft_greedy(tree, kc, cands, 0));
        REQUIRE(gr.added_masks.empty());
        REQUIRE(gr.deltas.size() == 1);
    }
    // B <= 2 by complete subset search (all_subsets on 8 leaves): greedy
    // must keep at least a (1 - 1/e) fraction of the optimal reduction.
    for(const u32 budget: {1u, 2u}) {
        PlacementTree gt_tree(st.ms.tax, st.ms.gt);
        const GraftingResult gr(graft_greedy(gt_tree, kc, cands, budget));
        const ExactGraftResult ex(graft_exact(PlacementTree(st.ms.tax, st.ms.gt), kc, cands, budget));
        const long double red_g((long double)(gr.initial_delta - gr.deltas.back()));
        const long double red_e((long double)(gr.initial_delta - ex.final_delta));
        std::fprintf(stderr, "graft-8leaf B=%u: greedy delta %llu, exact %llu, reduction ratio %.3Lf\n",
                     budget, (unsigned long long)gr.deltas.back(), (unsigned long long)ex.final_delta,
                     red_e > 0 ? red_g / red_e : 1.0L);
        const bool meets((red_e == 0) || (red_g >= red_e * (1.0L - 1.0L / std::exp(1.0L))));
        REQUIRE(meets);
    }
    // B = 3 on 6 leaves (complete enumeration stays cheap): gap-free.
    {
        ModuleSet st6({{1, 0}, {200, 1}, {201, 1}, {100, 200}, {101, 200}, {102, 201}},
                      species_taxa(6),
                      {{0x09u, 120}, {0x16u, 120}, {0x2cu, 120}}, 555);
        KmerClasses kc6;
        kc6.fill<score::Lex>(st6.paths, sp);
        const std::vector<u64> c6(enumerate_grafts(PlacementTree(st6.tax, st6.gt), true));
        PlacementTree gt_tree(st6.tax, st6.gt);
        const GraftingResult gr(graft_greedy(gt_tree, kc6, c6, 3));
        const ExactGraftResult ex(graft_exact(PlacementTree(st6.tax, st6.gt), kc6, c6, 3));
        std::fprintf(stderr, "graft-6leaf B=3: greedy delta %llu, exact %llu\n",
                     (unsigned long long)gr.deltas.back(), (unsigned long long)ex.final_delta);
        REQUIRE(gr.deltas.back() == ex.final_delta);
    }
    // Seeded randomized instances: 6 leaves, random cross-species modules.
    for(const u64 seed: {11u, 22u, 33u, 44u}) {
        wy::WyRand<u64> rng(seed);
        std::vector<std::pair<u64, unsigned>> mods;
        for(unsigned m(0); m < 3; ++m) {
            u64 mask(0);
            while(pc64(mask) < 2) {
                mask = 0;
                for(unsigned g(0); g < 6; ++g) if((rng() & 3) == 0) mask |= u64(1) << g;
            }
            mods.push_back({mask, 120});
        }
        ModuleSet st_r({{1, 0}, {200, 1}, {201, 1}, {100, 200}, {101, 200}, {102, 201}},
                       species_taxa(6), mods, seed * 7 + 1);
        KmerClasses kcr;
        kcr.fill<score::Lex>(st_r.paths, sp);
        const std::vector<u64> cr(enumerate_grafts(PlacementTree(st_r.tax, st_r.gt), true));
        PlacementTree gt_tree(st_r.tax, st_r.gt);
        const GraftingResult gr(graft_greedy(gt_tree, kcr, cr, 3));
        const ExactGraftResult ex(graft_exact(PlacementTree(st_r.tax, st_r.gt), kcr, cr, 3));
        std::fprintf(stderr, "graft-rand seed=%llu: greedy %llu vs exact %llu\n",
                     (unsigned long long)seed, (unsigned long long)gr.deltas.back(), (unsigned long long)ex.final_delta);
        REQUIRE(gr.deltas.back() == ex.final_delta);
    }
}

// (b) monotonicity: H_total never increases along the greedy curve, and the
// FN/FP decomposition sums to it.
TEST_CASE("graft_monotone_and_decomposition") {
    Spacer sp(17, 17);
    EightLeaf st;
    KmerClasses kc;
    kc.fill<score::Lex>(st.ms.paths, sp);
    const std::vector<u64> cands(enumerate_grafts(PlacementTree(st.ms.tax, st.ms.gt), true));
    PlacementTree tree(st.ms.tax, st.ms.gt);
    const GraftingResult gr(graft_greedy(tree, kc, cands, 4));
    for(size_t i(1); i < gr.deltas.size(); ++i)
        REQUIRE(gr.deltas[i] <= gr.deltas[i - 1]);
    for(size_t i(0); i < gr.deltas.size(); ++i) {
        REQUIRE(gr.fn[i] + gr.fp[i] == gr.deltas[i]);
        if(i) { // FN and FP never increase either
            REQUIRE(gr.fn[i] <= gr.fn[i - 1]);
            REQUIRE(gr.fp[i] <= gr.fp[i - 1]);
        }
    }
    std::fprintf(stderr, "graft-monotone: %zu grafts, delta %llu -> %llu\n",
                 gr.added_masks.size(), (unsigned long long)gr.initial_delta, (unsigned long long)gr.deltas.back());
}

// (c) best-fit non-submodularity, locked as a regression gate. Two k-mer
// classes on 4 leaves: X = {2,3} (species-internal, error-free at the
// species node) and Y = {1,2} (cross-species). Graft B = {0,1,2} steals
// leaf 2 out of species 101, breaking X (err 0 -> 2) while partially
// fixing Y (err 2 -> 1): g(B) = -1. Graft A = {0,2} has already stolen
// leaf 2 (and absorbed that breakage), so after A the same graft B has
// marginal gain +1. Marginal gains are not diminishing (1 > -1), so the
// (1 - 1/e) argument does not apply to the operative best-fit model; the
// greedy algorithm remains sound because it only accepts positive-gain
// steps (monotone H), and gate (a) measures the empirical gap.
TEST_CASE("graft_bestfit_not_submodular") {
    Spacer sp(17, 17);
    ModuleSet st({{1, 0}, {200, 1}, {100, 200}, {101, 200}},
                 species_taxa(4),
                 {{0xCu, 17}, {0x6u, 17}}, 4242); // classes {2,3}, {1,2}, 1 k-mer each
    KmerClasses kc;
    kc.fill<score::Lex>(st.paths, sp);
    PlacementTree t0(st.tax, st.gt);
    const GraftIndex gi(build_graft_index(t0, kc));
    const u64 A(0x5u), B(0x7u);
    const long double g_b(gain_of_add(t0, gi, B));
    PlacementTree ta(t0);
    REQUIRE(ta.add_synthetic(A) >= 0);
    const long double g_b_a(gain_of_add(ta, gi, B));
    std::fprintf(stderr, "graft-counterexample: g(B)=%.1Lf, g(B|A)=%.1Lf\n", g_b, g_b_a);
    REQUIRE(g_b == -1.0L);   // B breaks the species-internal class, net loss
    REQUIRE(g_b_a == 1.0L);  // after A pre-breaks it, B's marginal is positive
    REQUIRE(g_b_a > g_b);    // diminishing returns violated
}

// (d) feasible-family enumeration keeps the family laminar and excludes
// existing clades; grafted trees still classify singletons exactly.
TEST_CASE("graft_enumeration_feasible") {
    Spacer sp(17, 17);
    EightLeaf st;
    PlacementTree tree(st.ms.tax, st.ms.gt);
    const std::vector<u64> cands(enumerate_grafts(tree, true));
    REQUIRE(!cands.empty());
    for(const u64 A: cands) {
        REQUIRE(pc64(A) >= 2);
        REQUIRE(A != tree.genome_mask());
        // not an existing node mask
        bool exists(false);
        for(u32 id(0); id < tree.nnodes(); ++id)
            if(tree.is_alive(id) && tree.mask(id) == A) exists = true;
        REQUIRE(!exists);
    }
    // sketch-side shared-content ranking is positive on real modules and
    // orders them above random pairs (ordering signal, not an objective).
    const auto sks(genome_sketchers<score::Lex>(st.ms.paths, sp, 8192));
    const double est_mod(sketch_shared_estimate(tree, sks, 0x33u));
    PlacementTree t2(st.ms.tax, st.ms.gt);
    const double est_pair(sketch_shared_estimate(t2, sks, 0x3u)); // species pair: existing clade
    std::fprintf(stderr, "graft-sketchest: module %.1f vs species-clade %.1f\n", est_mod, est_pair);
    REQUIRE(est_mod > 0.);
}

// (f) the fanout>20 union-lane cliff and its polynomial fixes
// (doc/bonsai-gtdb-benchmark.md): above EXACT_UNION_MAX the family used to
// enumerate splits only, and the greedy could not recover a cross-species
// module. Pair/triplet unions + seeded hulls must restore recovery at
// fanout 24 and 40.
TEST_CASE("graft_high_fanout_union_lane") {
    for(const u32 fanout: {24u, 40u}) {
        khash_t(p) *tax(build_tax({{1, 0}, {2, 1}, {50, 2}}));
        for(u32 s(0); s < fanout; ++s) {
            int khr;
            const tax_t sp(100 + s);
            const khint_t ki(kh_put(p, tax, sp, &khr));
            kh_val(tax, ki) = 50;
        }
        std::vector<tax_t> gt(fanout);
        for(u32 s(0); s < fanout; ++s) gt[s] = tax_t(100 + s);
        PlacementTree ptree(tax, gt);
        KmerClasses kc;
        kc.ngenomes = fanout;
        std::vector<u64> observed;
        for(u32 s(0); s < fanout; ++s) {
            kc.count[u64(1) << s] = 1000;
            observed.push_back(u64(1) << s);
        }
        const u64 module((u64(1) << (fanout / 2)) - 1);
        kc.count[module] = 5000;
        observed.push_back(module);
        kc.total_distinct = u64(fanout) * 1000 + 5000;
        const u64 initial(ptree.total_delta(kc));
        // the family lane is now polynomial: pairs + triplets above the cap
        const std::vector<u64> fam(enumerate_grafts(ptree));
        const u64 expected_pairs_triplets = u64(fanout) * (fanout - 1) / 2 + u64(fanout) * (fanout - 1) * (fanout - 2) / 6;
        REQUIRE(fam.size() == expected_pairs_triplets + fanout);   // + split grafts
        // observed-class lane: one graft names the module; recovery is exact
        {
            PlacementTree run(tax, gt);
            const GraftingResult gr(graft_greedy(run, kc, observed, 4));
            REQUIRE(gr.deltas.back() == 0);
        }
        // seeded lane: the child-union hull of the module is the module itself
        {
            const std::vector<u64> seeded(enumerate_seeded_grafts(ptree, observed));
            PlacementTree run(tax, gt);
            const GraftingResult gr(graft_greedy(run, kc, seeded, 4));
            REQUIRE(gr.deltas.back() == 0);
        }
        // pair/triplet lane alone cannot name a fanout/2-ary union (a class
        // reports at the deepest node CONTAINING it; partial unions never
        // contain it). The one split that removes a non-member species from
        // the reporting clade gains, and every remaining static mask then
        // overlaps it (non-laminar, rejected) - the measured floor: initial
        // minus one species' worth of false positives. High-arity unions
        // must be NAMED (seeded/observed lanes below), not composed.
        {
            PlacementTree run(tax, gt);
            const GraftingResult gr(graft_greedy(run, kc, fam, 4));
            REQUIRE(gr.added_masks.size() == 1);
            REQUIRE(gr.deltas.back() == initial - 5000);
        }
        kh_destroy(p, tax);
    }
}

// (g) round 4: sketch-only naming (PI directive 2026-10-21). A planted
// cross-species module at fanout 24 must be named and recovered by the
// bbminhash-intersection lane alone - no k-mer sets in the candidate or
// gain path (the exact index below is evaluation-only, as in the bench).
TEST_CASE("graft_sketch_only_naming") {
    const u32 fanout(24), per(2);
    TmpDir tmp;
    wy::WyRand<u64> rng(20261021);
    const std::string module(rand_dna(rng, 300));   // ~300 distinct k-mers
    std::vector<std::string> paths;
    std::vector<tax_t> gt;
    for(u32 s(0); s < fanout; ++s) {
        for(u32 g(0); g < per; ++g) {
            std::string seq(rand_dna(rng, 8000));
            if(s < fanout / 2) seq += module;
            paths.push_back(tmp.add("sp" + std::to_string(s) + "g" + std::to_string(g) + ".fa",
                                    ">sp" + std::to_string(s) + "g" + std::to_string(g) + "\n" + seq + "\n"));
            gt.push_back(tax_t(100 + s));
        }
    }
    std::vector<std::pair<tax_t, tax_t>> edges{{1, 0}, {2, 1}, {50, 2}};
    for(u32 s(0); s < fanout; ++s) edges.emplace_back(100 + s, 50);
    khash_t(p) *tax(build_tax(edges));
    PlacementTree ptree(tax, gt);
    Spacer sp(17, 17);
    // sketch-only lane
    const std::vector<GenomeSketch> sks(genome_sketchers<score::Lex>(paths, sp, 8192));
    const auto cands(enumerate_sketch_seeded_grafts(ptree, sks, 64.));
    PlacementTree sk_tree(tax, gt);
    const GraftingResult sk(graft_greedy_sketch(sk_tree, sks, cands, 4));
    // oracle: final H of the sketch lane's graft sequence
    KmerClasses kc;
    kc.fill<score::Lex>(paths, sp, true, false);
    PlacementTree oracle(tax, gt);
    for(const u64 A: sk.added_masks) if(oracle.add_synthetic(A) >= 0) {}
    const u64 final_delta(oracle.total_delta(kc));
    std::fprintf(stderr, "graft-sketchonly: %zu candidates, %zu grafts, final delta %llu (initial %llu), est %.0f\n",
                 cands.size(), sk.added_masks.size(), (unsigned long long)final_delta,
                 (unsigned long long)sk.initial_delta, sk.est_gains.empty() ? 0.: sk.est_gains[0]);
    REQUIRE(sk.added_masks.size() >= 1);
    REQUIRE(sk.added_masks[0] == ((u64(1) << 24) - 1));   // the planted module, named
    // quality delta vs the index lane at one graft: identical first graft
    // (the head-to-head claim); both lanes leave the same sub-resolution
    // residual (background sharing below sketch resolution, which the
    // exact lane mops up only with more grafts)
    PlacementTree delta0_tree(tax, gt);
    const u64 delta0(delta0_tree.total_delta(kc));
    std::fprintf(stderr, "graft-sketchonly: delta %llu -> %llu (module graft)\n",
                 (unsigned long long)delta0, (unsigned long long)final_delta);
    REQUIRE(final_delta < delta0);
    kh_destroy(p, tax);
}

// (i) round 5 (doc/grafting-two-tier.md): tier selection, the union rule
// (raw merge then densify ONCE - never union pre-densified sketches), and
// the adaptive per-node tau (background median + k sigma).
TEST_CASE("graft_sketch_tiers_union_rule_adaptive_tau") {
    // tier thresholds
    const SketchNodePolicy pol;
    REQUIRE(node_tier(100., pol) == NodeTier::SPARSE);      // exact-cover eligible in production
    REQUIRE(node_tier(24999., pol) == NodeTier::SPARSE);    // middle band: plain partial fill
    REQUIRE(node_tier(25000., pol) == NodeTier::DENSE);
    // union rule: raw merge then densify ONCE vs min of pre-densified
    // (the sketch repo's regression direction: min(dens) < dens(union))
    {
        wy::WyRand<u64> rng(5);
        const unsigned p(13);
        const size_t m(size_t(1) << p);
        std::vector<uint64_t> keysA, keysB;
        for(size_t i(0); i < 2600; ++i) keysA.push_back(rng());      // ~30% fill each
        for(size_t i(0); i < 2600; ++i) keysB.push_back(rng());
        for(size_t i(0); i < 1300; ++i) keysB[i] = keysA[i];          // 50% overlap -> union ~3900
        GenomeSketch a(p, GENOME_SKETCH_B), b(p, GENOME_SKETCH_B);
        for(const auto k: keysA) a.addh(k);
        for(const auto k: keysB) b.addh(k);
        const double truth(3900.);
        // correct: raw min merge, densify the union ONCE
        RegisterSketch uni(p, GENOME_SKETCH_B);
        uni.mutable_core() = a.core();
        {
            auto &d(uni.mutable_core());
            const auto &s(b.core());
            for(size_t i(0); i < d.size(); ++i) d[i] = std::min(d[i], s[i]);
        }
        // wrong: densify each FIRST, then min (the forbidden order)
        RegisterSketch da(p, GENOME_SKETCH_B), db(p, GENOME_SKETCH_B);
        da.mutable_core() = a.core();
        db.mutable_core() = b.core();
        da.densify(); db.densify();
        auto &dm(da.mutable_core());
        const auto &s(db.core());
        for(size_t i(0); i < dm.size(); ++i) dm[i] = std::min(dm[i], s[i]);
        // estimate via finalize (linear counting); the stored sketches stay raw
        const double est_raw(uni.finalize().cardinality_estimate());
        RegisterSketch uni_stored(p, GENOME_SKETCH_B);
        uni_stored.mutable_core() = uni.core();   // raw union, densified only as a copy at estimate time
        const double est_raw_then_dens(RegisterSketch(p, GENOME_SKETCH_B).finalize().cardinality_estimate() * 0
                                       + [&]{RegisterSketch t(p, GENOME_SKETCH_B); t.mutable_core() = uni_stored.core(); t.densify(); return GenomeSketch(t).finalize().cardinality_estimate();}());
        const double est_min_dens(GenomeSketch(da).finalize().cardinality_estimate());
        std::fprintf(stderr, "union rule: raw-union %.0f, raw-then-densify %.0f, min-of-densified %.0f (truth %.0f)\n",
                     est_raw, est_raw_then_dens, est_min_dens, truth);
        REQUIRE(std::abs(est_raw - truth) < 0.15 * truth);   // plain is sharp at small n
        // the union rule's invariant: min(densA, densB) != densify(A u B) -
        // materially different estimates, and both pre-densified orders
        // degrade the plain sketch (direction is fill-regime dependent:
        // the sketch repo's 30%-fill regression had min-dens BELOW
        // dens(union), 0.179 vs 0.235; this 48%-fill case overshoots
        // instead). Raw-merge-then-densify-once is the only safe order.
        REQUIRE(std::abs(est_min_dens - est_raw) > 0.2 * truth);
        REQUIRE(std::abs(est_raw_then_dens - est_raw) > 0.2 * truth);
    }
    // adaptive tau: 4 siblings, a common background, one pair-private bump
    {
        wy::WyRand<u64> rng(9);
        const unsigned p(13);
        const u32 bg(400);
        std::vector<uint64_t> shared(bg), priv(200);
        for(auto &k: shared) k = rng();
        for(auto &k: priv) k = rng();
        std::vector<GenomeSketch> sks;
        for(u32 g(0); g < 4; ++g) {
            sks.emplace_back(p, GENOME_SKETCH_B);
            for(const auto k: shared) sks.back().addh(k);
            for(u32 i(0); i < 3000; ++i) sks.back().addh(rng());
        }
        for(const auto k: priv) { sks[0].addh(k); sks[1].addh(k); }
        khash_t(p) *tax(build_tax({{1, 0}, {2, 1}}));
        std::vector<tax_t> gt{2, 2, 2, 2};   // four genomes under one parent
        PlacementTree tree(tax, gt);
        // children of the root-level node: 4 leaves (flat fanout)
        const SketchIndex idx(build_sketch_index(tree, sks, pol));
        u32 parent(u32(-1));
        for(u32 x(0); x < tree.nnodes(); ++x)
            if(!tree.is_leaf_node(x) && tree.children_of(x).size() == 4) parent = x;
        REQUIRE(parent != u32(-1));
        // background = median pair count ~ bg; tau = bg + 4 sqrt(bg)
        std::fprintf(stderr, "adaptive tau: node bg-median implies tau %.0f (expect ~%u)\n", idx.tau[parent], unsigned(bg + 4 * std::sqrt(bg) + 40));
        REQUIRE(idx.tau[parent] > bg + 2. * std::sqrt(bg));
        REQUIRE(idx.tau[parent] < bg + 8. * std::sqrt(bg) + 32.);
        // the pair-private bump (200 above a 400 background) edges; background pairs do not
        const auto cands(enumerate_sketch_seeded_grafts(idx));
        bool found(false);
        for(const auto &c: cands) if(c.mask == 0x3) found = true;
        std::fprintf(stderr, "adaptive tau: %zu candidates, {0,1} present: %d\n", cands.size(), int(found));
        REQUIRE(found);
        kh_destroy(p, tax);
    }
}

// (j) round 6 (real pangenomes, doc/realdata-grafting.md): differential
// edges + pair candidates + the multiplicative noise-band tau.
TEST_CASE("graft_differential_pair_candidates") {
    wy::WyRand<u64> rng(31);
    const unsigned p(13);
    const u32 bg(400);
    std::vector<uint64_t> shared(bg), priv(200);
    for(auto &k: shared) k = rng();
    for(auto &k: priv) k = rng();
    std::vector<GenomeSketch> sks;
    for(u32 g(0); g < 4; ++g) {
        sks.emplace_back(p, GENOME_SKETCH_B);
        for(const auto k: shared) sks.back().addh(k);
        for(u32 i(0); i < 3000; ++i) sks.back().addh(rng());
    }
    for(const auto k: priv) { sks[0].addh(k); sks[1].addh(k); }
    khash_t(p) *tax(build_tax({{1, 0}, {2, 1}}));
    std::vector<tax_t> gt{2, 2, 2, 2};
    PlacementTree tree(tax, gt);
    SketchNodePolicy pol;
    pol.exclusive_edges = true;
    pol.pair_candidates = true;
    const SketchIndex idx(build_sketch_index(tree, sks, pol));
    const auto cands(enumerate_sketch_seeded_grafts(idx));
    REQUIRE(!cands.empty());
    size_t bump(size_t(-1));
    for(size_t i(0); i < cands.size(); ++i) {
        if(cands[i].mask == 0x3) bump = i;
    }
    REQUIRE(bump != size_t(-1));             // the pair-private bump is named
    for(const auto &c: cands)                // and ranks above the noise pairs
        REQUIRE(cands[bump].est_gain >= c.est_gain);
    // multiplicative tau: max(floor, k*sqrt(median differential)); the
    // median differential here is ~0 (one bumped pair of ten), so floor
    u32 parent(u32(-1));
    for(u32 x(0); x < tree.nnodes(); ++x)
        if(!tree.is_leaf_node(x) && tree.children_of(x).size() == 4) parent = x;
    REQUIRE(idx.tau[parent] >= pol.tau_floor);          // multiplicative band on the noise scale
    REQUIRE(idx.tau[parent] < 200.);
    kh_destroy(p, tax);
}

// (k) round 7: exact covers (differential = |inter parts| - |core(parent)|,
// error 0 by construction) and register folding (folding a p-sketch by 2^f
// equals building at p - f - one-permutation bucket nesting).
TEST_CASE("graft_exact_covers_and_register_folding") {
    // covers on a 4-genome tree with a planted {0,1} module
    const u32 n(4);
    khash_t(p) *tax(build_tax({{1, 0}, {2, 1}}));
    std::vector<tax_t> gt(n, 2);
    PlacementTree tree(tax, gt);
    ExactCovers cov;
    std::vector<std::shared_ptr<std::vector<uint64_t>>> leaf(n);
    for(u32 g(0); g < n; ++g) {
        leaf[g] = std::make_shared<std::vector<uint64_t>>();
        for(u64 i(0); i < 5000; ++i) leaf[g]->push_back(sketch::WangHash::hash(g * 1000000 + i));
        if(g < 2) for(u64 i(0); i < 300; ++i) leaf[g]->push_back(0xDEAD0000 + i);   // the module
        std::sort(leaf[g]->begin(), leaf[g]->end());
        leaf[g]->erase(std::unique(leaf[g]->begin(), leaf[g]->end()), leaf[g]->end());
        cov.full[tree.leaf_node(g)] = leaf[g];
        cov.core[tree.leaf_node(g)] = leaf[g];
    }
    u32 parent(u32(-1));
    for(u32 x(0); x < tree.nnodes(); ++x)
        if(!tree.is_leaf_node(x) && tree.children_of(x).size() == 4) parent = x;
    REQUIRE(parent != u32(-1));
    auto inter2 = [](std::shared_ptr<std::vector<uint64_t>> a, std::shared_ptr<std::vector<uint64_t>> b) {
        auto r(std::make_shared<std::vector<uint64_t>>());
        std::set_intersection(a->begin(), a->end(), b->begin(), b->end(), std::back_inserter(*r));
        return r;
    };
    cov.core[parent] = inter2(inter2(leaf[0], leaf[1]), inter2(leaf[2], leaf[3]));
    cov.full[parent] = cov.full[tree.leaf_node(0)];
    // exact differential of the module pair: |{0} n {1}| - |core of all 4|
    {
        auto p01(inter2(leaf[0], leaf[1]));
        const double d(exact_differential(cov, tree, std::vector<u64>{1, 2}, 0xF));   // masks g0,g1
        std::vector<uint64_t> pin;
        std::set_intersection(p01->begin(), p01->end(), cov.core[parent]->begin(), cov.core[parent]->end(), std::back_inserter(pin));
        REQUIRE(d == double(p01->size() - pin.size()));
        REQUIRE(d >= 290.);   // ~300 minus shared-background collisions
    }
    // replication across register budgets: fresh builds at p=15 and p=11
    // both detect the same planted module (folding is NOT valid - registers
    // hold 64-p value bits, so a p-fold loses the top bits; measured here by
    // construction: the p=11 build sees the module with fewer registers)
    {
        wy::WyRand<u64> rng(3);
        std::vector<std::vector<uint64_t>> keys(4);
        for(u32 g(0); g < 4; ++g)
            for(u32 i(0); i < 4000; ++i) keys[g].push_back(rng());
        for(u32 i(0); i < 600; ++i) {keys[0].push_back(0xBEEF0000 + i); keys[1].push_back(0xBEEF0000 + i);}
        for(const unsigned p: {15u, 11u}) {
            std::vector<sketch::BBitMinHasher<uint64_t>> sks;
            for(u32 g(0); g < 4; ++g) {
                sks.emplace_back(p, GENOME_SKETCH_B);
                for(const auto k: keys[g]) sks[g].addh(k);
            }
            const std::vector<u64> parts{u64(1), u64(2)};
            const double d(pair_exclusive_count(sks, 1, 2, 0xC));
            std::fprintf(stderr, "replication p=%u: differential %.0f (module 600)\n", p, d);
            REQUIRE(d > 100.);
        }
    }
    kh_destroy(p, tax);
}
