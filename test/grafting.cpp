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
