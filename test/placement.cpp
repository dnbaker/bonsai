#include "test/catch.hpp"

#include "placement.h"
#include "tx.h"
#include "util.h"

#include <map>
#include <random>
#include <string>
#include <unordered_set>
#include <vector>

using namespace bns;

namespace {

khash_t(p) *build_taxmap(std::initializer_list<std::pair<tax_t, tax_t>> pairs) {
    khash_t(p) *ret(kh_init(p));
    int khr;
    for(const auto &pair: pairs) {
        const auto ki(kh_put(p, ret, pair.first, &khr));
        kh_val(ret, ki) = pair.second;
    }
    return ret;
}

// Two-genera / 16-species taxonomy over 32 genomes (genome 2i, 2i+1 -> species 100+i).
khash_t(p) *build_synth_taxonomy() {
    khash_t(p) *ret(kh_init(p));
    int khr;
    auto put = [&](tax_t k, tax_t v) {
        const khint_t ki(kh_put(p, ret, k, &khr)); // two-step: kh_val(h, kh_put(...)) = v miscompiles under g++-15 -O1 -march=native
        kh_val(ret, ki) = v;
    };
    put(1, 0); put(200, 1); put(201, 1);
    for(unsigned i(0); i < 8; ++i) put(100 + i, 200);
    for(unsigned i(8); i < 16; ++i) put(100 + i, 201);
    return ret;
}

std::string rand_dna(std::mt19937_64 &rng, size_t len) {
    static const char bases[] = "ACGT";
    std::string ret;
    ret.reserve(len);
    while(ret.size() < len) ret.push_back(bases[(rng() >> 33) & 3]); // mt19937_64 low bits are low-quality; use high bits
    return ret;
}

struct TmpDir {
    std::string dir;
    std::vector<std::string> files;
    TmpDir(): dir([] {
        char tmpl[] = "/tmp/bns_place_XXXXXX";
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

void require_monotone(const GreedyResult &gr) {
    for(size_t i(1); i < gr.curve.size(); ++i)
        REQUIRE(gr.curve[i].unweighted <= gr.curve[i - 1].unweighted + 1e-12);
}

} // namespace

// Hand-checkable acceptance gate for the metric and the add semantics.
TEST_CASE("placement_metric_exact") {
    khash_t(p) *tax(build_taxmap({{1, 0}, {100, 1}, {200, 1}}));
    const std::vector<tax_t> gt{100, 100, 200}; // g0, g1 in species 100; g2 in species 200
    PlacementTree tree(tax, gt);
    REQUIRE(tree.implied(0b011) == 0b011); // species block classifies exactly
    REQUIRE(tree.implied(0b101) == 0b111); // cross-species pair lifts to the root
    REQUIRE(tree.err(0b101) == 1);
    REQUIRE(tree.err(0b001) == 0);

    KmerClasses kc;
    kc.ngenomes = 3;
    kc.total_distinct = 11;
    kc.count = {{0b011, 3}, {0b101, 1}, {0b001, 2}, {0b100, 5}};
    kc.weight = {{0b011, 3}, {0b101, 1}, {0b001, 2}, {0b100, 5}};
    REQUIRE(tree.total_delta(kc) == 1);
    REQUIRE(std::abs(tree.unweighted_error(kc) - 1. / 11.) < 1e-12);

    // Adding {g0,g2} fixes the one {g0,g2} k-mer but scatters the {g0,g1}
    // class to the root: net gain 0, so the greedy must decline it.
    PlacementTree copy(tree);
    const GreedyResult gr(greedy_synthetic_nodes(copy, kc, 2));
    REQUIRE(gr.added.empty());
    REQUIRE(copy.total_delta(kc) == 1);

    // Without the competing species class, the same node is a clean win.
    KmerClasses kc2;
    kc2.ngenomes = 3;
    kc2.total_distinct = 4;
    kc2.count = {{0b101, 4}};
    kc2.weight = {{0b101, 4}};
    const GreedyResult gr2(greedy_synthetic_nodes(tree, kc2, 2));
    REQUIRE(gr2.added.size() == 1);
    REQUIRE(gr2.added_masks.front() == 0b101);
    REQUIRE(tree.total_delta(kc2) == 0);
    // g0 and g2 moved; species 100 shrank to {g1}.
    REQUIRE(tree.implied(0b001) == 0b001);
    REQUIRE(tree.implied(0b010) == 0b010);
    REQUIRE(tree.err(0b011) == 1); // species class degraded by the theft

    // The full genome set is already the root: no node to add.
    PlacementTree t4(tax, gt);
    REQUIRE(t4.add_synthetic(0b111) < 0);

    // Adoption-only addition (whole species subtrees move) is exactly undone
    // by removal; masks return to their pre-add values.
    khash_t(p) *tax4(build_taxmap({{1, 0}, {100, 1}, {200, 1}, {300, 1}}));
    PlacementTree t5(tax4, std::vector<tax_t>{100, 100, 200, 300});
    const u64 pre_mask(t5.implied(0b0111));
    const i64 id(t5.add_synthetic(0b0111)); // {g0,g1,g2}: adopts species 100 and 200
    REQUIRE(id >= 0);
    REQUIRE(t5.implied(0b0111) == 0b0111);
    REQUIRE(t5.remove_synthetic(u32(id)));
    REQUIRE(t5.implied(0b0111) == pre_mask);
    REQUIRE(t5.implied(0b0011) == 0b0011);
    REQUIRE(t5.implied(0b0100) == 0b0100);
    kh_destroy(p, tax4);
    kh_destroy(p, tax);
}

/*
Structured case: 16 species x 2 genomes under 2 genera. Three modules of
8/4/2 species (16/8/4 genomes) share k-mer blocks and cut across genera;
the big module's block is repeated 3x per genome to separate the weighted
curve. Species pairs share species blocks (already exact in the taxonomy).
Expected initial total error: 50 * ((32-16) + (32-8) + (32-4)) = 3400,
removed by 3 synthetic nodes. (60bp blocks give 44 exact kmers at k=17.)
*/
TEST_CASE("greedy_curve_structured") {
    TmpDir tmp;
    std::mt19937_64 rng(1337);
    enum {NSPECIES = 16, NGEN = 32};
    khash_t(p) *tax(build_synth_taxonomy());
    auto in_m16 = [](unsigned sp) {return sp < 4 || (sp >= 8 && sp < 12);};
    std::string species_blocks[NSPECIES];
    for(unsigned i(0); i < NSPECIES; ++i) species_blocks[i] = rand_dna(rng, 60);
    const std::string m16_block(rand_dna(rng, 60)), m8_block(rand_dna(rng, 60)), m4_block(rand_dna(rng, 60));
    std::vector<std::string> paths;
    std::vector<tax_t> gt;
    for(unsigned g(0); g < NGEN; ++g) {
        const unsigned sp(g >> 1);
        std::string seq(rand_dna(rng, 1300));
        seq += species_blocks[sp];
        if(in_m16(sp)) seq += m16_block + rand_dna(rng, 17) + m16_block + rand_dna(rng, 17) + m16_block;
        else if(sp == 4 || sp == 5 || sp == 12 || sp == 13) seq += m8_block;
        else if(sp == 6 || sp == 14) seq += m4_block;
        const std::string name("g" + std::to_string(g));
        paths.push_back(tmp.add(name + ".fa", ">" + name + "\n" + seq + "\n"));
        gt.push_back(100 + sp);
    }
    PlacementTree tree(tax, gt);
    Spacer spacer(17, 17);
    KmerClasses kc;
    kc.fill<score::Lex>(paths, spacer);
    std::fprintf(stderr, "structured: initial unweighted %g weighted %g delta %llu over %zu distinct kmers\n",
                 tree.unweighted_error(kc), tree.weighted_error(kc), (unsigned long long)tree.total_delta(kc), (size_t)kc.total_distinct);
    // The 16-bit boundary where repeated module copies meet their random
    // separators yields boundary 17-mers (16 fixed block bases + 1 random
    // base) that collide across genomes 1/4 of the time. Those classes are
    // subsets of the M16 mask and are partially absorbed by its node, but a
    // junk floor remains; assertions below account for it.
    const u64 m16_mask(0x0000000000ff00ff), m8_mask(0x000000000f000f00), m4_mask(0x0000000030003000);
    const u64 initial_delta(tree.total_delta(kc));
    REQUIRE(initial_delta >= 44 * (16 + 24 + 28));      // 2992 of designed module error at k=17
    REQUIRE(initial_delta < 2992 + 4000);               // plus bounded boundary junk

    PlacementTree greedy_tree(tree);
    const GreedyResult gr(greedy_synthetic_nodes(greedy_tree, kc, 10));
    require_monotone(gr);
    for(const auto &p: gr.curve)
        std::fprintf(stderr, "structured: budget %u unweighted %.9g weighted %.9g delta %llu\n",
                     p.budget, p.unweighted, p.weighted, (unsigned long long)p.total_delta);
    // Exactly the three module nodes are added, and error falls by >2/3.
    for(auto m: gr.added_masks) std::fprintf(stderr, "structured: added node mask=%016llx pc=%u\n", (unsigned long long)m, pc64(m));
    REQUIRE(gr.added.size() >= 3);
    {
        const std::set<u64> have(gr.added_masks.begin(), gr.added_masks.end());
        for(auto m: {m16_mask, m8_mask, m4_mask}) {INFO("module node missing"); REQUIRE(have.count(m) == 1);}
    }
    REQUIRE(gr.curve.back().total_delta * 3 <= initial_delta);
    for(size_t i(1); i < gr.curve.size(); ++i)
        REQUIRE(gr.curve[i].weighted <= gr.curve[i - 1].weighted + 1e-12);

    // Surplus beyond the greedy set: the only remaining candidates are the
    // species bitmaps, which are already nodes, so trim must be a no-op.
    const double pre_trim(greedy_tree.unweighted_error(kc));
    const TrimResult tr(trim_synthetic_nodes(greedy_tree, kc, gr.added, 0.));
    REQUIRE(tr.removed == 0);
    REQUIRE(tr.kept == gr.added.size());
    REQUIRE(tr.final_unweighted == pre_trim);

    // Removing one adopted module node reopens error, but stays below the
    // initial error since the other two modules remain fixed.
    const u64 pre_remove(greedy_tree.total_delta(kc));
    REQUIRE(greedy_tree.remove_synthetic(gr.added.front()));
    const u64 post_remove(greedy_tree.total_delta(kc));
    REQUIRE(post_remove > pre_remove);
    REQUIRE(post_remove < initial_delta);
    kh_destroy(p, tax);
}

// Null control 1: pure random genomes have no shared k-mers.
TEST_CASE("greedy_null_random") {
    TmpDir tmp;
    std::mt19937_64 rng(99);
    enum {NGEN = 32};
    khash_t(p) *tax(build_synth_taxonomy());
    std::vector<std::string> paths;
    std::vector<tax_t> gt;
    for(unsigned g(0); g < NGEN; ++g) {
        const std::string name("g" + std::to_string(g));
        paths.push_back(tmp.add(name + ".fa", ">" + name + "\n" + rand_dna(rng, 1600) + "\n"));
        gt.push_back(100 + (g >> 1));
    }
    PlacementTree tree(tax, gt);
    Spacer spacer(17, 17);
    KmerClasses kc;
    kc.fill<score::Lex>(paths, spacer);
    std::fprintf(stderr, "null-random: initial unweighted %g delta %llu classes %zu\n",
                 tree.unweighted_error(kc), (unsigned long long)tree.total_delta(kc), kc.count.size());
    REQUIRE(tree.unweighted_error(kc) == 0.);
    const GreedyResult gr(greedy_synthetic_nodes(tree, kc, 10));
    REQUIRE(gr.added.empty());
    kh_destroy(p, tax);
}

/*
Null control 2: cross-genus genome pairs share short blocks. Representing
a pair requires stealing leaves and destroying shared structure, so the
greedy plateaus quickly and the residual stays high.
*/
TEST_CASE("greedy_null_pairs") {
    TmpDir tmp;
    std::mt19937_64 rng(2024);
    enum {NGEN = 32, NPAIRS = 40};
    khash_t(p) *tax(build_synth_taxonomy());
    std::vector<std::string> paths;
    std::vector<tax_t> gt;
    for(unsigned g(0); g < NGEN; ++g) {
        const std::string name("g" + std::to_string(g));
        paths.push_back(tmp.add(name + ".fa", ">" + name + "\n" + rand_dna(rng, 1500) + "\n"));
        gt.push_back(100 + (g >> 1));
    }
    std::vector<std::pair<unsigned, unsigned>> pairs;
    for(unsigned p(0); p < NPAIRS; ++p) {
        const unsigned i(rng() % NGEN), j(rng() % NGEN);
        if(i == j || (i >> 1) == (j >> 1) || (i >> 4) == (j >> 4)) {--p; continue;} // cross-species, cross-genus
        pairs.emplace_back(i, j);
        const std::string block(rand_dna(rng, 30)); // 14 exact 17-mers per copy
        for(auto path: {paths[i], paths[j]}) {
            std::FILE *fp(std::fopen(path.data(), "a"));
            std::fputs((">" + std::to_string(p) + "_block\n" + block + "\n").data(), fp);
            std::fclose(fp);
        }
    }
    PlacementTree tree(tax, gt);
    Spacer spacer(17, 17);
    KmerClasses kc;
    kc.fill<score::Lex>(paths, spacer);
    const u64 initial_delta(tree.total_delta(kc));
    std::fprintf(stderr, "null-pairs: %zu pairs, initial unweighted %g delta %llu over %zu distinct\n",
                 pairs.size(), tree.unweighted_error(kc), (unsigned long long)initial_delta, (size_t)kc.total_distinct);
    REQUIRE(initial_delta == pairs.size() * 14 * 30);
    const GreedyResult gr(greedy_synthetic_nodes(tree, kc, 64));
    require_monotone(gr);
    for(const auto &p: gr.curve)
        std::fprintf(stderr, "null-pairs: budget %u unweighted %.9g delta %llu\n", p.budget, p.unweighted, (unsigned long long)p.total_delta);
    REQUIRE(gr.added.size() >= 1);
    REQUIRE(gr.added.size() < pairs.size());
    REQUIRE(gr.curve.back().total_delta > initial_delta / 3);

    // Surplus: force every remaining candidate in, then trim at zero tolerance.
    PlacementTree surplus_tree(tree);
    std::vector<u32> surplus(gr.added);
    std::unordered_set<u64> forced(gr.added_masks.begin(), gr.added_masks.end());
    for(const auto &pair: kc.count) {
        if(pc64(pair.first) < 2 || forced.find(pair.first) != forced.end()) continue;
        const i64 id(surplus_tree.add_synthetic(pair.first));
        if(id >= 0) {surplus.push_back(u32(id)); forced.insert(pair.first);}
    }
    const double pre_trim(surplus_tree.unweighted_error(kc));
    const TrimResult tr(trim_synthetic_nodes(surplus_tree, kc, surplus, 0.));
    std::fprintf(stderr, "null-pairs: surplus %zu pre-trim %g removed %u kept %u final %g (greedy endpoint %g)\n",
                 surplus.size(), pre_trim, tr.removed, tr.kept, tr.final_unweighted, gr.curve.back().unweighted);
    REQUIRE(tr.removed >= 1);
    REQUIRE(tr.final_unweighted <= pre_trim + 1e-12);
    kh_destroy(p, tax);
}

// End-to-end with round 1's TaxonomyReformation as the tree source.
TEST_CASE("placement_reformation_compose") {
    TmpDir tmp;
    std::mt19937_64 rng(555);
    enum {NSPECIES = 4, NGEN = 8};
    const std::string taxpath(tmp.add("tax.txt", "100|\t1|\n101|\t1|\n102|\t1|\n103|\t1|\n"));
    khash_t(p) *taxmap(build_parent_map(taxpath.data()));
    std::string names;
    std::vector<std::string> paths;
    const std::string species_blocks[NSPECIES] = {rand_dna(rng, 50), rand_dna(rng, 50), rand_dna(rng, 50), rand_dna(rng, 50)};
    const std::string module_block(rand_dna(rng, 50));
    for(unsigned g(0); g < NGEN; ++g) {
        const unsigned sp(g >> 1);
        const std::string name("g" + std::to_string(g));
        std::string seq(rand_dna(rng, 400) + species_blocks[sp]);
        if(sp == 0 || sp == 2) seq += module_block; // module = species 100 and 102
        const std::string path(tmp.add(name + ".fa", ">" + name + "\n" + seq + "\n"));
        paths.push_back(path);
        names += name + "\t" + std::to_string(100 + sp) + "\n";
    }
    const std::string namepath(tmp.add("names.txt", names));
    TaxonomyReformation tr(namepath.data(), paths, taxmap);
    const std::string nm(tmp.dir + "/nm.gz");
    tr.write_name_map(nm.data());
    gzFile fp(gzopen(nm.data(), "rb"));
    REQUIRE(fp != nullptr);
    std::map<std::string, tax_t> path_to_new;
    char buf[2048];
    while(gzgets(fp, buf, sizeof(buf))) {
        char *tab(std::strchr(buf, '\t'));
        REQUIRE(tab != nullptr);
        *tab = '\0';
        path_to_new[buf] = std::atoi(tab + 1);
    }
    gzclose(fp);
    std::vector<tax_t> new_gt;
    for(const auto &path: paths) new_gt.push_back(path_to_new.at(path));
    PlacementTree tree(tr.pmap(), new_gt);
    Spacer spacer(17, 17);
    KmerClasses kc;
    kc.fill<score::Lex>(paths, spacer);
    std::fprintf(stderr, "reformation-compose: initial delta %llu\n", (unsigned long long)tree.total_delta(kc));
    const GreedyResult gr(greedy_synthetic_nodes(tree, kc, 4));
    for(const auto &p: gr.curve)
        std::fprintf(stderr, "reformation-compose: budget %u unweighted %.9g delta %llu\n", p.budget, p.unweighted, (unsigned long long)p.total_delta);
    REQUIRE(gr.added.size() == 1);
    REQUIRE(gr.curve.back().total_delta == 0);
    kh_destroy(p, taxmap);
}
