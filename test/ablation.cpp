#include "test/catch.hpp"

#include "sketch_select.h"
#include "util.h"

#include <random>
#include <string>
#include <vector>

using namespace bns;

namespace {

std::string rand_dna(std::mt19937_64 &rng, size_t len) {
    static const char bases[] = "ACGT";
    std::string ret;
    ret.reserve(len);
    while(ret.size() < len) ret.push_back(bases[(rng() >> 33) & 3]);
    return ret;
}

struct TmpDir {
    std::string dir;
    std::vector<std::string> files;
    TmpDir(): dir([] {
        char tmpl[] = "/tmp/bns_ablate_XXXXXX";
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

khash_t(p) *build_synth_taxonomy() {
    khash_t(p) *ret(kh_init(p));
    int khr;
    auto put = [&](tax_t k, tax_t v) {
        const khint_t ki(kh_put(p, ret, k, &khr));
        kh_val(ret, ki) = v;
    };
    put(1, 0); put(200, 1); put(201, 1);
    for(unsigned i(0); i < 8; ++i) put(100 + i, 200);
    for(unsigned i(8); i < 16; ++i) put(100 + i, 201);
    return ret;
}

// 16 species x 2 genomes under 2 genera, NO species-level sharing, and two
// modules that cross taxonomic boundaries at genome level: one genome from
// each of species 0-3 (module A) and 8-11... (module B uses species 4-7's
// second genomes), so neither module is a union of taxa.
struct CrossingSet {
    TmpDir tmp;
    std::vector<std::string> paths;
    std::vector<tax_t> gt;
    khash_t(p) *tax;
    const u64 ma_mask = 0x0000000000000055; // genomes 0,2,4,6 (species 0-3, genus 200)
    const u64 mb_mask = 0x0000000000AA0000; // genomes 17,19,21,23 (species 8-11, genus 201)
    CrossingSet(): tax(build_synth_taxonomy()) {
        std::mt19937_64 rng(4242);
        const std::string ma_block(rand_dna(rng, 60)), mb_block(rand_dna(rng, 60));
        for(unsigned g(0); g < 32; ++g) {
            std::string seq(rand_dna(rng, 1400));
            if(g % 2 == 0 && g < 8) seq += ma_block;
            if(g % 2 == 1 && g >= 16 && g < 24) seq += mb_block;
            const std::string name("g" + std::to_string(g));
            paths.push_back(tmp.add(name + ".fa", ">" + name + "\n" + seq + "\n"));
            gt.push_back(100 + (g >> 1));
        }
    }
    ~CrossingSet() {kh_destroy(p, tax);}
};

} // namespace

TEST_CASE("taxonomy_crossing_ablation") {
    CrossingSet st;
    Spacer sp(17, 17);
    KmerClasses kc;
    kc.fill<score::Lex>(st.paths, sp);
    const u64 kUnremovable(0);

    SECTION("unrestricted grafting removes all module error") {
        PlacementTree tree(st.tax, st.gt);
        const u64 initial(tree.total_delta(kc));
        std::fprintf(stderr, "ablation-crossing: initial delta %llu\n", (unsigned long long)initial);
        // Each module lifts only to its genus (16 genomes), so the error is
        // ~44 k-mers x (16 - 4) per module, plus rare block/random collisions.
        REQUIRE(initial >= 2 * 44 * 11);
        REQUIRE(initial <= 2 * 44 * 13);
        const GreedyResult gr(greedy_synthetic_nodes(tree, kc, 8));
        for(const auto &p: gr.curve)
            std::fprintf(stderr, "ablation unrestricted: budget %u delta %llu\n", p.budget, (unsigned long long)p.total_delta);
        const std::set<u64> have(gr.added_masks.begin(), gr.added_masks.end());
        REQUIRE(have.count(st.ma_mask) == 1);
        REQUIRE(have.count(st.mb_mask) == 1);
        REQUIRE(gr.curve.back().total_delta <= initial / 6);
        (void)kUnremovable;
    }
    SECTION("taxonomy-respecting grafting cannot represent the modules") {
        PlacementTree tree(st.tax, st.gt);
        const u64 initial(tree.total_delta(kc));
        std::vector<u64> cands;
        for(const auto &pair: kc.count) if(pc64(pair.first) >= 2) cands.push_back(pair.first);
        const GreedyResult gr(greedy_synthetic_nodes(tree, kc, 8, cands, GraftPolicy::respect_taxonomy));
        std::fprintf(stderr, "ablation respect_taxonomy: added %zu, residual %llu of %llu (cost %llu)\n",
                     gr.added.size(), (unsigned long long)gr.curve.back().total_delta, (unsigned long long)initial, (unsigned long long)initial);
        REQUIRE(gr.added.empty());
        REQUIRE(gr.curve.back().total_delta == initial);
    }
}

// The round-2 structured set has modules that ARE unions of species, so
// respecting the taxonomy costs nothing there.
TEST_CASE("taxonomy_respecting_structured") {
    TmpDir tmp;
    std::mt19937_64 rng(1337);
    khash_t(p) *tax(build_synth_taxonomy());
    std::string species_blocks[16];
    for(unsigned i(0); i < 16; ++i) species_blocks[i] = rand_dna(rng, 60);
    const std::string m16_block(rand_dna(rng, 60)), m8_block(rand_dna(rng, 60)), m4_block(rand_dna(rng, 60));
    auto in_m16 = [](unsigned sp) {return sp < 4 || (sp >= 8 && sp < 12);};
    std::vector<std::string> paths;
    std::vector<tax_t> gt;
    for(unsigned g(0); g < 32; ++g) {
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
    Spacer sp(17, 17);
    KmerClasses kc;
    kc.fill<score::Lex>(paths, sp);
    PlacementTree t1(tax, gt);
    const GreedyResult free_run(greedy_synthetic_nodes(t1, kc, 10));
    PlacementTree t2(tax, gt);
    std::vector<u64> cands;
    for(const auto &pair: kc.count) if(pc64(pair.first) >= 2) cands.push_back(pair.first);
    const GreedyResult respect_run(greedy_synthetic_nodes(t2, kc, 10, cands, GraftPolicy::respect_taxonomy));
    std::fprintf(stderr, "ablation structured: unrestricted %zu nodes final %llu; respect_taxonomy %zu nodes final %llu\n",
                 free_run.added.size(), (unsigned long long)free_run.curve.back().total_delta,
                 respect_run.added.size(), (unsigned long long)respect_run.curve.back().total_delta);
    REQUIRE(respect_run.added_masks == free_run.added_masks); // unions of species: no cost
    REQUIRE(respect_run.curve.back().total_delta == free_run.curve.back().total_delta);
    kh_destroy(p, tax);
}
