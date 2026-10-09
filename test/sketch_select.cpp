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
        char tmpl[] = "/tmp/bns_sketch_XXXXXX";
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

void require_monotone(const GreedyResult &gr) {
    for(size_t i(1); i < gr.curve.size(); ++i)
        REQUIRE(gr.curve[i].unweighted <= gr.curve[i - 1].unweighted + 1e-12);
}

// The round-2 structured set: 16 species x 2 genomes, 2 genera, module
// blocks of 16/8/4 genomes crossing genera (big module repeated 3x).
struct StructuredSet {
    TmpDir tmp;
    std::vector<std::string> paths;
    std::vector<tax_t> gt;
    khash_t(p) *tax;
    const u64 m16_mask = 0x0000000000ff00ff, m8_mask = 0x000000000f000f00, m4_mask = 0x0000000030003000;
    StructuredSet(): tax(build_synth_taxonomy()) {
        std::mt19937_64 rng(1337);
        std::string species_blocks[16];
        for(unsigned i(0); i < 16; ++i) species_blocks[i] = rand_dna(rng, 60);
        const std::string m16_block(rand_dna(rng, 60)), m8_block(rand_dna(rng, 60)), m4_block(rand_dna(rng, 60));
        auto in_m16 = [](unsigned sp) {return sp < 4 || (sp >= 8 && sp < 12);};
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
    }
    ~StructuredSet() {kh_destroy(p, tax);}
};

} // namespace

// (a) MinHash/Jaccard sketch clusters as candidate generators.
TEST_CASE("sketch_candidates_structured") {
    StructuredSet st;
    Spacer sp(17, 17);
    KmerClasses kc;
    kc.fill<score::Lex>(st.paths, sp);
    PlacementTree exact_tree(st.tax, st.gt);
    const u64 initial_delta(exact_tree.total_delta(kc));

    // Module edges sit at J ~ 0.0154 (44 shared of ~2.9k union); species
    // edges at the same 0.0154; random genomes at ~1e-4. 0.008 separates
    // them at 4+ sigma with 4096 bottom-k hashes.
    const auto sketches(genome_sketches<score::Lex>(st.paths, sp, 4096));
    const auto cands(sketch_cluster_candidates(sketches, 0.008));
    std::fprintf(stderr, "sketch-structured: %zu sketch candidates, exact initial delta %llu\n", cands.size(), (unsigned long long)initial_delta);
    unsigned nmodules(0);
    for(const auto m: cands)
        nmodules += (m == st.m16_mask || m == st.m8_mask || m == st.m4_mask);
    REQUIRE(nmodules == 3); // the module masks are recovered

    PlacementTree sketch_tree(st.tax, st.gt);
    const GreedyResult gr(greedy_synthetic_nodes(sketch_tree, kc, 10, cands));
    require_monotone(gr);
    for(const auto &p: gr.curve)
        std::fprintf(stderr, "sketch-structured: budget %u unweighted %.9g delta %llu\n", p.budget, p.unweighted, (unsigned long long)p.total_delta);
    REQUIRE(gr.curve.back().total_delta * 4 <= initial_delta); // knee reproduced
    REQUIRE(std::set<u64>(gr.added_masks.begin(), gr.added_masks.end()) == std::set<u64>{st.m16_mask, st.m8_mask, st.m4_mask});
}

TEST_CASE("sketch_candidates_nulls") {
    Spacer sp(17, 17);
    SECTION("pure random: nothing clusters") {
        TmpDir tmp;
        std::mt19937_64 rng(99);
        std::vector<std::string> paths;
        std::vector<tax_t> gt;
        for(unsigned g(0); g < 32; ++g) {
            const std::string name("g" + std::to_string(g));
            paths.push_back(tmp.add(name + ".fa", ">" + name + "\n" + rand_dna(rng, 1600) + "\n"));
            gt.push_back(100 + (g >> 1));
        }
        khash_t(p) *tax(build_synth_taxonomy());
        const auto sketches(genome_sketches<score::Lex>(paths, sp, 2048));
        const auto cands(sketch_cluster_candidates(sketches, 0.005)); // even a low bar
        std::fprintf(stderr, "sketch-null-random: %zu candidates\n", cands.size());
        REQUIRE(cands.empty());
        PlacementTree tree(tax, gt);
        KmerClasses kc;
        kc.fill<score::Lex>(paths, sp);
        const GreedyResult gr(greedy_synthetic_nodes(tree, kc, 10, cands));
        REQUIRE(gr.added.empty());
        kh_destroy(p, tax);
    }
    SECTION("random cross-genus pairs: conservative at 0.02, plateau at 0.005") {
        TmpDir tmp;
        std::mt19937_64 rng(2024);
        std::vector<std::string> paths;
        std::vector<tax_t> gt;
        for(unsigned g(0); g < 32; ++g) {
            const std::string name("g" + std::to_string(g));
            paths.push_back(tmp.add(name + ".fa", ">" + name + "\n" + rand_dna(rng, 1500) + "\n"));
            gt.push_back(100 + (g >> 1));
        }
        for(unsigned p(0); p < 40; ++p) {
            const unsigned i(rng() % 32), j(rng() % 32);
            if(i == j || (i >> 1) == (j >> 1) || (i >> 4) == (j >> 4)) {--p; continue;}
            const std::string block(rand_dna(rng, 30));
            for(auto path: {paths[i], paths[j]}) {
                std::FILE *fp(std::fopen(path.data(), "a"));
                std::fputs((">" + std::to_string(p) + "_block\n" + block + "\n").data(), fp);
                std::fclose(fp);
            }
        }
        khash_t(p) *tax(build_synth_taxonomy());
        PlacementTree tree(tax, gt);
        KmerClasses kc;
        kc.fill<score::Lex>(paths, sp);
        const u64 initial_delta(tree.total_delta(kc));
        const auto sketches(genome_sketches<score::Lex>(paths, sp, 2048));
        for(const double thresh: {0.02, 0.005}) {
            const auto cands(sketch_cluster_candidates(sketches, thresh));
            PlacementTree t2(tax, gt);
            const GreedyResult gr(greedy_synthetic_nodes(t2, kc, 64, cands));
            std::fprintf(stderr, "sketch-null-pairs: threshold %g -> %zu candidates, added %zu, residual %llu of %llu\n",
                         thresh, cands.size(), gr.added.size(), (unsigned long long)gr.curve.back().total_delta, (unsigned long long)initial_delta);
            require_monotone(gr);
            REQUIRE(gr.added.size() < 40);
            REQUIRE(gr.curve.back().total_delta > initial_delta / 3);
        }
        kh_destroy(p, tax);
    }
}

// (b) Hash-subsampled disagreement metric vs exact.
TEST_CASE("subsampled_metric") {
    StructuredSet st;
    Spacer sp(17, 17);
    KmerClasses exact;
    exact.fill<score::Lex>(st.paths, sp);
    PlacementTree baseline_tree(st.tax, st.gt);
    const GreedyResult baseline(greedy_synthetic_nodes(baseline_tree, exact, 10));
    const double exact_initial(baseline.curve.front().unweighted);
    require_monotone(baseline);

    for(const double p: {1. / 16, 1. / 256}) {
        double kept(0.);
        KmerClasses sampled(subsample_kmers<score::Lex>(st.paths, sp, p, 777, true, &kept));
        PlacementTree t2(st.tax, st.gt);
        const GreedyResult gr(greedy_synthetic_nodes(t2, sampled, 10));
        const double distortion(std::abs(gr.curve.front().unweighted - exact_initial) / exact_initial);
        std::fprintf(stderr, "subsampled: p=%g kept %.4f initial-estimate %.9g (exact %.9g, distortion %.4f), placements %s\n",
                     p, kept, gr.curve.front().unweighted, exact_initial, distortion,
                     (gr.added_masks == baseline.added_masks ? "identical" : "DIFFER"));
        REQUIRE(gr.added_masks == baseline.added_masks); // structural k-mers forced kept: placements invariant
        REQUIRE(gr.curve.back().total_delta == baseline.curve.back().total_delta);
        REQUIRE(distortion < 0.20);
        REQUIRE(kept < p + 0.20); // structural copies are kept exactly
    }
    SECTION("without forced structure, placements distort") {
        KmerClasses sampled(subsample_kmers<score::Lex>(st.paths, sp, 1. / 16, 777, false));
        PlacementTree t2(st.tax, st.gt);
        const GreedyResult gr(greedy_synthetic_nodes(t2, sampled, 10));
        std::fprintf(stderr, "subsampled-no-force: added %zu (exact %zu), final delta %llu (exact %llu)\n",
                     gr.added_masks.size(), baseline.added_masks.size(), (unsigned long long)gr.curve.back().total_delta, (unsigned long long)baseline.curve.back().total_delta);
        REQUIRE(gr.added_masks != baseline.added_masks);
    }
}
