/*
Grafting bench (doc/grafting_bench.md): H_total vs budget curve with FN/FP
decomposition, greedy vs exact optimum gap on small instances, baselines
(no grafts / random grafts / feasible-family / observed-class / sketch
candidates), and a wall-clock breakdown of the objective-evaluation
pipeline at POC scale (32 genomes). Single thread; g++-15 -O3 -march=native.

Fixture A: cross-genus module (16 genomes across two genera) — the
tree-layered sketch path dilutes cross-genus agreement at the root layer
and cannot propose the full module (documented limitation); observed-class
candidates reach zero error.
Fixture B: genus-local modules — both the sketch path and the feasible
union/split family propose them.
*/
#include "encoder.h"
#include "spacer.h"
#include "grafting.h"
#include "sketch_select.h"
#include "util.h"
#include "aesctr/wy.h"

#include <algorithm>
#include <chrono>
#include <cmath>
#include <cstdio>
#include <numeric>
#include <string>
#include <vector>

using namespace bns;

namespace {

double now_s() {
    return std::chrono::duration<double>(std::chrono::steady_clock::now().time_since_epoch()).count();
}

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
        char tmpl[] = "/tmp/bns_gbench_XXXXXX";
        return std::string(mkdtemp(tmpl));
    }()) {}
    ~TmpDir() {
        for(const auto &f: files) std::remove(f.data());
        std::remove(dir.data());
    }
    std::string add(const std::string &name, const std::string &content) {
        files.push_back(dir + "/" + name);
        std::FILE *fp(std::fopen(files.back().data(), "w"));
        std::fputs(content.data(), fp);
        std::fclose(fp);
        return files.back();
    }
};

khash_t(p) *build_synth_tax() {
    khash_t(p) *ret(kh_init(p));
    int khr;
    auto put = [&](tax_t k, tax_t v) {const khint_t ki(kh_put(p, ret, k, &khr)); kh_val(ret, ki) = v;};
    put(1, 0); put(200, 1); put(201, 1);
    for(unsigned i(0); i < 16; ++i) put(100 + i, i < 8 ? 200 : 201);
    return ret;
}

struct Fixture {
    TmpDir tmp;
    std::vector<std::string> paths;
    std::vector<tax_t> gt;
    khash_t(p) *tax;
    explicit Fixture(u64 seed, bool genus_local): tax(build_synth_tax()) {
        wy::WyRand<u64> rng(seed);
        std::string species_blocks[16];
        for(unsigned i(0); i < 16; ++i) species_blocks[i] = rand_dna(rng, 60);
        const std::string mA_block(rand_dna(rng, 120)), mB_block(rand_dna(rng, 120)), mC_block(rand_dna(rng, 120));
        // A: 16 genomes across both genera (genus_local = false: species
        // 0-3 and 8-11) or 8 genomes inside genus 200 (species 0-3).
        auto carries_A = [&](unsigned sp) {return genus_local ? (sp < 4) : (sp < 4 || (sp >= 8 && sp < 12));};
        auto carries_B = [&](unsigned sp) {return genus_local ? (sp >= 4 && sp < 8) : (sp == 4 || sp == 5 || sp == 12 || sp == 13);};
        auto carries_C = [&](unsigned sp) {return genus_local ? (sp >= 8 && sp < 10) : (sp == 6 || sp == 14);};
        for(unsigned g(0); g < 32; ++g) {
            const unsigned sp(g >> 1);
            std::string content(">g\n" + rand_dna(rng, 1300) + "\n");
            content += ">sp\n" + species_blocks[sp] + "\n";
            if(carries_A(sp)) content += ">mA\n" + mA_block + "\n";
            else if(carries_B(sp)) content += ">mB\n" + mB_block + "\n";
            else if(carries_C(sp)) content += ">mC\n" + mC_block + "\n";
            paths.push_back(tmp.add("g" + std::to_string(g) + ".fa", content));
            gt.push_back(100 + sp);
        }
    }
    ~Fixture() {kh_destroy(p, tax);}
};

void print_curve(const char *tag, const GraftingResult &gr, const KmerClasses &kc) {
    std::printf("%-18s", tag);
    for(size_t i(0); i < gr.deltas.size(); ++i)
        std::printf("  B=%zu H=%.5f fn=%llu fp=%llu", i,
                    kc.total_distinct ? double(gr.deltas[i]) / kc.total_distinct : 0.,
                    (unsigned long long)gr.fn[i], (unsigned long long)gr.fp[i]);
    std::printf("\n");
}

void run_fixture(const char *name, u64 seed, bool genus_local) {
    Spacer sp(17, 17);
    std::printf("\n== %s ==\n", name);
    Fixture st(seed, genus_local);
    KmerClasses kc;
    double t0(now_s());
    kc.fill<score::Lex>(st.paths, sp);
    const double fill_t(now_s() - t0);
    PlacementTree base(st.tax, st.gt);
    std::printf("k-mer index fill: %.3f s (%zu classes, %llu distinct, initial delta %llu)\n",
                fill_t, kc.count.size(), (unsigned long long)kc.total_distinct,
                (unsigned long long)base.total_delta(kc));
    const std::vector<u64> fam(enumerate_grafts(base, false));
    std::vector<u64> class_cands;
    for(const auto &pair: kc.count) if(pc64(pair.first) >= 2) class_cands.push_back(pair.first);
    for(const char *tag: {"feasible-set", "observed-classes"}) {
        const std::vector<u64> &cands(std::string(tag) == "feasible-set" ? fam : class_cands);
        PlacementTree tree(st.tax, st.gt);
        t0 = now_s();
        const GraftingResult gr(graft_greedy(tree, kc, cands, 8));
        const double g_t(now_s() - t0);
        print_curve(tag, gr, kc);
        std::printf("%-18s wall %.3f s (%zu grafts)\n", tag, g_t, gr.added_masks.size());
    }
    { // random grafts: unconditionally apply 3 random class candidates
        double acc(0.);
        const unsigned trials(20);
        for(u64 seed2(1); seed2 <= trials; ++seed2) {
            wy::WyRand<u64> rng(seed2);
            PlacementTree tree(st.tax, st.gt);
            for(unsigned i(0); i < 3; ++i) tree.add_synthetic(class_cands[rng() % class_cands.size()]);
            acc += double(tree.total_delta(kc));
        }
        std::printf("random-3 mean final delta %.1f\n", acc / trials);
    }
    { // sketch lane: BBMH tree-layered candidates -> exact-verified greedy
        t0 = now_s();
        const auto sks(genome_sketchers<score::Lex>(st.paths, sp, 8192));
        const double sk_t(now_s() - t0);
        t0 = now_s();
        const auto comps(sketch_cluster_candidates_tree(sks, base, 0.008));
        std::vector<u64> cands;
        for(const auto &c: comps) {
            u64 m(0);
            for(const auto g: c) m |= u64(1) << g;
            cands.push_back(m);
        }
        const double cg_t(now_s() - t0);
        PlacementTree tree(st.tax, st.gt);
        t0 = now_s();
        const GraftingResult gr(graft_greedy(tree, kc, cands, 8));
        const double g_t(now_s() - t0);
        print_curve("sketch-cands", gr, kc);
        std::printf("sketch lane: sketch build %.3f s, candidates %.4f s (%zu), greedy %.3f s\n",
                    sk_t, cg_t, cands.size(), g_t);
        // ranking quality: within-candidate shared-k-mer estimate (sketch
        // inclusion-exclusion) vs exact first-step gain, top-5 overlap
        auto within_est = [&](u64 A) {
            double sum(0.);
            GenomeSketch uni(sks[__builtin_ctzll(A)].getp(), GENOME_SKETCH_B);
            for(u32 g(0); g < 32; ++g) if(A & (u64(1) << g)) {
                sum += sks[g].cardinality_estimate();
                uni += sks[g];
            }
            return sum - uni.cardinality_estimate();
        };
        std::vector<std::pair<double, u64>> est;
        for(const u64 A: class_cands) est.push_back({within_est(A), A});
        std::sort(est.begin(), est.end(), [](const auto &a, const auto &b){return a.first > b.first;});
        const GraftIndex gi(build_graft_index(base, kc));
        std::vector<std::pair<long double, u64>> exact;
        for(const u64 A: class_cands) exact.push_back({gain_of_add(base, gi, A), A});
        std::sort(exact.begin(), exact.end(), [](const auto &a, const auto &b){return a.first > b.first;});
        size_t overlap(0);
        for(size_t i(0); i < 5 && i < est.size(); ++i)
            for(size_t j(0); j < 5 && j < exact.size(); ++j) if(est[i].second == exact[j].second) ++overlap;
        std::printf("ranking: top-5 sketch within-est vs top-5 exact gain overlap %zu/25\n", overlap);
    }
}

} // namespace

int main() {
    run_fixture("A: cross-genus module", 1337, false);
    run_fixture("B: genus-local modules", 1337, true);
    std::printf("\n== greedy vs exact (8-leaf module instance) ==\n");
    Spacer sp(17, 17);
    TmpDir tmp2;
    wy::WyRand<u64> rng(777);
    khash_t(p) *tax8(kh_init(p));
    {
        int khr;
        auto put = [&](tax_t k, tax_t v) {const khint_t ki(kh_put(p, tax8, k, &khr)); kh_val(tax8, ki) = v;};
        put(1, 0); put(200, 1); put(201, 1);
        for(unsigned i(0); i < 4; ++i) put(100 + i, i < 2 ? 200 : 201);
    }
    std::vector<std::string> paths8;
    std::vector<tax_t> gt8;
    const std::string b0(rand_dna(rng, 120)), b1(rand_dna(rng, 120)), b2(rand_dna(rng, 120));
    for(unsigned g(0); g < 8; ++g) {
        std::string content(">g\n" + rand_dna(rng, 1400) + "\n");
        if(g < 4 && (g & 1) == 0) content += ">m\n" + b0 + "\n";   // {0,2}
        if(g >= 4 && (g & 1) == 1) content += ">m\n" + b1 + "\n";  // {5,7}
        if((g == 1) || (g == 4)) content += ">m\n" + b2 + "\n";    // {1,4}
        paths8.push_back(tmp2.add("g" + std::to_string(g) + ".fa", content));
        gt8.push_back(100 + (g >> 1));
    }
    KmerClasses kc8;
    kc8.fill<score::Lex>(paths8, sp);
    const std::vector<u64> c8(enumerate_grafts(PlacementTree(tax8, gt8), true));
    for(const u32 budget: {1u, 2u, 3u}) {
        PlacementTree gt_tree(tax8, gt8);
        const GraftingResult gr(graft_greedy(gt_tree, kc8, c8, budget));
        const ExactGraftResult ex(graft_exact(PlacementTree(tax8, gt8), kc8, c8, budget));
        const long double rg((long double)(gr.initial_delta - gr.deltas.back()));
        const long double re((long double)(gr.initial_delta - ex.final_delta));
        std::printf("B=%u: greedy %llu exact %llu reduction-ratio %.3Lf\n", budget,
                    (unsigned long long)gr.deltas.back(), (unsigned long long)ex.final_delta,
                    re > 0 ? rg / re : 1.0L);
    }
    kh_destroy(p, tax8);
    return 0;
}
