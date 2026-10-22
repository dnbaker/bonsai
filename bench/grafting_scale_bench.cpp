/*
Real-scale grafting bench over a GTDB subset (doc/bonsai-gtdb-benchmark.md).

Measures the sketch-lane candidate machinery against the exhaustive
per-genome index fill as the genome count grows:

  Stage A (fill projection): cumulative wall time of the LcaMap per-genome
    set fill (fill_set_genome) over a growing genome prefix - the stage the
    sketch lane replaces.
  Stage B (sketch lane): parallel per-genome BBMinHash<b=16> build over all
    genomes, then tree-layered sibling-clade candidates at a threshold grid.
  Stage C (dilution): sibling-pair clade Jaccards recorded per tree depth -
    where cross-clade agreement survives to the comparison layer, and where
    it dilutes below threshold.
  Stage D (sub-clade H vs B): PlacementTree + KmerClasses over the largest
    genus sub-clade (<=64 genomes, u64 masks), sketch + observed-class
    candidates, exact-verified greedy curve.

The candidate path runs on a mask-free adjacency tree, so it is not subject
to the 64-genome PlacementTree cap. Wall clock only.

Usage: grafting_scale_bench <taxmap> <names> <genomes.list> [flags]
  -k: k-mer size [17]
  -s: sketch registers [8192]
  -p: threads for the sketch fill [8]
  -N: sub-clade max genomes [32]
  -c: append CSV rows (section,metric,value) to this file
*/
#include "hashmap.h"
#include "encoder.h"
#include "spacer.h"
#include "feature_min.h"
#include "grafting.h"
#include "sketch_select.h"
#include "util.h"

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

std::vector<std::string> split(const std::string &s, char delim) {
    std::vector<std::string> ret;
    size_t pos(0);
    while(pos <= s.size()) {
        const size_t next(s.find(delim, pos));
        ret.push_back(s.substr(pos, next == std::string::npos ? std::string::npos: next - pos));
        if(next == std::string::npos) break;
        pos = next + 1;
    }
    return ret;
}

void csv_row(std::FILE *fp, const char *section, const char *metric, double value) {
    if(fp) std::fprintf(fp, "%s,%s,%.6g\n", section, metric, value);
}

/*
Mask-free taxonomic tree over an arbitrary genome count, exposing the
subset of the PlacementTree interface the layered candidate generator
uses, so sketch_cluster_candidates_tree runs at any scale.
*/
struct LightTaxTree {
    struct Node {
        u32 parent;
        std::vector<u32> children;
        tax_t taxid;
        bool leaf;
    };
    std::vector<Node> nodes;
    std::vector<u32> leaf_of;
    u32 root_;

    LightTaxTree(const khash_t(p) *tax, const std::vector<tax_t> &genome_taxa) {
        bns::umap<tax_t, u32> tmap;
        std::vector<tax_t> pending(genome_taxa.begin(), genome_taxa.end());
        root_ = u32(-1);
        while(!pending.empty()) {
            const tax_t t(pending.back());
            if(tmap.find(t) != tmap.end()) {pending.pop_back(); continue;}
            const khint_t ki(kh_get(p, tax, t));
            if(ki == kh_end(tax)) RUNTIME_ERROR(ks::sprintf("Tax ID %u missing from taxonomy.", t).data());
            const tax_t par(kh_val(tax, ki));
            if(par == 0 || t == 1) {
                pending.pop_back();
                if(tmap.find(t) == tmap.end()) {
                    if(t == 1) {
                        if(root_ == u32(-1)) {
                            root_ = nodes.size();
                            nodes.push_back(Node{u32(-1), {}, 1, false});
                        }
                        tmap[1] = root_;
                    } else {
                        // non-1 root entries hang under node 1
                        if(root_ == u32(-1)) RUNTIME_ERROR("Super-root 1 missing while a root-level taxon appeared.");
                        const u32 id(nodes.size());
                        nodes.push_back(Node{root_, {}, t, false});
                        nodes[root_].children.push_back(id);
                        tmap[t] = id;
                    }
                }
            } else if(tmap.find(par) == tmap.end()) {
                pending.push_back(par);
            } else {
                pending.pop_back();
                const u32 id(nodes.size());
                nodes.push_back(Node{tmap[par], {}, t, false});
                nodes[tmap[par]].children.push_back(id);
                tmap[t] = id;
            }
        }
        if(root_ == u32(-1)) RUNTIME_ERROR("Taxonomy has no root.");
        leaf_of.resize(genome_taxa.size());
        for(u32 i(0); i < genome_taxa.size(); ++i) {
            const u32 id(nodes.size());
            nodes.push_back(Node{tmap.at(genome_taxa[i]), {}, 0, true});
            nodes[tmap.at(genome_taxa[i])].children.push_back(id);
            leaf_of[i] = id;
        }
    }
    u32 ngenomes() const {return u32(leaf_of.size());}
    u32 nnodes() const {return u32(nodes.size());}
    u32 root() const {return root_;}
    u32 leaf_node(u32 g) const {return leaf_of[g];}
    bool is_leaf_node(u32 x) const {return nodes[x].leaf;}
    u32 parent_of(u32 x) const {return nodes[x].parent;}
    const std::vector<u32> &children_of(u32 x) const {return nodes[x].children;}
    u32 depth(u32 x) const {
        u32 d(0);
        while(x != u32(-1) && nodes[x].parent != u32(-1)) ++d, x = nodes[x].parent;
        return d;
    }
    u64 descendant_genomes(u32 x) const {
        u64 ret(0);
        std::function<void(u32)> dfs = [&](u32 y) {
            if(nodes[y].leaf) {++ret; return;}
            for(const auto c: nodes[y].children) dfs(c);
        };
        dfs(x);
        return ret;
    }
};

struct SketchFillHelper {
    const Spacer &sp;
    const std::vector<std::string> &paths;
    unsigned p, b;
    std::vector<GenomeSketch> *sks;
};
void sketch_fill_fn(void *data_, long index, int tid) {
    (void)tid;
    SketchFillHelper &h(*static_cast<SketchFillHelper *>(data_));
    Encoder<score::Lex> enc(0, 0, h.sp, nullptr, true);
    GenomeSketch bb(h.p, h.b);
    enc.for_each([&](u64 v) {bb.addh(v);}, h.paths[index].data());
    // Each index is processed exactly once, so += into its own slot is
    // race-free with no lock.
    (*(h.sks))[u32(index)] += bb;
}

} // namespace

int main(int argc, char *argv[]) {
    int c, k(17), nthreads(8), subclade_max(32);
    size_t sketch_size(8192);
    const char *csv_path(nullptr);
    while((c = getopt(argc, argv, "k:s:p:N:c:h?")) >= 0) {
        switch(c) {
            case 'h': case '?': goto usage;
            case 'k': k = std::atoi(optarg); break;
            case 's': sketch_size = std::strtoull(optarg, nullptr, 10); break;
            case 'p': nthreads = std::atoi(optarg); break;
            case 'N': subclade_max = std::atoi(optarg); break;
            case 'c': csv_path = optarg; break;
        }
    }
    if(argc - optind < 3) {
        usage:
        std::fprintf(stderr, "Usage: %s <taxmap> <names> <genomes.list> [-k 17] [-s 8192] [-p 8] [-N 32] [-c out.csv]\n", *argv);
        return EXIT_FAILURE;
    }
    std::FILE *csv(csv_path ? std::fopen(csv_path, "a"): nullptr);
    khash_t(p) *taxmap(build_parent_map(argv[optind]));
    std::vector<std::string> paths;
    {
        khash_t(name) *name_hash(build_name_hash(argv[optind + 1]));
        std::FILE *fp(std::fopen(argv[optind + 2], "r"));
        if(!fp) RUNTIME_ERROR("Could not open genome list");
        char *line(nullptr);
        size_t cap(0);
        while(getline(&line, &cap, fp) >= 0) {
            while(*line && std::isspace(line[std::strlen(line) - 1])) line[std::strlen(line) - 1] = 0;
            if(*line == 0 || *line == '#') continue;
            std::string seqid;
            const tax_t tid(first_header_taxid(line, name_hash, &seqid));
            if(tid == tax_t(-1)) {
                std::fprintf(stderr, "[W] %s not in name map; skipping.\n", line);
                continue;
            }
            paths.emplace_back(line);
        }
        std::free(line);
        std::fclose(fp);
        destroy_name_hash(name_hash);
    }
    const u32 n(u32(paths.size()));
    std::fprintf(stderr, "[I] %u genomes\n", n);
    if(n == 0) return EXIT_FAILURE;

    khash_t(name) *name_hash(build_name_hash(argv[optind + 1]));
    std::vector<tax_t> genome_taxa(n);
    for(u32 i(0); i < n; ++i) {
        std::string seqid;
        genome_taxa[i] = first_header_taxid(paths[i].data(), name_hash, &seqid);
    }
    destroy_name_hash(name_hash);

    LightTaxTree tree(taxmap, genome_taxa);
    std::fprintf(stderr, "[I] tree: %u nodes\n", tree.nnodes());

    Spacer sp(k, k);

    // ---- Stage A: per-genome index fill projection ----
    {
        std::fprintf(stderr, "[I] stage A: LcaMap per-genome fill projection\n");
        const unsigned grid[] {32, 64, 128, 256, 512, 1024};
        unsigned gi(0);
        double cum(0.);
        khash_t(all) *set(kh_init(all));
        const u32 an(std::min<u32>(n, 1024));
        for(u32 i(0); i < an; ++i) {
            kh_clear(all, set);
            const double g0(now_s());
            LcaMap::fill<score::Lex>(paths[i].data(), sp, set, i, nullptr, true, nullptr);
            cum += now_s() - g0;
            while(gi < sizeof(grid) / sizeof(grid[0]) && i + 1 == std::min(grid[gi], an)) {
                const unsigned g(std::min(grid[gi], an));
                std::fprintf(stderr, "[I]   fill %4u genomes: cumulative %8.1f s (%6.3f s/genome)\n", g, cum, cum / g);
                char metric[64];
                std::snprintf(metric, sizeof(metric), "fill_seconds_n%u", g);
                csv_row(csv, "fill", metric, cum);
                ++gi;
            }
        }
        kh_destroy(all, set);
    }

    // ---- Stage B: sketch lane at scale ----
    std::vector<GenomeSketch> sks;
    {
        const unsigned p(sketch::ilog2(sketch_size));
        std::fprintf(stderr, "[I] stage B: BBMH b=%u p=%u build over %u genomes (%d threads)\n", GENOME_SKETCH_B, p, n, nthreads);
        double t0(now_s());
        sks.assign(n, GenomeSketch(p, GENOME_SKETCH_B));
        SketchFillHelper helper{sp, paths, p, GENOME_SKETCH_B, &sks};
        {
            ForPool pool(nthreads);
            pool.forpool(&sketch_fill_fn, &helper, n);
        }
        const double bt(now_s() - t0);
        std::fprintf(stderr, "[I]   sketch build %.1f s (%.3f s/genome)\n", bt, bt / n);
        csv_row(csv, "sketch", "build_seconds", bt);
        for(const double thr: {0.002, 0.008, 0.02}) {
            std::vector<std::pair<u32, double>> jaccards;
            t0 = now_s();
            const auto comps(sketch_cluster_candidates_tree(sks, tree, thr, thr == 0.008 ? &jaccards: nullptr));
            const double ct(now_s() - t0);
            u32 multi(0);
            for(const auto &comp: comps) multi += comp.size() >= 2;
            std::fprintf(stderr, "[I]   threshold %.3f: %zu components in %.1f s\n", thr, comps.size(), ct);
            char metric[64];
            std::snprintf(metric, sizeof(metric), "candidates_seconds_t%g", thr);
            csv_row(csv, "sketch", metric, ct);
            std::snprintf(metric, sizeof(metric), "candidates_count_t%g", thr);
            csv_row(csv, "sketch", metric, double(comps.size()));
            if(thr == 0.008 && !jaccards.empty()) {
                // ---- Stage C: dilution by depth ----
                std::fprintf(stderr, "[I] stage C: %zu sibling-clade pairs; per-depth Jaccard\n", jaccards.size());
                std::vector<std::vector<double>> by_depth(16);
                for(const auto &pr: jaccards) {
                    const u32 d(std::min(tree.depth(pr.first), 15u));
                    by_depth[d].push_back(pr.second);
                }
                for(u32 d(0); d < by_depth.size(); ++d) {
                    const auto &v(by_depth[d]);
                    if(v.empty()) continue;
                    const double mx(*std::max_element(v.begin(), v.end()));
                    const double mean(std::accumulate(v.begin(), v.end(), 0.) / v.size());
                    u32 above(0);
                    for(const double j: v) above += j >= 0.002;
                    std::fprintf(stderr, "[I]   depth %2u: pairs %7zu  max %.4f  mean %.5f  >=0.002: %u\n",
                                 d, v.size(), mx, mean, above);
                    char metric[64];
                    std::snprintf(metric, sizeof(metric), "depth%u_max", d);
                    csv_row(csv, "dilution", metric, mx);
                    std::snprintf(metric, sizeof(metric), "depth%u_mean", d);
                    csv_row(csv, "dilution", metric, mean);
                    std::snprintf(metric, sizeof(metric), "depth%u_pairs", d);
                    csv_row(csv, "dilution", metric, double(v.size()));
                }
            }
        }
    }

    // ---- Stage D: sub-clade H vs B ----
    {
        std::fprintf(stderr, "[I] stage D: sub-clade grafting curve (max %d genomes)\n", subclade_max);
        // genus = depth-6 node (d=1,p=2,c=3,o=4,f=5,g=6) with the most
        // descendant genomes.
        u32 best(u32(-1)), best_count(0);
        for(u32 x(0); x < tree.nnodes(); ++x) {
            if(tree.depth(x) != 6 || tree.is_leaf_node(x)) continue;
            const u64 dg(tree.descendant_genomes(x));
            if(dg >= 2 && dg > best_count) {best = x; best_count = u32(std::min<u64>(dg, UINT32_MAX));}
        }
        if(best == u32(-1)) {
            std::fprintf(stderr, "[W] no genus-level sub-clade found; skipping stage D\n");
        } else {
            std::vector<u32> members;
            std::function<void(u32)> dfs = [&](u32 y) {
                if(tree.is_leaf_node(y)) {
                    const u32 g(u32(std::find(tree.leaf_of.begin(), tree.leaf_of.end(), y) - tree.leaf_of.begin()));
                    members.push_back(g);
                    return;
                }
                for(const auto ch: tree.children_of(y)) dfs(ch);
            };
            dfs(best);
            std::fprintf(stderr, "[I]   genus node %u (taxid %u): %zu genomes\n", best, tree.nodes[best].taxid, members.size());
            if(members.size() > u32(subclade_max)) {
                members.resize(subclade_max); // arbitrary but deterministic: tree order
            }
            std::vector<std::string> sub_paths;
            std::vector<tax_t> sub_taxa;
            for(const auto g: members) {
                sub_paths.push_back(paths[g]);
                sub_taxa.push_back(genome_taxa[g]);
            }
            const u32 m(u32(sub_paths.size()));
            PlacementTree ptree(taxmap, sub_taxa);
            KmerClasses kc;
            double t0(now_s());
            kc.fill<score::Lex>(sub_paths, sp, true, false);
            const double fill_t(now_s() - t0);
            const u64 delta0(ptree.total_delta(kc));
            std::fprintf(stderr, "[I]   sub-clade fill %.1f s: %zu classes, %llu distinct, initial delta %llu (H %.5f)\n",
                         fill_t, kc.count.size(), (unsigned long long)kc.total_distinct,
                         (unsigned long long)delta0, kc.total_distinct ? double(delta0) / kc.total_distinct: 0.);
            csv_row(csv, "subclade", "fill_seconds", fill_t);
            csv_row(csv, "subclade", "ngenomes", double(m));
            csv_row(csv, "subclade", "classes", double(kc.count.size()));
            csv_row(csv, "subclade", "initial_delta", double(delta0));
            // candidate lanes: observed classes + tree-layered sketch
            std::vector<u64> cands;
            for(const auto &pair: kc.count) if(pc64(pair.first) >= 2) cands.push_back(pair.first);
            const size_t nclasses(cands.size());
            t0 = now_s();
            const auto sks_sub(genome_sketchers<score::Lex>(sub_paths, sp));
            const auto comps(sketch_cluster_candidates_tree(sks_sub, ptree, 0.008));
            for(const auto &comp: comps) {
                u64 mask(0);
                for(const auto g: comp) mask |= u64(1) << g;
                cands.push_back(mask);
            }
            const double cand_t(now_s() - t0);
            std::fprintf(stderr, "[I]   candidates: %zu classes + %zu sketch in %.2f s\n", nclasses, cands.size() - nclasses, cand_t);
            PlacementTree run_tree(taxmap, sub_taxa);
            t0 = now_s();
            const GraftingResult gr(graft_greedy(run_tree, kc, cands, 8));
            const double greedy_t(now_s() - t0);
            for(size_t i(0); i < gr.deltas.size(); ++i) {
                std::fprintf(stderr, "[I]   B=%zu H=%.5f fn=%llu fp=%llu\n", i,
                             kc.total_distinct ? double(gr.deltas[i]) / kc.total_distinct: 0.,
                             (unsigned long long)gr.fn[i], (unsigned long long)gr.fp[i]);
                char metric[64];
                std::snprintf(metric, sizeof(metric), "H_B%zu", i);
                csv_row(csv, "subclade", metric, kc.total_distinct ? double(gr.deltas[i]) / kc.total_distinct: 0.);
            }
            std::fprintf(stderr, "[I]   greedy %zu grafts in %.2f s\n", gr.added_masks.size(), greedy_t);
            csv_row(csv, "subclade", "greedy_seconds", greedy_t);
            csv_row(csv, "subclade", "grafts", double(gr.added_masks.size()));
            csv_row(csv, "subclade", "candidates", double(cands.size()));
            csv_row(csv, "subclade", "scanned", double(gr.scanned));
            csv_row(csv, "subclade", "bound_skipped", double(gr.bound_skipped));
            csv_row(csv, "subclade", "exact_evals", double(gr.exact_evals));
            // Stage E: the feasible union/split family lane (candidate-scan
            // cost assessment). Union grafts are 2^fanout-2 per node, capped
            // at fanout 20 - the count itself is the finding at GTDB shapes.
            t0 = now_s();
            const std::vector<u64> fam(enumerate_grafts(ptree));
            const double enum_t(now_s() - t0);
            std::fprintf(stderr, "[I]   stage E: feasible family %zu grafts enumerated in %.2f s\n", fam.size(), enum_t);
            csv_row(csv, "family", "candidates", double(fam.size()));
            csv_row(csv, "family", "enum_seconds", enum_t);
            if(fam.size() <= 400000) {
                PlacementTree fam_tree(taxmap, sub_taxa);
                t0 = now_s();
                const GraftingResult fr(graft_greedy(fam_tree, kc, fam, 8));
                const double fam_t(now_s() - t0);
                for(size_t i(0); i < fr.deltas.size(); ++i) {
                    std::fprintf(stderr, "[I]   family B=%zu H=%.5f\n", i,
                                 kc.total_distinct ? double(fr.deltas[i]) / kc.total_distinct: 0.);
                }
                std::fprintf(stderr, "[I]   family greedy %zu grafts in %.2f s: scanned %llu bound-skipped %llu exact %llu\n",
                             fr.added_masks.size(), fam_t, (unsigned long long)fr.scanned,
                             (unsigned long long)fr.bound_skipped, (unsigned long long)fr.exact_evals);
                {
                    std::string rs;
                    for(double r: fr.round_seconds) rs += std::to_string(r) + " ";
                    std::fprintf(stderr, "[I]   family round-s: %s\n", rs.data());
                    u32 maxf(0);
                    for(u32 x(0); x < ptree.nnodes(); ++x)
                        maxf = std::max<u32>(maxf, u32(ptree.children_of(x).size()));
                    std::fprintf(stderr, "[I]   placement tree: %u nodes, max fanout %u\n", ptree.nnodes(), maxf);
                }
                csv_row(csv, "family", "greedy_seconds", fam_t);
                csv_row(csv, "family", "scanned", double(fr.scanned));
                csv_row(csv, "family", "bound_skipped", double(fr.bound_skipped));
                csv_row(csv, "family", "exact_evals", double(fr.exact_evals));
                csv_row(csv, "family", "grafts", double(fr.added_masks.size()));
            } else {
                std::fprintf(stderr, "[I]   family lane skipped over the 400k candidate cap (union grafts are 2^fanout)\n");
            }
        }
    }

    // ---- Stage F: high-fanout lanes (GTDB-like genus) ----
    // A genus with >20 species children is common in GTDB. Round 2 measured
    // the cliff: above fanout 20 the union lane was dropped entirely and the
    // greedy stalled at ~8% of initial delta. Round 3 lanes, all polynomial:
    //   family   - enumerate_grafts: pairs + triplets of children above the cap
    //   seeded   - enumerate_seeded_grafts: child-union hulls of the observed classes
    //   observed - the class masks themselves (fanout-independent)
    for(const u32 fanout: {16u, 24u, 40u}) {
        const u32 ngen(fanout);   // one genome per species: fanout genomes under one genus
        khash_t(p) *big_tax(kh_init(p));
        int khr;
        kh_put(p, big_tax, 1, &khr); kh_val(big_tax, kh_get(p, big_tax, 1)) = 0;
        kh_put(p, big_tax, 2, &khr); kh_val(big_tax, kh_get(p, big_tax, 2)) = 1;
        for(u32 s(0); s < fanout; ++s) {
            const tax_t sp_tax(100 + s);
            kh_put(p, big_tax, sp_tax, &khr);
            kh_val(big_tax, kh_get(p, big_tax, sp_tax)) = 2;
        }
        std::vector<tax_t> gt(ngen);
        for(u32 s(0); s < fanout; ++s) gt[s] = tax_t(100 + s);
        PlacementTree big_tree(big_tax, gt);
        // classes: each species' k-mers unique (clean) + one cross-species
        // module spanning species 0..fanout/2-1 (the graftable signal)
        KmerClasses bkc;
        bkc.ngenomes = ngen;
        std::vector<u64> observed;
        for(u32 s(0); s < fanout; ++s) {
            bkc.count[u64(1) << s] = 1000;
            observed.push_back(u64(1) << s);
        }
        bkc.count[(u64(1) << (fanout / 2)) - 1] = 5000;
        observed.push_back((u64(1) << (fanout / 2)) - 1);
        bkc.total_distinct = u64(fanout) * 1000 + 5000;
        u32 maxf(0);
        for(u32 x(0); x < big_tree.nnodes(); ++x)
            maxf = std::max<u32>(maxf, u32(big_tree.children_of(x).size()));
        std::fprintf(stderr, "[I] stage F fanout=%2u: %u-node tree (max fanout %u)\n", fanout, big_tree.nnodes(), maxf);
        const auto lane = [&](const char *name, const std::vector<u64> &cands, u32 budget) {
            if(cands.empty()) return;
            PlacementTree run_tree(big_tax, gt);
            const double tg0(now_s());
            const GraftingResult fr(graft_greedy(run_tree, bkc, cands, budget));
            const double tgf(now_s() - tg0);
            std::fprintf(stderr, "[I]   %-8s %8zu cands B<=%u: delta %llu -> %llu in %.2f s (scanned %llu exact %llu)\n",
                         name, cands.size(), budget, (unsigned long long)fr.initial_delta,
                         (unsigned long long)fr.deltas.back(), tgf,
                         (unsigned long long)fr.scanned, (unsigned long long)fr.exact_evals);
            char metric[64];
            std::snprintf(metric, sizeof(metric), "%s_delta_f%u_b%u", name, fanout, budget);
            csv_row(csv, "fanoutlane", metric, double(fr.deltas.back()));
            std::snprintf(metric, sizeof(metric), "%s_seconds_f%u_b%u", name, fanout, budget);
            csv_row(csv, "fanoutlane", metric, tgf);
        };
        double t0f(now_s());
        const std::vector<u64> famf(enumerate_grafts(big_tree));
        const double enumf(now_s() - t0f);
        t0f = now_s();
        const std::vector<u64> seedf(enumerate_seeded_grafts(big_tree, observed));
        const double enums(now_s() - t0f);
        std::fprintf(stderr, "[I]   family %zu grafts (%.3f s), seeded %zu grafts (%.3f s)\n",
                     famf.size(), enumf, seedf.size(), enums);
        csv_row(csv, "fanout", "family_candidates", double(famf.size()));
        csv_row(csv, "fanout", "seeded_candidates", double(seedf.size()));
        lane("family", famf, 4);
        lane("seeded", seedf, 4);
        lane("observed", observed, 4);
        lane("famB16", famf, 16);   // pairs/triplets compound hierarchically with budget
        lane("seedB16", seedf, 16);
        kh_destroy(p, big_tax);
    }

    // ---- Stage G: sketch-only vs index-seeded naming, head-to-head ----
    // PI directive 2026-10-21 (round 4): naming from bbminhash intersection
    // cardinalities alone - no k-mer sets in candidate generation or gain
    // path. The exact lane below builds the oracle KmerClasses ONLY to
    // evaluate final H and the estimator bias; the sketch lane never sees it.
    bool excl_mode(false);
    int r7_mode(0);                    // 1 covers, 2 replication, 3 multiseed
    const ExactCovers *r7_covers(nullptr);
    auto head2head = [&](khash_t(p) *htax, const std::vector<std::string> &paths, const std::vector<tax_t> &taxa,
                         const char *label, double tau) {
        std::fprintf(stderr, "[I] stage G [%s]: %zu genomes, tau=%g\n", label, paths.size(), tau);
        PlacementTree ptree(htax, taxa);
        // oracle (evaluation + index lane only)
        KmerClasses kc;
        double t0(now_s());
        kc.fill<score::Lex>(paths, sp, true, false);
        const double index_fill(now_s() - t0);
        // index-seeded lane: observed classes -> hulls -> exact greedy
        std::vector<u64> class_masks;
        for(const auto &pair: kc.count) if(pc64(pair.first) >= 2) class_masks.push_back(pair.first);
        t0 = now_s();
        const std::vector<u64> seeded(enumerate_seeded_grafts(ptree, class_masks));
        PlacementTree idx_tree(htax, taxa);
        const GraftingResult idx(graft_greedy(idx_tree, kc, seeded, 8));
        const double idx_wall(now_s() - t0);
        // sketch-only lane
        t0 = now_s();
        const std::vector<GenomeSketch> sks(genome_sketchers<score::Lex>(paths, sp, sketch_size));
        SketchNodePolicy skpol;
        skpol.exclusive_edges = excl_mode;
        skpol.pair_candidates = excl_mode;
        SketchIndex skidx(build_sketch_index(ptree, sks, skpol));
        if(r7_mode == 1) skidx.covers = r7_covers;
        std::vector<SketchGraft> sk_cands(enumerate_sketch_seeded_grafts(skidx, tau));
        if(r7_mode == 2) {   // replication: fresh builds at p=13/16/19 must agree
            const unsigned ptop(std::min(19u, sketch::ilog2(sketch_size)));
            std::vector<std::vector<GenomeSketch>> budget(3);   // ptop, -3, -6
            for(unsigned b(0); b < 3; ++b) {
                const unsigned p(ptop > 6 * b ? ptop - 3 * b: ptop);
                for(const auto &path: paths) {
                    budget[b].emplace_back(p, GENOME_SKETCH_B);
                    Encoder<score::Lex> enc(0, 0, sp, nullptr, true);
                    enc.for_each([&](u64 km) {budget[b].back().addh(km);}, path.data());
                }
            }
            std::vector<SketchGraft> kept;
            for(const auto &c: sk_cands) {
                if(c.parts.size() != 2) {kept.push_back(c); continue;}
                std::vector<double> ests;
                for(const auto &set: budget) {
                    const auto A(clade_sketch(set, c.parts[0]));
                    const auto B(clade_sketch(set, c.parts[1]));
                    std::vector<const uint64_t *> cx;
                    for(u64 mm(((paths.size() == 64 ? ~0ull: (u64(1) << paths.size()) - 1)) & ~(c.parts[0] | c.parts[1])); mm; mm &= mm - 1)
                        cx.push_back(set[__builtin_ctzll(mm)].core().data());
                    ests.push_back(differential_from_sketches(A.core().data(), B.core().data(), cx, set[0].core().size(), set[0].getp()));
                }
                std::fprintf(stderr, "[I]     repl mask %06llx: %.0f %.0f %.0f\n",
                             (unsigned long long)c.mask, ests[0], ests[1], ests[2]);
                const double mx(std::max({ests[0], ests[1], ests[2]})), mn(std::min({ests[0], ests[1], ests[2]}));
                if(mx > 0. && mn > 0. && mx <= 3. * mn) kept.push_back(c);
            }
            sk_cands.swap(kept);
        }
        if(r7_mode == 3) {   // multi-seed CV (R=3, p=min(19, sketch p))
            const unsigned p(std::min(19u, sketch::ilog2(sketch_size)));
            using SeedSk = sketch::BBitMinHasher<uint64_t, sketch::SeededHash<sketch::WangHash>>;
            std::vector<std::vector<SeedSk>> seeded;
            for(unsigned seed(1); seed <= 3; ++seed) {
                seeded.emplace_back();
                for(const auto &path: paths) {
                    seeded.back().emplace_back(p, GENOME_SKETCH_B, uint64_t(seed));
                    Encoder<score::Lex> enc(0, 0, sp, nullptr, true);
                    enc.for_each([&](u64 km) {seeded.back().back().addh(km);}, path.data());
                }
            }
            std::vector<SketchGraft> kept;
            for(const auto &c: sk_cands) {
                if(c.parts.size() != 2) {kept.push_back(c); continue;}
                std::vector<double> ests;
                for(unsigned seed(1); seed <= 3; ++seed) {
                    const auto &set(seeded[seed - 1]);
                    const auto A(clade_sketch(set, c.parts[0]));
                    const auto B(clade_sketch(set, c.parts[1]));
                    std::vector<const uint64_t *> cx;
                    for(u64 mm(((paths.size() == 64 ? ~0ull: (u64(1) << paths.size()) - 1)) & ~(c.parts[0] | c.parts[1])); mm; mm &= mm - 1)
                        cx.push_back(set[__builtin_ctzll(mm)].core().data());
                    ests.push_back(differential_from_sketches(A.core().data(), B.core().data(), cx, size_t(1) << p, p));
                }
                const double mean((ests[0] + ests[1] + ests[2]) / 3.);
                const double sd(std::sqrt((std::pow(ests[0] - mean, 2) + std::pow(ests[1] - mean, 2) + std::pow(ests[2] - mean, 2)) / 3.));
                const double cv(mean > 0. ? sd / mean: 9.);
                std::fprintf(stderr, "[I]     mseed mask %06llx: %.0f %.0f %.0f (cv %.2f)\n",
                             (unsigned long long)c.mask, ests[0], ests[1], ests[2], cv);
                if(mean > 0. && cv <= 0.5) kept.push_back(c);
            }
            sk_cands.swap(kept);
        }
        for(const auto &c: sk_cands)
            std::fprintf(stderr, "[I]     sk cand mask %06llx parts %zu est %.0f node %u\n",
                         (unsigned long long)c.mask, c.parts.size(), c.est_gain, c.node);
        PlacementTree sk_tree(htax, taxa);
        const GraftingResult sk(graft_greedy_sketch(sk_tree, sks, sk_cands, 8, r7_covers));
        const double sk_wall(now_s() - t0);
        // oracle H curves for both graft sequences
        auto replay = [&](const std::vector<u64> &grafts) {
            PlacementTree t(htax, taxa);
            std::vector<double> hs{kc.total_distinct ? double(t.total_delta(kc)) / kc.total_distinct: 0.};
            for(const u64 A: grafts) {
                if(t.add_synthetic(A) >= 0) {}
                hs.push_back(kc.total_distinct ? double(t.total_delta(kc)) / kc.total_distinct: 0.);
            }
            return hs;
        };
        std::fprintf(stderr, "[I]   index grafts:");
        for(const u64 A: idx.added_masks) std::fprintf(stderr, " %02llx", (unsigned long long)A);
        std::fprintf(stderr, "\n[I]   sketch grafts:");
        for(const u64 A: sk.added_masks) std::fprintf(stderr, " %02llx", (unsigned long long)A);
        std::fprintf(stderr, "\n");
        const auto hi(replay(idx.added_masks)), hsk(replay(sk.added_masks));
        std::fprintf(stderr, "[I]   index : %zu hulls, %zu grafts, fill %.2f s + greedy %.2f s | H:", seeded.size(), idx.added_masks.size(), index_fill, idx_wall);
        for(const double h: hi) std::fprintf(stderr, " %.5f", h);
        std::fprintf(stderr, "\n[I]   sketch: %zu cands, %zu grafts, %.2f s | H:", sk_cands.size(), sk.added_masks.size(), sk_wall);
        for(const double h: hsk) std::fprintf(stderr, " %.5f", h);
        std::fprintf(stderr, "\n");
        // estimator bias/variance: sketch est vs exact marginal per accepted graft
        if(hsk.size() > 1) {
            std::fprintf(stderr, "[I]   sketch est vs exact marginal:");
            for(size_t i(0); i < sk.added_masks.size() && i + 1 < hsk.size(); ++i) {
                const double exact_marg((hsk[i] - hsk[i + 1]) * kc.total_distinct);
                std::fprintf(stderr, " %.0f/%.0f", sk.est_gains[i], exact_marg);
                char metric[64];
                std::snprintf(metric, sizeof(metric), "%s_est_f%u", label, unsigned(i));
                csv_row(csv, "skbias", metric, sk.est_gains[i]);
                std::snprintf(metric, sizeof(metric), "%s_exact_f%u", label, unsigned(i));
                csv_row(csv, "skbias", metric, exact_marg);
            }
            std::fprintf(stderr, "\n");
        }
        char metric[64];
        // Stage I (real-data inventory, BNS_INVENTORY=1): for each accepted
        // sketch-lane graft, extract the exclusive k-mer content (present in
        // every carrier, absent from every other genome) and locate it on the
        // first carrier's records. Analysis only - the lane never sees this.
        if(std::getenv("BNS_INVENTORY") && paths.size() <= 64) {
            auto genome_kmers = [&](u32 g) {
                std::vector<uint64_t> v;
                Encoder<score::Lex> enc(0, 0, sp, nullptr, true);
                enc.for_each([&](uint64_t km) {v.push_back(km);}, paths[g].data());
                std::sort(v.begin(), v.end());
                v.erase(std::unique(v.begin(), v.end()), v.end());
                return v;
            };
            u32 gi(0);
            for(const u64 A: sk.added_masks) {
                // exclusive set: intersect carriers, subtract non-carriers
                std::vector<uint64_t> mod(genome_kmers(u32(__builtin_ctzll(A))));
                for(u64 m(A & (A - 1)); m; m &= m - 1) {
                    const auto o(genome_kmers(u32(__builtin_ctzll(m))));
                    std::vector<uint64_t> inter;
                    std::set_intersection(mod.begin(), mod.end(), o.begin(), o.end(), std::back_inserter(inter));
                    mod.swap(inter);
                }
                for(u64 m(((paths.size() == 64 ? ~0ull: (u64(1) << paths.size()) - 1)) & ~A); m; m &= m - 1) {
                    const auto o(genome_kmers(u32(__builtin_ctzll(m))));
                    std::vector<uint64_t> diff;
                    std::set_difference(mod.begin(), mod.end(), o.begin(), o.end(), std::back_inserter(diff));
                    mod.swap(diff);
                }
                // locate on the first carrier: per-record positions
                std::fprintf(stderr, "[I]   graft %u: carriers", gi);
                for(u64 m(A); m; m &= m - 1) {
                    const u32 g(__builtin_ctzll(m));
                    std::string base(paths[g]);
                    base = base.substr(base.find_last_of('/') + 1);
                    std::fprintf(stderr, " %s", base.data());
                }
                std::fprintf(stderr, " (%zu exclusive k-mers) -> spans:", mod.size());
                std::FILE *fa(std::fopen(paths[u32(__builtin_ctzll(A))].data(), "r"));
                if(fa) {
                    std::string rec, seq, name;
                    auto flush = [&]() {
                        if(seq.empty()) return;
                        Encoder<score::Lex> enc(0, 0, sp, nullptr, true);
                        size_t pos(0);
                        u64 pstart(u64(-1));
                        auto mark = [&](u64 km) {
                            const bool in(std::binary_search(mod.begin(), mod.end(), km));
                            if(in && pstart == u64(-1)) pstart = pos;
                            if(!in && pstart != u64(-1)) {
                                if(pos + 16 > pstart + 1000)
                                    std::fprintf(stderr, " [%s %llu-%llu (%llu bp)]", name.data(),
                                                 (unsigned long long)pstart, (unsigned long long)pos + 16,
                                                 (unsigned long long)(pos + 16 - pstart));
                                pstart = u64(-1);
                            }
                            ++pos;
                        };
                        enc.for_each(mark, seq.data(), seq.size());
                        if(pstart != u64(-1) && pos + 16 > pstart + 1000)
                            std::fprintf(stderr, " [%s %llu-%llu (%llu bp)]", name.data(),
                                         (unsigned long long)pstart, (unsigned long long)pos + 16,
                                         (unsigned long long)(pos + 16 - pstart));
                    };
                    int c;
                    while((c = std::fgetc(fa)) != EOF) {
                        if(c == '>') {
                            flush();
                            seq.clear();
                            rec.clear();
                            while((c = std::fgetc(fa)) != EOF && c != '\n') rec.push_back(char(c));
                            name = rec.substr(0, rec.find(' '));
                        } else if(c != '\n' && c != '\r') seq.push_back(char(c));
                    }
                    flush();
                    std::fclose(fa);
                }
                std::fprintf(stderr, "\n");
                ++gi;
            }
        }
        std::snprintf(metric, sizeof(metric), "%s_H_sketch", label);
        csv_row(csv, "head2head", metric, hsk.back());
        std::snprintf(metric, sizeof(metric), "%s_H_index", label);
        csv_row(csv, "head2head", metric, hi.back());
        std::snprintf(metric, sizeof(metric), "%s_sketch_seconds", label);
        csv_row(csv, "head2head", metric, sk_wall);
        std::snprintf(metric, sizeof(metric), "%s_index_seconds", label);
        csv_row(csv, "head2head", metric, index_fill + idx_wall);
    };

    auto head2head_excl = [&](khash_t(p) *htax, const std::vector<std::string> &paths, const std::vector<tax_t> &taxa,
                              const char *label) {
        excl_mode = true;
        head2head(htax, paths, taxa, label, -1.);
        excl_mode = false;
    };

    // G1: the round-4 probe (single genus, planted cross-species module)
    if(std::getenv("BNS_PROBE_DIR")) {
        const std::string pdir(std::getenv("BNS_PROBE_DIR"));
        khash_t(p) *ptax(build_parent_map((pdir + "/taxmap.txt").data()));
        std::vector<std::string> ppaths;
        std::vector<tax_t> ptaxa;
        {
            khash_t(name) *pnh(build_name_hash((pdir + "/names.txt").data()));
            std::FILE *fp(std::fopen((pdir + "/genomes.list").data(), "r"));
            if(!fp) RUNTIME_ERROR("could not open probe genomes.list");
            char *line(nullptr);
            size_t cap(0);
            while(getline(&line, &cap, fp) >= 0) {
                while(*line && std::isspace(line[std::strlen(line) - 1])) line[std::strlen(line) - 1] = 0;
                if(*line == 0 || *line == '#') continue;
                std::string seqid;
                const tax_t tid(first_header_taxid(line, pnh, &seqid));
                if(tid == tax_t(-1)) continue;
                ppaths.emplace_back(line);
                ptaxa.push_back(tid);
            }
            std::free(line);
            std::fclose(fp);
            destroy_name_hash(pnh);
        }
        for(const double tau: {16., 64., 256., 512.}) head2head(ptax, ppaths, ptaxa, "probe", tau);
        kh_destroy(p, ptax);
    } else {
        std::fprintf(stderr, "[I] stage G: BNS_PROBE_DIR not set; probe head-to-head skipped\n");
    }
    // G2: the main dataset's largest genus sub-clade (realistic shape)
    if(n > 1) {
        u32 best(u32(-1)), best_count(0);
        for(u32 x(0); x < tree.nnodes(); ++x) {
            if(tree.depth(x) != 6 || tree.is_leaf_node(x)) continue;
            const u64 dg(tree.descendant_genomes(x));
            if(dg >= 2 && dg > best_count) {best = x; best_count = u32(std::min<u64>(dg, UINT32_MAX));}
        }
        if(best != u32(-1)) {
            std::vector<u32> members;
            std::function<void(u32)> dfs = [&](u32 y) {
                if(tree.is_leaf_node(y)) {
                    members.push_back(u32(std::find(tree.leaf_of.begin(), tree.leaf_of.end(), y) - tree.leaf_of.begin()));
                    return;
                }
                for(const auto ch: tree.children_of(y)) dfs(ch);
            };
            dfs(best);
            if(members.size() > 64) members.resize(64);
            std::vector<std::string> sub_paths;
            std::vector<tax_t> sub_taxa;
            for(const auto g: members) {
                sub_paths.push_back(paths[g]);
                sub_taxa.push_back(genome_taxa[g]);
            }
            for(const double tau: {16., 64., 256., 512.}) head2head(taxmap, sub_paths, sub_taxa, "gtdb_genus", tau);
            // real pangenomes: raw-count edges are core-genome dominated;
            // rerun the adaptive arm with exclusive-edge mode (pair-private)
            head2head_excl(taxmap, sub_paths, sub_taxa, "gtdb_genus_excl");
            // round 7 arms (one per run, BNS_MODE): covers | replication | multiseed
            if(const char *m = std::getenv("BNS_MODE")) {
                const std::string mode(m);
                auto build_covers = [&](const PlacementTree &pt, const std::vector<std::string> &pths) {
                    ExactCovers cov;
                    for(u32 g(0); g < pths.size(); ++g) {
                        auto v(std::make_shared<std::vector<uint64_t>>());
                        Encoder<score::Lex> enc(0, 0, sp, nullptr, true);
                        enc.for_each([&](u64 km) {v->push_back(km);}, pths[g].data());
                        std::sort(v->begin(), v->end());
                        v->shrink_to_fit();
                        cov.full[pt.leaf_node(g)] = v;
                        cov.core[pt.leaf_node(g)] = v;   // single genome: core == full
                    }
                    std::function<void(u32)> up = [&](u32 x) {
                        if(pt.is_leaf_node(x)) return;
                        for(const u32 c: pt.children_of(x)) up(c);
                        const auto &ch(pt.children_of(x));
                        if(ch.size() == 1) {
                            const auto it(cov.full.find(ch[0]));
                            const auto jt(cov.core.find(ch[0]));
                            cov.full[x] = it == cov.full.end() ? nullptr: it->second;
                            cov.core[x] = jt == cov.core.end() ? nullptr: jt->second;
                            return;
                        }
                        std::vector<std::shared_ptr<std::vector<uint64_t>>> sets, isets;
                        for(const u32 c: ch) sets.push_back(cov.full[c]);
                        for(const u32 c: ch) isets.push_back(cov.core[c]);
                        std::sort(sets.begin(), sets.end(), [](const auto &a, const auto &b) {return a->size() < b->size();});
                        std::sort(isets.begin(), isets.end(), [](const auto &a, const auto &b) {return a->size() < b->size();});
                        auto merged(std::make_shared<std::vector<uint64_t>>(*sets[0]));
                        auto inter(std::make_shared<std::vector<uint64_t>>(*isets[0]));
                        std::vector<uint64_t> tmp;
                        for(size_t i(1); i < sets.size(); ++i) {
                            tmp.clear();
                            std::set_union(merged->begin(), merged->end(), sets[i]->begin(), sets[i]->end(), std::back_inserter(tmp));
                            merged->swap(tmp);
                            tmp.clear();
                            std::set_intersection(inter->begin(), inter->end(), isets[i]->begin(), isets[i]->end(), std::back_inserter(tmp));
                            inter->swap(tmp);
                        }
                        merged->shrink_to_fit();
                        inter->shrink_to_fit();
                        cov.full[x] = merged;
                        cov.core[x] = inter;
                    };
                    up(pt.root());
                    return cov;
                };
                ExactCovers cov;
                const std::string label(std::string("r7_") + mode);
                if(mode == "covers") {
                    cov = build_covers(PlacementTree(taxmap, sub_taxa), sub_paths);
                    std::fprintf(stderr, "[I] exact covers: %zu bytes over %zu nodes\n", cov.bytes(), cov.full.size());
                }
                r7_mode = mode == "covers" ? 1: mode == "replication" ? 2: mode == "multiseed" ? 3: 0;
                r7_covers = mode == "covers" ? &cov: nullptr;
                excl_mode = true;
                head2head(taxmap, sub_paths, sub_taxa, label.data(), -1.);
                excl_mode = false;
            }
        }
    }

    // ---- Stage H: adaptive tau vs absolute tau across background regimes ----
    // Round 4's measured failure: a 120 bp module is named at 8 kb backgrounds
    // and lost in the 30 kb chance-sharing floor under an absolute tau. The
    // adaptive per-node threshold (doc/grafting-two-tier.md) must flip the
    // 30 kb case to named. Probe dirs come from BNS_H_PROBES (colon-separated).
    if(const char *dirs = std::getenv("BNS_H_PROBES")) {
        for(const std::string pdir: split(dirs, ':')) {
            if(pdir.empty()) continue;
            khash_t(p) *ptax(build_parent_map((pdir + "/taxmap.txt").data()));
            std::vector<std::string> ppaths;
            std::vector<tax_t> ptaxa;
            {
                khash_t(name) *pnh(build_name_hash((pdir + "/names.txt").data()));
                std::FILE *fp(std::fopen((pdir + "/genomes.list").data(), "r"));
                if(!fp) RUNTIME_ERROR("could not open probe genomes.list");
                char *line(nullptr);
                size_t cap(0);
                while(getline(&line, &cap, fp) >= 0) {
                    while(*line && std::isspace(line[std::strlen(line) - 1])) line[std::strlen(line) - 1] = 0;
                    if(*line == 0 || *line == '#') continue;
                    std::string seqid;
                    const tax_t tid(first_header_taxid(line, pnh, &seqid));
                    if(tid == tax_t(-1)) continue;
                    ppaths.emplace_back(line);
                    ptaxa.push_back(tid);
                }
                std::free(line);
                std::fclose(fp);
                destroy_name_hash(pnh);
            }
            PlacementTree ptree(ptax, ptaxa);
            KmerClasses kc;
            kc.fill<score::Lex>(ppaths, sp, true, false);
            // p must scale with the input so n/m stays near or below 1:
            // 8 kb genomes resolve at p=13, 30 kb genomes need p=16 (the
            // round-4 failure was register resolution, not tau - measured in
            // doc/grafting-two-tier.md). Sweep both to show it.
            for(const unsigned p_reg: {13u, 16u}) {
                const std::vector<GenomeSketch> sks(genome_sketchers<score::Lex>(ppaths, sp, size_t(1) << p_reg));
                const SketchIndex idx(build_sketch_index(ptree, sks));
                std::fprintf(stderr, "[I] stage H [%s] p=%u: %zu genomes\n",
                             pdir.substr(pdir.find_last_of('/') + 1).data(), p_reg, ppaths.size());
                auto run_arm = [&](const char *name, const std::vector<SketchGraft> &cands) {
                    PlacementTree t(ptax, ptaxa);
                    const GraftingResult gr(graft_greedy_sketch(t, sks, cands, 8));
                    PlacementTree o(ptax, ptaxa);
                    for(const u64 A: gr.added_masks) if(o.add_synthetic(A) >= 0) {}
                    const double h0(kc.total_distinct ? double(o.total_delta(kc)) / kc.total_distinct: 0.);
                    std::fprintf(stderr, "[I]   %-9s: %2zu cands, %zu grafts, final H %.5f\n", name, cands.size(), gr.added_masks.size(), h0);
                    char metric[96];
                    std::snprintf(metric, sizeof(metric), "%s_p%u_%s_H", pdir.substr(pdir.find_last_of('/') + 1).data(), p_reg, name);
                    csv_row(csv, "stageH", metric, h0);
                };
                {
                    const auto c16(enumerate_sketch_seeded_grafts(idx, 16.));
                    for(const auto &c: c16) std::fprintf(stderr, "[I]     abs16 cand est %.0f tau %.0f mask %06llx\n", c.est_gain, c.tau, (unsigned long long)c.mask);
                    run_arm("abs16", c16);
                    const auto cad(enumerate_sketch_seeded_grafts(idx));
                    for(const auto &c: cad) std::fprintf(stderr, "[I]     adapt cand est %.0f tau %.0f mask %06llx\n", c.est_gain, c.tau, (unsigned long long)c.mask);
                    run_arm("adaptive", cad);
                }
            }
            kh_destroy(p, ptax);
        }
    }

    if(csv) std::fclose(csv);
    kh_destroy(p, taxmap);
    return EXIT_SUCCESS;
}
