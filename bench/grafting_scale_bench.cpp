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
        std::unordered_map<tax_t, u32> tmap;
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
        }
    }

    if(csv) std::fclose(csv);
    kh_destroy(p, taxmap);
    return EXIT_SUCCESS;
}
