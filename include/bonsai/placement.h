#ifndef __BNS_PLACEMENT_H
#define __BNS_PLACEMENT_H

/*
Approximation of full k-mer -> target-genome-bitmap encodings by
k-mer -> tree-placement encodings (LCA/placement node implies the set of
genomes under it). Provides the approximation metric, a greedy
add-synthetic-nodes loop with surplus-and-trim, and error-vs-budget
curves. See TODO.md items 1, 2, and the synthetic core of item 3.
*/

#include "util.h"
#include "feature_min.h"

#include <algorithm>
#include <functional>
#include <unordered_map>
#include <unordered_set>
#include <vector>

namespace bns {

INLINE unsigned pc64(u64 x) {return __builtin_popcountll(x);}

// Whether synthetic nodes may pull genomes out of existing taxa.
// unrestricted: adopt whole subtrees, then steal remaining member leaves.
// respect_taxonomy: the grouped set must be a union of whole existing
// subtrees (any depth); candidates that would split taxa are rejected.
enum class GraftPolicy { unrestricted, respect_taxonomy };

// Distinct k-mer target sets: bitmap over genomes -> number of distinct
// k-mers with that bitmap and their total occurrence count (weight).
struct KmerClasses {
    u32 ngenomes;
    u64 total_distinct;                       // sum of count over all bitmaps
    std::unordered_map<u64, u32> count;       // bitmap -> #distinct k-mers
    std::unordered_map<u64, u64> weight;      // bitmap -> total occurrences
    KmerClasses(): ngenomes(0), total_distinct(0) {}
    template<typename ScoreType>
    void fill(const std::vector<std::string> &paths, const Spacer &sp, bool canon=true) {
        ngenomes = paths.size();
        LOG_ASSERT(ngenomes <= 64);
        std::unordered_map<u64, u64> kmer_bitmap, kmer_weight;
        for(u32 i(0); i < ngenomes; ++i) {
            khash_t(64) *counts(kh_init(64));
            fill_count_genome<ScoreType>(paths[i].data(), sp, counts, i, nullptr, canon);
            const u64 bit(u64(1) << i);
            for(khiter_t ki(0); ki != kh_end(counts); ++ki) {
                if(!kh_exist(counts, ki)) continue;
                kmer_bitmap[kh_key(counts, ki)] |= bit;
                kmer_weight[kh_key(counts, ki)] += kh_val(counts, ki);
            }
            kh_destroy(64, counts);
        }
        for(const auto &pair: kmer_bitmap) {
            ++count[pair.second];
            weight[pair.second] += kmer_weight[pair.first];
        }
        for(const auto &pair: count) total_distinct += pair.second;
    }
};

/*
Explicit placement tree. Every node's mask is at all times the union of
its descendant genome leaves, so the implied target set of a k-mer with
true genome set S is the mask of the deepest node containing S (the
placement/LCA node). Synthetic nodes group arbitrary genome sets: maximal
descendant subtrees fully inside the set are adopted, remaining member
leaves are stolen (their old ancestors shrink), so masks stay truthful.
*/
class PlacementTree {
    struct PNode {
        u64 mask;
        u32 parent;              // u32(-1) for the root
        std::vector<u32> children;
        tax_t taxid;             // original taxonomy id for taxonomy nodes, 0 for leaves/synthetic
        bool leaf;
        bool synthetic;
        bool alive;
    };
    std::vector<PNode> nodes_;
    std::vector<u32> leaf_;      // genome index -> leaf node id
    u32 root_;
    u32 ngenomes_;

    // Recompute masks from x upward; drop emptied non-root internal nodes.
    void fixup(u32 x) {
        while(x != u32(-1) && x != root_ && !nodes_[x].leaf && nodes_[x].children.empty()) {
            const u32 p(nodes_[x].parent);
            auto &pc(nodes_[p].children);
            pc.erase(std::find(pc.begin(), pc.end(), x));
            nodes_[x].alive = false;
            x = p;
        }
        while(x != u32(-1)) {
            u64 m(0);
            for(const auto c: nodes_[x].children) m |= nodes_[c].mask;
            if(nodes_[x].leaf || m == nodes_[x].mask) break; // leaves keep their bit; unchanged => ancestors unchanged
            nodes_[x].mask = m;
            x = nodes_[x].parent;
        }
    }

public:
    // tax: child -> parent map (root's parent 0). genome_taxa[i]: taxid of genome i.
    PlacementTree(const khash_t(p) *tax, const std::vector<tax_t> &genome_taxa) {
        ngenomes_ = genome_taxa.size();
        LOG_ASSERT(ngenomes_ <= 64);
        root_ = u32(-1);
        std::unordered_map<tax_t, u32> tmap;
        std::vector<tax_t> pending;
        for(const auto t: genome_taxa) if(tmap.find(t) == tmap.end()) pending.push_back(t);
        while(!pending.empty()) { // expand to ancestor closure, root last
            const tax_t t(pending.back());
            const khint_t ki(kh_get(p, tax, t));
            if(ki == kh_end(tax)) RUNTIME_ERROR(ks::sprintf("Tax ID %u missing from taxonomy.", t).data());
            const tax_t par(kh_val(tax, ki));
            if(par == 0) {
                pending.pop_back();
                // A map may hold several parent-0 entries (e.g. a reformed
                // taxonomy with super-root 1 plus an old root); kraken
                // expects id 1 to be the root, so hang every other
                // parent-0 node under 1.
                const bool has_superroot(kh_get(p, tax, 1) != kh_end(tax));
                if(has_superroot && t != 1 && tmap.find(1) == tmap.end()) {
                    const u32 rid(nodes_.size());
                    nodes_.push_back(PNode{0, u32(-1), {}, 1, false, false, true});
                    tmap[1] = rid;
                    root_ = rid;
                }
                if(tmap.find(t) == tmap.end()) {
                    const u32 pid(has_superroot && t != 1 ? tmap[1]: u32(-1));
                    const u32 id(nodes_.size());
                    nodes_.push_back(PNode{0, pid, {}, t, false, false, true});
                    if(pid != u32(-1)) nodes_[pid].children.push_back(id);
                    tmap[t] = id;
                    if(pid == u32(-1)) {
                        if(root_ != u32(-1)) RUNTIME_ERROR("Taxonomy has multiple roots.");
                        root_ = id;
                    }
                }
            } else if(tmap.find(par) == tmap.end()) {
                pending.push_back(par);
            } else {
                pending.pop_back();
                if(tmap.find(t) == tmap.end()) {
                    const u32 id(nodes_.size());
                    nodes_.push_back(PNode{0, tmap[par], {}, t, false, false, true});
                    nodes_[tmap[par]].children.push_back(id);
                    tmap[t] = id;
                }
            }
        }
        if(root_ == u32(-1)) RUNTIME_ERROR("Taxonomy has no root (no taxid with parent 0).");
        leaf_.resize(ngenomes_);
        for(u32 i(0); i < ngenomes_; ++i) {
            const u32 id(nodes_.size());
            nodes_.push_back(PNode{u64(1) << i, tmap.at(genome_taxa[i]), {}, 0, true, false, true});
            nodes_[tmap.at(genome_taxa[i])].children.push_back(id);
            leaf_[i] = id;
        }
        std::function<void(u32)> post = [&](u32 x) {
            u64 m(0);
            for(const auto c: nodes_[x].children) {
                post(c);
                m |= nodes_[c].mask;
            }
            if(!nodes_[x].leaf) nodes_[x].mask = m;
        };
        post(root_);
    }

    u32 lca(u64 set) const {
        u32 cur(root_);
        for(;;) {
            u32 nxt(u32(-1));
            for(const auto c: nodes_[cur].children)
                if(nodes_[c].alive && (nodes_[c].mask & set) == set) {nxt = c; break;}
            if(nxt == u32(-1)) return cur;
            cur = nxt;
        }
    }
    u64 implied(u64 set) const {return nodes_[lca(set)].mask;}
    u64 err(u64 set) const {return u64(pc64(implied(set) ^ set));}

    // Add a synthetic node grouping `mask`; returns its id, or -1 if the set
    // is already a node or (under GraftPolicy::respect_taxonomy) the set is
    // not a union of whole existing subtrees. Nothing outside the
    // adopted/stolen subtrees moves.
    i64 add_synthetic(u64 mask, GraftPolicy policy = GraftPolicy::unrestricted) {
        if(pc64(mask) < 2) return -1;
        const u32 N(lca(mask));
        if(nodes_[N].mask == mask) return -1;
        std::vector<u32> stack;
        for(const auto c: nodes_[N].children) stack.push_back(c);
        std::vector<u32> adopt;
        while(!stack.empty()) {
            const u32 d(stack.back());
            stack.pop_back();
            if(!nodes_[d].alive) continue;
            if((nodes_[d].mask & ~mask) == 0) adopt.push_back(d);
            else for(const auto c: nodes_[d].children) stack.push_back(c);
        }
        u64 covered(0);
        for(const auto d: adopt)
            if(policy == GraftPolicy::unrestricted || !nodes_[d].leaf) covered |= nodes_[d].mask;
        if(policy == GraftPolicy::respect_taxonomy && (mask & ~covered) != 0)
            return -1; // would need to pull genomes out of existing taxa
        const u32 x(nodes_.size());
        nodes_.push_back(PNode{mask, N, {}, 0, false, true, true});
        nodes_[N].children.push_back(x);
        std::vector<u32> touched;
        for(const auto d: adopt) {
            if(policy == GraftPolicy::respect_taxonomy && nodes_[d].leaf) continue;
            const u32 p(nodes_[d].parent);
            auto &pc(nodes_[p].children);
            pc.erase(std::find(pc.begin(), pc.end(), d));
            nodes_[d].parent = x;
            nodes_[x].children.push_back(d);
            touched.push_back(p);
        }
        for(u64 remaining(mask & ~covered); remaining; remaining &= remaining - 1) {
            const u32 l(leaf_[__builtin_ctzll(remaining)]);
            const u32 p(nodes_[l].parent);
            auto &pc(nodes_[p].children);
            pc.erase(std::find(pc.begin(), pc.end(), l));
            nodes_[l].parent = x;
            nodes_[x].children.push_back(l);
            touched.push_back(p);
        }
        touched.push_back(N);
        u64 check(0);
        for(const auto c: nodes_[x].children) check |= nodes_[c].mask;
        if(check != mask) RUNTIME_ERROR("Synthetic node mask mismatch.");
        for(const auto t: touched) fixup(t);
        nodes_[x].mask = mask; // fixup(N) may have recomputed through x before its children settled
        return x;
    }

    bool remove_synthetic(u32 id) {
        if(id >= nodes_.size() || !nodes_[id].alive || !nodes_[id].synthetic) return false;
        const u32 p(nodes_[id].parent);
        for(const auto c: nodes_[id].children) {
            nodes_[c].parent = p;
            nodes_[p].children.push_back(c);
        }
        auto &pc(nodes_[p].children);
        pc.erase(std::find(pc.begin(), pc.end(), id));
        nodes_[id].alive = false;
        fixup(p);
        return true;
    }

    // Singleton bitmaps always classify exactly (implied == the genome's leaf),
    // so they contribute zero error and are skipped in the sums.
    u64 total_delta(const KmerClasses &kc) const {
        u64 total(0);
        for(const auto &pair: kc.count) {
            if(pc64(pair.first) < 2) continue;
            total += pair.second * err(pair.first);
        }
        return total;
    }
    double unweighted_error(const KmerClasses &kc) const {
        return kc.total_distinct ? double(total_delta(kc)) / kc.total_distinct : 0.;
    }
    double weighted_error(const KmerClasses &kc) const {
        long double num(0), den(0);
        for(const auto &pair: kc.count) {
            den += pair.second * kc.weight.at(pair.first);
            if(pc64(pair.first) < 2) continue;
            num += pair.second * err(pair.first) * kc.weight.at(pair.first);
        }
        return den ? double(num / den) : 0.;
    }

    u64 mask(u32 id) const {return nodes_[id].mask;}
    u32 parent_of(u32 id) const {return nodes_[id].parent;}
    tax_t node_taxid(u32 id) const {return nodes_[id].taxid;}
    bool is_synthetic(u32 id) const {return nodes_[id].synthetic;}
    bool is_leaf_node(u32 id) const {return nodes_[id].leaf;}
    bool is_alive(u32 id) const {return nodes_[id].alive;}
    u32 leaf_node(u32 genome) const {return leaf_[genome];}
    const std::vector<u32> &children_of(u32 id) const {return nodes_[id].children;}
    u32 nnodes() const {return u32(nodes_.size());}
    u32 ngenomes() const {return ngenomes_;}
    u32 root() const {return root_;}
    u64 genome_mask() const {return ngenomes_ == 64 ? ~u64(0): (u64(1) << ngenomes_) - 1;}
    u32 nsynthetic() const {
        u32 ret(0);
        for(const auto &n: nodes_) ret += n.alive && n.synthetic;
        return ret;
    }
};

struct CurvePoint {
    u32 budget;
    double unweighted, weighted;
    u64 total_delta;
};

struct GreedyResult {
    std::vector<CurvePoint> curve;
    std::vector<u32> added;
    std::vector<u64> added_masks;
};

/*
Greedy surplus addition over an explicit candidate mask list: repeatedly
add the candidate with the largest exact error reduction until the
budget is met or no candidate improves. Never increases the unweighted
error (each step is accepted only on positive exact gain).
*/
inline GreedyResult greedy_synthetic_nodes(PlacementTree &tree, const KmerClasses &kc, u32 budget, const std::vector<u64> &candidates, GraftPolicy policy) {
    GreedyResult ret;
    ret.curve.push_back(CurvePoint{0, tree.unweighted_error(kc), tree.weighted_error(kc), tree.total_delta(kc)});
    std::unordered_set<u64> used;
    while(ret.added.size() < budget) {
        const double base(tree.unweighted_error(kc));
        double best_gain(1e-9);
        i64 best(-1);
        for(u32 ci(0); ci < candidates.size(); ++ci) {
            const u64 m(candidates[ci]);
            if(used.find(m) != used.end()) continue;
            PlacementTree tmp(tree);
            if(tmp.add_synthetic(m, policy) < 0) {used.insert(m); continue;}
            const double gain(base - tmp.unweighted_error(kc));
            // Store the candidate index: a mask with bit 63 set is negative as
            // i64 and would be mistaken for "no candidate".
            if(gain > best_gain) {best_gain = gain; best = i64(ci);}
        }
        if(best < 0) break;
        const u64 m = candidates[u32(best)];
        const i64 id(tree.add_synthetic(m, policy));
        LOG_ASSERT(id >= 0);
        used.insert(m);
        ret.added.push_back(u32(id));
        ret.added_masks.push_back(m);
        ret.curve.push_back(CurvePoint{u32(ret.added.size()), tree.unweighted_error(kc), tree.weighted_error(kc), tree.total_delta(kc)});
    }
    return ret;
}

inline GreedyResult greedy_synthetic_nodes(PlacementTree &tree, const KmerClasses &kc, u32 budget, const std::vector<u64> &candidates) {
    return greedy_synthetic_nodes(tree, kc, budget, candidates, GraftPolicy::unrestricted);
}

/*
Greedy surplus addition: repeatedly add the candidate (a distinct observed
non-singleton k-mer bitmap) with the largest exact error reduction until
the budget is met or no candidate improves. Never increases the
unweighted error (each step is accepted only on positive exact gain).
*/
inline GreedyResult greedy_synthetic_nodes(PlacementTree &tree, const KmerClasses &kc, u32 budget) {
    std::vector<u64> candidates;
    candidates.reserve(kc.count.size());
    for(const auto &pair: kc.count) if(pc64(pair.first) >= 2) candidates.push_back(pair.first);
    return greedy_synthetic_nodes(tree, kc, budget, candidates);
}

struct TrimResult {
    u32 kept, removed;
    double final_unweighted;
    std::vector<u32> nodes;
};

/*
Trim surplus synthetic nodes whose removal degrades the unweighted error
by at most rel_tolerance times the error at trim start. Nodes are
considered in reverse addition order.
*/
inline TrimResult trim_synthetic_nodes(PlacementTree &tree, const KmerClasses &kc, const std::vector<u32> &surplus_ids, double rel_tolerance) {
    TrimResult ret{0, 0, 0., {}};
    double cur(tree.unweighted_error(kc));
    const double tol(rel_tolerance * std::max(cur, 1e-9));
    ret.final_unweighted = cur;
    for(auto it(surplus_ids.rbegin()); it != surplus_ids.rend(); ++it) {
        const u32 id(*it);
        if(id == u32(-1)) continue;
        PlacementTree tmp(tree);
        if(!tmp.remove_synthetic(id)) {ret.nodes.push_back(id); continue;}
        const double newe(tmp.unweighted_error(kc));
        if(newe - cur <= tol) {
            tree = tmp;
            cur = newe;
            ++ret.removed;
        } else ret.nodes.push_back(id);
    }
    ret.kept = u32(ret.nodes.size());
    ret.final_unweighted = cur;
    return ret;
}

} // namespace bns

#endif // __BNS_PLACEMENT_H
