# Bonsai grafting as optimization: formal objective

## Problem

Given the taxonomic tree T over leaves L (target entities) and the multiset of
k-mers with their true occurrence bitmaps S_x ⊆ L (the exhaustive answer):

- **True structure**: each k-mer x ↦ S_x (arbitrary subset of L).
- **Tree effective structure**: T's encoding can only produce **subtree-closed
  bitmaps** — the laminar family of clades {leaves(v) : v ∈ T}. Each k-mer is
  reported at some node v(x); its effective bitmap is leaves(v(x))).

**Design goal**: choose ≤ B grafted (synthetic) nodes to add to T,
minimizing the total Hamming distance

    H_total = Σ_x [ H(S_x, leaves(v_T'(x))) ]

i.e. false-negatives |S_x \ leaves(v)| + false-positives |leaves(v) \ S_x|,
where T' is the augmented tree and v_T'(x) is x's reporting node in it.

Two assignment models, and the guarantee story differs between them:

1. **Fixed-slot model** (matches the hash-to-leaf encoding): each k-mer's slot
   in T is fixed by hashing; a graft with leaf-set A directly absorbs the
   k-mers whose slots fall in A, upgrading their effective bitmap from
   leaves(parent position) to A. Benefit of graft set G = Σ_x benefit_x(best
   graft covering x's slot). This is maximum-coverage-shaped: if each
   graft's per-k-mer benefit is nonnegative and grafts interact only through
   best-of selection, the total benefit is **monotone submodular** → greedy
   gives (1 − 1/e) ≈ 63% of optimal for uniform graft costs, with budget B.
2. **Best-fit model** (k-mers free to report at their Hamming-nearest node):
   adding a node can *steal* k-mers from other grafts → interactions are not
   obviously submodular; greedy is still the right heuristic but the
   (1−1/e) proof needs an additional assumption (e.g. grafts at disjoint
   scales / laminar-compatible graft sets) or it fails — produce a small
   counterexample if so, and state which restriction restores the bound.

## Feasible grafts (keep the family laminar)

- **Split graft**: carve A ⊊ clade(v) out of an existing clade (splits v into
  A and clade(v)\A).
- **Union graft**: A = union of some children of an existing node (a new
  internal node grouping existing clades).
Both keep the family laminar, so the encoding machinery (kraken2/krakenuniq
export, tree walk) is unchanged.

## Algorithm (greedy)

Repeat B times: evaluate every feasible graft A at its optimal placement,
score by ΔH_total (exact incremental evaluation on the k-mer index), take the
max. With sketches: estimate S_x and |leaves(v) ∩ S_x| from the per-genome
BBMinHash<b=16> sketches (this is where the sketch lane feeds in — Hamming
distance between a clade bitmap and a k-mer's true bitmap becomes register-set
operations), exact-verify the top candidates on the k-mer index.

## Guarantee program (what the theory deliverable must contain)

1. Prove (1 − 1/e) under the fixed-slot model with stated assumptions; state
   the assumptions in terms of the actual encoding (which model bonsai uses —
   the implementer must nail this from the code, it determines everything).
2. If best-fit is the true model: either find the restoring restriction or
   exhibit a counterexample to submodularity; in the latter case give the
   weaker guarantee that does hold (e.g. greedy = 1/2·OPT for non-submodular
   monotone via the standard forward-backward argument if it applies, or
   nothing beyond monotone improvement — stated honestly).
3. Empirical gap check: brute-force optimum (small instances, B ≤ 3, |L| ≤ 8)
   vs greedy — is the real gap near (1−1/e), near 1.0, or does the
   counterexample geometry occur in practice?

## Evaluation

- Metric: H_total per k-mer (mean), FN/FP decomposition, vs B (graft budget
  curve — the money figure: Hamming distance vs #new nodes).
- Baselines: no grafts; random feasible grafts; greedy-without-sketch
  (exact counting).
- Validation end-to-end: kraken2 classification accuracy on a labeled
  metagenome before/after grafting at matched B.

## Implementation and results (2026-10-10)

### Operative model

The objective evaluation in placement.h (`implied`/`err`) reports each
k-mer class at the deepest node whose mask contains it, and grafts adopt
whole subtrees and steal member leaves (`add_synthetic`,
`GraftPolicy::unrestricted`). This is the **best-fit model**: a graft can
move classes that report at other grafts. The fixed-slot model corresponds
to hash-to-leaf encoding with reporting slots frozen in the original tree.

### Fixed-slot guarantee (proof)

Under fixed-slot semantics, each k-mer class x keeps a fixed leaf slot
s(x); a graft A absorbs the classes with s(x) in A and upgrades their
effective bitmap from the parent-position mask P(s(x)) to A. Define
b(x, A) = c_x * [H(S_x, P(s(x))) - H(S_x, A)] for s(x) in A (0 otherwise),
with b >= 0 by placing each graft optimally and rejecting non-improving
placements. For a graft set G let

    f(G) = sum_x max(0, max_{A in G, s(x) in A} b(x, A))

i.e. each class counts its best covering graft once. Then:

- f is monotone (a max over a growing set with a nonnegative basis);
- f is submodular: each h_A(G) = b(x, A) * 1[A in G] is modular, and for
  any two sets S, T the argmax achievers c, d of F(S cup T), F(S cap T)
  satisfy F(S cup T) + F(S cap T) = h_c(S cup T) + h_d(S cap T) <= F(S) +
  F(T) (c contributes to S or T; d to both); summation over x preserves
  submodularity.
- Therefore (Nemhauser-Wolsey-Fisher) greedy with uniform costs and budget
  B achieves >= (1 - 1/e) of the optimal reduction.

The coverage ILP below is the exact reference for this model.

### Best-fit: counterexample and restriction

The best-fit gain is **not submodular**. Minimal counterexample (locked in
test/grafting.cpp): 4 leaves, species {0,1} and {2,3}; k-mer classes
X = {2,3} (species-internal, placed exactly) and Y = {1,2} (cross-species).
With B = {0,1,2}:

- g(B | empty) = -1: B steals leaf 2 from species 101, breaking X
  (err 0 -> 2) while partially fixing Y (err 2 -> 1).
- g(B | A={0,2}) = +1: A has already stolen leaf 2 and absorbed the
  breakage, so B's marginal only sees Y's improvement.

Marginal gains increased (1 > -1), so diminishing returns fail. The
destruction channel is stealing that breaks an exactly-placed class.
Randomized searches back both sides: 5.45M pair-marginal checks over
random sparse-class instances (n <= 6) found zero violations, so the
regime is dense overlapping structure with exactly-placed classes sharing
leaves with cross-taxon structure — which is also where the greedy-exact
gap appears in practice (below).

**Restoring restriction.** If grafts must be unions of whole existing
subtrees and mutually laminar (no stealing: every graft contains every
class it touches in full, e.g. GraftPolicy::respect_taxonomy over
whole-subtree sets), the destruction channel closes, each class moves at
most from its current reporting node to a tighter graft, and the
fixed-slot proof applies verbatim (slot = current reporting node). Under
that restriction the (1 - 1/e) bound holds.

**What holds in general.** Monotone improvement only (greedy accepts
strictly positive exact gains), hence H_total is non-increasing per
graft. No constant-factor guarantee is claimed for unrestricted best-fit.

### Exact reference and gap measurements

`graft_exact` (include/bonsai/grafting.h) is branch-and-bound over graft
subsets in index order under the true tree dynamics. No lookahead bound is
used: marginals are not submodular, so single-gain bounds are not
admissible. ILP formulation (fixed-slot coverage form; the best-fit
objective is not linear in z, which is why the exact reference is the B&B
over subsets):

    variables: z_A in {0,1} for each feasible graft A;
               y_{x,A} in {0,1} (class x assigned to graft A)
    maximize  sum_{x,A} b(x,A) y_{x,A}
    s.t.      sum_A y_{x,A} <= 1        (each class counted once, best-of)
              y_{x,A} <= z_A
              sum_A z_A <= B
              z, y binary

Measured greedy-vs-exact reduction ratios (B&B over all feasible subsets,
8 leaves unless noted): B=1: 1.000; B=2: 0.833 on the module fixture and
0.900 on a second 8-leaf instance (the gap is exactly the destruction
geometry); B=3: 1.000; gap-free on 6 leaves and on all seeded randomized
instances. Greedy stays above the (1 - 1/e) threshold on every tested
instance despite the bound not applying.

### Evaluation summary (bench: doc/grafting_bench.md)

- H vs B curve: on the 32-genome module fixtures the observed-class
  candidate lane drives mean Hamming distance per k-mer from 0.168 to 0.099
  (B=1), 0.040 (B=2), 0 (B=3) on the cross-genus fixture, and
  0.069 -> 0.040 -> 0.020 -> 0 on the genus-local fixture; error is
  exclusively false positives at baseline (classes sit inside implied
  clades); false negatives only arise through stealing-induced shrinkage.
- Candidate lanes: the feasible union/split family and the BBMinHash
  tree-layered sketch lane both match the observed-class lane end-to-end
  on genus-local modules, and both provably miss cross-genus modules
  (sketch dilution at the top tree layer; non-taxon shapes absent from the
  family) - see doc/grafting_bench.md.
- Time: at POC scale the exhaustive k-mer index fill dominates total wall
  clock (0.1 s vs <= 0.05 s for sketch build + tree-layered candidates +
  all greedy scoring); candidate evaluation is never the bottleneck at this
  scale. The sketch lane is what removes the index requirement at scale.
