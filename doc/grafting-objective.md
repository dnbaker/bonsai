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
