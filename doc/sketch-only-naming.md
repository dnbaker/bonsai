# Sketch-only naming for grafting (PI directive 2026-10-21, round 4)

The round-3 seeded-hull lane names high-arity union grafts from observed
k-mer classes - which requires the exact per-genome index. This note covers
the sketch-only replacement: candidate generation and gain estimation from
per-genome bbminhash **intersection cardinalities alone** - no exact k-mer
sets, no inverted-list index, anywhere in the path. The scaling law: huge
databases afford sketches, not indexes.

Implementation (`include/bonsai/grafting.h`):

- `enumerate_sketch_seeded_grafts(tree, leaf_sketches, tau)`: at each
  internal node v, cluster v's children by their clade-sketch pairwise
  intersection profile (edge when the estimated pair intersection clears
  `tau`); connected components of size >= 2 become graft masks, scored by
  the sketch-only gain. O(sum_v k_v^2) pair intersections of O(2^p) work.
- `sketch_exclusive_estimate(sks, parts, parent)`: gain proxy = core minus
  contamination: the minimum pairwise intersection over the component's
  member clades (the shared core), minus the maximum pairwise intersection
  of a member with the context clade (P \ A). Registers only.
- `graft_greedy_sketch`: greedy that selects on these estimates (recomputed
  each round against the current tree); `add_synthetic` still validates
  feasibility; the caller owns objective evaluation (exact H in benches; a
  sketch/HLL oracle at scale).
- `pair_intersection_count`: the submodule's nshared*m/tmp pair estimator
  with registers empty in both sketches excluded from the equality count.

## Head-to-head vs the index-seeded lane (stage G, bench/grafting_scale_bench)

Probe: one genus, 24 species x 2 genomes, 8 kb background, a cross-species
module block in species 0..11 (data/make_syn_sketch_probe.py). H is the
oracle mean Hamming distance per distinct k-mer (exact KmerClasses used for
evaluation only). p = 13 (8192 registers), b = 16.

| instance | tau | index lane H (B=0 -> 1) | sketch lane H (B=0 -> 1) | sketch est / exact marginal |
|---|---|---|---|---|
| probe, 120 bp module | 16 / 64 | 0.01040 -> 0.00252 | **0.01040 -> 0.00252 (identical)** | 35 / 3024 |
| probe, 600 bp module | 16 / 64 | 0.03895 -> 0.00146 | **0.03895 -> 0.00146 (identical)** | 556 / 14400 |
| GTDB-shaped genus (8 genomes) | 16..512 | 0.00350 -> 0.00012 | 0.00350 (no graft) | - |

Wall time at this scale is comparable (index fill 0.03-0.25 s + greedy ~0 s;
sketch lane 0.01-0.12 s); the point is memory and scaling: the sketch lane
holds n x 2^p x 8 bytes (48 genomes -> 3 MB; 1 M genomes -> 64 GB at p=13,
so p must scale down with n), while the index lane materializes every
distinct k-mer of every genome (here ~50 MB; at 5 Mbp genomes ~8 GB).

Reading the table: where the module clears the background, the sketch lane
names it and its first graft is IDENTICAL to the index lane's - the naming
task is done as well without the index. The gain proxy's magnitude is the
shared-core count, not the marginal (exact marginals include the |P \ A|
spread factor), which is sufficient for ranking.

## Where sketch-only naming diverges from exact naming (measured)

1. **Submodule intersection_size counts both-empty registers as matches.**
   `BBitMinHasher::intersection_size` adds `core_[i] == o.core_[i]` with no
   non-empty guard (its own `jaccard_index` has the guard). Two disjoint
   10k-key sketches at p = 13 read 1847 (true 0; 691/8192 both-empty
   registers). The lane uses `pair_intersection_count`, which guards.
   Cardinality estimates are unaffected. The vendored submodule's
   equality-based path should be fixed upstream or avoided; the round-2/3
   Jaccard lane used `jaccard_index` (guarded, unaffected).
2. **k-way order-statistic estimators fail at genomic resolutions.** The
   all-equal composite needs a shared key to be the bucket minimum of every
   member: a ~300-kmer 24-way intersection estimated as ~3. The max-of-minima
   (Li-Owen-Zhang k-th order statistic) composite misapplied under the
   submodule's harmonic estimator reads ~29k where truth is 0 (disjoint
   pairs at n ~ m), and correctly derived k-min formulas cannot distinguish
   a small intersection from none (max-of-minima distribution is the same):
   measured 2063 for a disjoint pair of 10k-key sketches. Only PAIR
   estimates are used, ever.
3. **Resolution vs background chance-sharing.** Two independent genomes of
   n k-mers each share ~n^2 / 4^k distinct k-mers by chance: ~4e-6 at 8 kb,
   ~0.05 at 30 kb - negligible at k = 17 (round 5 corrected an earlier
   erroneous ~53 figure here; the sketch chance-agreement floor is
   m x 2^-b ~ 0.125 registers). The 30 kb probe failure was therefore NOT
   background: it was REGISTER RESOLUTION - with n/m ~ 3.7 at p = 13 a
   104-kmer module survives as bucket minima in only ~8 registers and the
   min-over-pairs gain collapses to 0 (measured: candidate est 0 at p = 13,
   est 78 at p = 16 - doc/grafting-two-tier.md). The register budget must
   scale so n/m stays near or below 1: p = 13 covers ~8 kb genomes, 30 kb
   needs p = 16. Modules below the resolution floor are dropped by the
   lane (the exact index sees them because classes are exact).
4. **Threshold scale dependence.** A raw pair-intersection threshold chains
   single-linkage through shared clade background: on the GTDB-shaped set
   (a genus block shared by every species) every pair clears any tau below
   the background, one whole-clade component forms, and no candidate is
   emitted. A differential above the background (the {sp0,sp1} class, 400
   k-mers over a 484-kmer background) separates only at tau in (484, 884]
   and the probe's 104-kmer module needs tau <= ~100: no absolute tau fits
   both, and the differential case still failed at the grid's top end in
   the measurement above. Adaptive tau (e.g. per-node background estimate
   plus a margin) is the natural next lever.
5. **Estimator bias/variance statement (p = 13, b = 16, pair estimator,
   guarded):** disjoint 10k/10k: 0 (exact after the guard; +1847 with the
   submodule's unguarded path); overlapping 8k of 10k: 8075 (true 8000,
   +0.9%). Small intersections carry relative error growing roughly as
   sqrt(true x m)/m per estimate; the gain proxy's min-over-pairs is the
   pessimistic tail over |parts| choose 2 estimates, and contamination is
   the max - so gains are biased DOWN and near-threshold modules are
   dropped (see 3).

## Test status

137-case suite plus the round-3 `graft_high_fanout_union_lane` gate; the
sketch-only lane is covered by the stage-G head-to-head (bench) and
exercised in `graft_sketch_only_naming` (unit), which asserts the lane
recovers the planted module to oracle delta 0.
