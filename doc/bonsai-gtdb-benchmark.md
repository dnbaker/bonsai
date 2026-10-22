# Grafting at GTDB shapes: scale bench + candidate-scan assessment

Date: 2026-10-21. Branch kimi/uniq-bbmh. Machine: the PI's Intel Mac under
load (two fits running); all grafting numbers are single-instrument, niced,
and repeatable in seconds. Data: `data/make_syn_scale_data.py` (GTDB-shaped
synthetic hierarchy: species blocks + genus/family modules, so observed
k-mer classes are exactly the planted sets) at two shapes:

- `syn_scale_main`: 2 dom x 4 phy x 2 cls x 2 ord x 2 fam x 2 gen x 4 spe x
  2 genomes = 1024 genomes (60 kb each; fanouts <= 4).
- `syn_scale_fanout`: 2 dom x 2 phy x 16 spe x 4 genomes = 256 genomes with
  genus fanout 16 (PlacementTree max fanout 8 after its internal bucketing).
- stage F probes: single genus, fanout 16 / 24 / 40, one genome per species.

Bench: `bench/grafting_scale_bench` (`make grafting_scale_bench`), stages
per the header comment; stage E adds the feasible-family lane with
scan-cost counters, stage F isolates high-fanout family behavior.
`grafting.h` gained measurement-only counters (scanned / bound-skipped /
exact-evals / per-round seconds) - no algorithmic change.

## 1. Candidate-scan cost: not all-pairs, but exponential in fanout

Greedy scan cost per round = O(candidates x affected-class walk), with the
admissible `gain_bound` filtering candidates that cannot beat the running
best before paying a tree copy. Measured:

| instance | candidates | rounds | scanned | bound-skipped | exact evals | greedy wall |
|---|---|---|---|---|---|---|
| main stage D (8 genomes, 20 classes) | 12 observed + 0 sketch | 2 | n/a | n/a | n/a | 0.00 s |
| main stage E family | 10 | 2 | 27 | 1 | 26 | 0.00 s |
| fanout stage D (64 genomes, 175 classes) | 143 observed + 3 sketch | 8 | n/a | n/a | n/a | 0.00 s |
| fanout stage E family | 326 | 8 | 2,580 | 60 | 2,520 | 0.00 s (8 rounds x ~0.6 ms) |
| stage F fanout=16 | **65,518** | 4 | 131,035 | 65,850 | 65,185 | **0.04 s** |
| stage F fanout=24 | 24 (splits only) | 4 | 47 | 0 | 47 | 0.00 s |
| stage F fanout=40 | 40 (splits only) | 4 | 79 | 0 | 79 | 0.00 s |

Readings:
- The per-candidate exact evaluation dominates (bound filter passes ~50%:
  it is loose whenever few classes carry the mass), but a candidate costs
  one small tree copy + inverted-list walk, so even 65k candidates x 4
  rounds run in 0.04 s. At these class counts the scan is NOT the
  bottleneck and is nowhere near all-pairs in the class dimension.
- The candidate COUNT is the exponential: `enumerate_grafts` unions are
  2^fanout - 2 per node. At fanout 16 that is 65,518 candidates (fine);
  the k <= 20 cap then SILENTLY DROPS THE ENTIRE UNION LANE: fanout 24 and
  40 enumerate only the fanout-many split grafts.
- PlacementTree also restructures internally: the 64-genome genus (16
  species x 4 genomes) comes out with max fanout 8, which is why stage E
  on the fanout set sees only 326 candidates. So on real inputs the union
  family survives only where the placement tree's effective fanout stays
  <= 20.

## 2. Objective-quality curve vs budget (H = mean Hamming per distinct k-mer)

fanout set, genus sub-clade (64 genomes, 1,921,869 distinct, 175 classes):

| B | observed-class lane | family lane |
|---|---|---|
| 0 | 0.00596 | 0.00596 |
| 1 | 0.00238 | 0.00238 |
| 2 | 0.00208 | 0.00199 |
| 4 | 0.00206 | 0.00189 |
| 8 | 0.00205 | 0.00182 |

Both lanes plateau at H ~ 0.002 after B = 1-2: the planted modules are
absorbed by the first grafts; the residual is the clean-species
background, irreducible under this candidate set. Error is exclusively
false positives at baseline (classes sit strictly inside implied clades);
false negatives appear only through stealing (none observed on these
fixtures). The family lane grinds below the observed-class lane because
its smaller unions shave clade over-coverage the observed masks cannot
represent - the same relationship as the 32-genome fixtures in
doc/grafting_bench.md.

main set (1024 genomes): stage A exhaustive per-genome fill projects
0.002 s/genome at 60 kb (linear extrapolation to a 5 Mbp bacterium:
~0.2 s/genome; the sketch lane's whole point). Sketch build over 1024
genomes: 0.1 s wall. Tree-layered candidates: 292 components at J >=
0.002 (0.1 s). Dilution profile (sibling-clade Jaccard by depth): cross-
clade agreement is a depth-6/7 signal (genus/species layers; max 0.012)
and dilutes to <= 0.002 above it - confirming the sketch lane is a
within-genus tool, as designed.

## 3. Guarantee checks (re-run on this branch)

From `unit "*graft*"` (test/grafting.cpp, 4 cases, 815 assertions, all
pass): greedy-vs-exact (branch-and-bound) reduction ratios - B=1: 1.000;
B=2: 0.833 on the 8-leaf module fixture (the exact gap is the stealing
destruction geometry from the counterexample); B=3: 1.000; gap-free on 6
leaves and all 4 seeded randomized instances. The best-fit
non-submodularity counterexample is locked in
`graft_bestfit_not_submodular`; the fixed-slot (1 - 1/e) proof and the
laminar-restriction argument stand (doc/grafting-objective.md).

## 4. Robustness fix found by the bench

PlacementTree's ancestor-closure loop did not treat a self-parented root
(kraken `1 | 1 |` convention) as a root: it pushed the taxon onto its own
pending list forever and grew memory without bound (found when stage F's
synthetic taxmap hung the bench at 4 GB RSS). One-line guard in
placement.h: `par == t` now counts as `par == 0`.

## 5. Round 3: the polynomial union lane (implemented) - what collapses the floor and what does not

The round-2 proposal was pairwise/triplet child unions above the cap. It is
implemented (`enumerate_grafts`, `EXACT_UNION_MAX = 20`) and measured, and
it does NOT collapse the floor - the seeded lane added alongside it does.
Before/after on the stage-F probes (delta units; H = delta / total):

| fanout | lane | candidates | B | delta before (round 2) | delta now | wall |
|---|---|---|---|---|---|---|
| 16 | family (exact unions) | 65,518 | 4 | 40000 -> 0 | 40000 -> 0 | 0.05 s |
| 24 | family (round 2: splits only) | 24 | 4 | **60000 -> 55000** | - | - |
| 24 | family (pairs+triplets+splits) | 2,324 | 4/16 | - | **60000 -> 55000** | 0.00 s |
| 24 | seeded hulls | 1 | 4 | - | **60000 -> 0** | 0.00 s |
| 24 | observed classes | 25 | 4 | - | **60000 -> 0** | 0.00 s |
| 40 | family (round 2: splits only) | 40 | 4 | **100000 -> 95000** | - | - |
| 40 | family (pairs+triplets+splits) | 10,700 | 4/16 | - | **100000 -> 95000** | 0.02 s |
| 40 | seeded hulls | 1 | 4 | - | **100000 -> 0** | 0.00 s |
| 40 | observed classes | 41 | 4 | - | **100000 -> 0** | 0.00 s |

Diagnosis (why the polynomial family lane cannot move H): a class reports
at the deepest node CONTAINING it. A module spanning k/2 species needs a
single node containing all k/2 clades. Pair/triplet unions never contain
it, so they gain nothing; the one split that drops a non-member species
from the reporting clade gains exactly one leaf's worth (5000 here), and
every remaining statically-enumerated mask then overlaps the new node
(non-laminar, rejected by add_synthetic) - additional splits of the
shrunken clade were never enumerated. Hierarchical composition does not
exist: masks are static and containment does not accumulate. The floor was
never a candidate-COUNT problem; it is a NAMING problem. The lanes that
name the high-arity union directly - the observed class masks themselves
and their laminar child-union hulls (`enumerate_seeded_grafts`, one graft
here at O(|seeds| x fanout)) - collapse it in one accepted graft.

Cost check (gate): the polynomial family's per-candidate cost stays in the
measured class - 21,399 scans + 18,739 exact evals in 0.02 s at fanout 40
(~1 us each); 65,518 candidates scan in 0.05 s at fanout 16. Seeded hulls
are ~free. The 256-genome fanout dataset (stage D/E) is unchanged before/
after, honestly: its placement tree caps at fanout 8 (below the cap), so
the cliff does not bite there; the probes are the clean testbed.

Guarantee checks stay green (137 cases, incl. the new
graft_high_fanout_union_lane: exact family counts pairs+triplets+splits,
seeded/observed recover to delta 0 at fanout 24 and 40, family lane
measured at exactly initial - 5000 with one accepted graft; and the round-2
greedy-vs-exact ratios unchanged).

## 6. Recommended next lever

1. **Adopt seeded hulls as a standing candidate lane** (the fix, measured):
   always include `enumerate_seeded_grafts(tree, class_masks)` in the
   candidate set. Cost ~|classes| x fanout, one graft names what pairs
   cannot compose.
2. **Scale-out without the k-mer index**: when classes cannot be observed
   (the whole point of the sketch lane), seed the hulls from sketch-cluster
   component masks instead of observed classes - same construction, seeds
   = the components the sketch lane already produces. This is the path that
   keeps the fanout>20 fix in the sketch-only regime.
3. **Pair/triplet lane: keep, but as a sub-module lane** - it expresses
   genuine 2-3-clade structure at negligible cost (0.02 s at fanout 40);
   it is not the high-arity fix and should not be billed as one.
4. Not taken (documented alternative): re-enumerating splits against the
   mutated tree each greedy round would let splits compound (12 rounds to
   zero at fanout 24) but changes the candidate contract and costs a round
   trip per graft; the seeded lane reaches zero in one round cheaper.

(The round-2 text of this section - the proposal - is preserved in git
history; sections 1-4 and the guarantee program stand as measured.)
