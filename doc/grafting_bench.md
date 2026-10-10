# Grafting bench: H vs budget, baselines, exact gap, timing

Date: 2026-10-10. Machine: x86-64, single thread, `g++-15 -O3 -march=native`.
Bench: `bench/grafting_bench.cpp` (`make grafting_bench`). Both fixtures:
32 genomes (16 species x 2 genomes, 2 genera), 1300 bp random background,
60 bp species blocks, 120 bp module blocks, k = 17 canonical k-mers.
H = total_delta / total_distinct (mean Hamming distance per distinct
k-mer); fn/fp are the false-negative / false-positive decompositions of
total_delta.

## Fixture A: cross-genus modules (mA spans species 0-3 and 8-11)

Initial: 51 classes, 42104 distinct k-mers, delta 7072 (H = 0.16797).
Candidate lanes:

| lane | B=0 | B=1 | B=2 | B=3 |
|---|---|---|---|---|
| no grafts | 0.16797 (fp 7072) | | | |
| feasible union/split family (492) | 0.16797 | no grafts taken | | |
| observed classes (19) | 0.16797 | 0.09880 (fp 4160) | 0.03952 (fp 1664) | 0.00000 |
| BBMH tree-layered sketch (7) | 0.16797 | no grafts taken | | |
| random 3 grafts (20 seeds) | | | | mean delta 6614 |

The feasible family cannot represent cross-taxon masks; the sketch lane
misses mA because cross-genus agreement (pair J ~ 0.036) dilutes to
J ~ 0.002 at the root-layer clade comparison. fn stays 0 throughout: all
baseline error is false positives (classes sit strictly inside implied
clades).

## Fixture B: genus-local modules (mA = species 0-3, mB = species 4-7,
## mC = species 8-9)

Initial: delta 2912 (H = 0.06916). All three lanes are identical:

| lane | B=0 | B=1 | B=2 | B=3 |
|---|---|---|---|---|
| feasible family / observed classes / sketch lane | 0.06916 | 0.03952 | 0.01976 | 0.00000 |

(random 3 grafts: mean delta 2579; fp follows 2912 -> 1664 -> 832 -> 0.)

## Greedy vs exact optimum (8-leaf module instance, B&B over all subsets)

| B | greedy delta | exact delta | reduction ratio |
|---|---|---|---|
| 1 | 208 | 208 | 1.000 |
| 2 | 104 | 0 | 0.900 |
| 3 | 0 | 0 | 1.000 |

Unit-test instances (test/grafting.cpp): 8-leaf module fixture B=2 ratio
0.833; 6-leaf B=3 and four seeded randomized instances gap-free.

## Wall clock ( Fixture A / B )

| stage | A | B |
|---|---|---|
| exhaustive k-mer index fill | 0.094 s | 0.112 s |
| feasible-family enumeration | 0.0002 s | 0.0002 s |
| observed-class candidates | ~0 s | ~0 s |
| BBMH sketch build (32 genomes) | 0.012 s | 0.010 s |
| tree-layered candidates | 0.047 s (7) | 0.052 s (9) |
| greedy scoring (observed-class lane) | 0.001 s | 0.001 s |
| greedy scoring (feasible family, 492 cands) | 0.009 s | 0.032 s |

Where the objective-evaluation time goes at this scale: the exhaustive
index fill dominates everything else combined (~10x); candidate generation
(sketch build + per-node merges/finalizes along the 63-node tree) is next;
gain evaluation itself is milliseconds (incremental, affected-classes
only). At real scale the index fill is precisely what cannot be afforded,
which is the sketch lane's role: register-only candidate discovery, with
exact verification restricted to finalists. The sketch within-candidate
ranking estimate (inclusion-exclusion on b = 16 register cardinalities)
overlaps the exact top-5 first-step gains by 3-4 of 25 ordered pairs on
these fixtures; it is an ordering heuristic only - the greedy verifies
every accepted graft exactly.
