# Real-data grafting test (round 6): the D6331 Escherichia coli pangenome

Date: 2026-10-22, overnight run. Branch kimi/uniq-bbmh.

## Dataset (provenance)

The five ZymoBIOMICS D6331 E. coli strain references, already on disk in
the repo (no downloads): `test/zymo/refs/D6331.refseq/genomes/
Escherichia_coli_{B1109,b2207,B3008,B766,JM109}.fasta` (~4.7-5.2 Mbp each;
chromosome + plasmid records; 2019 RefSeq builds used by the D6331 mock).
A well-known commercial reference set and a genuine mini-pangenome: five
same-species isolates with shared core, strain-specific accessory content,
and plasmids (mobile elements by construction). Fixture builder:
`data/make_ecoli_pangenome.py` (taxonomy chain root -> ... -> species 562
-> five strain leaves 910001-910005); 14 tree nodes, 5 leaves.

## Parameters

k = 17 canonical k-mers; round-5 two-tier policy with the differential-edge
extensions; B = 8 graft budget; the exact KmerClasses oracle is used for
EVALUATION and the k-mer inventory only.

**CORRECTION (round 7): the round-6 runs reported "p = 21" but actually
built p = 13 sketches.** BSD getopt stops at the first positional argument,
so every invocation with flags after the dataset paths silently used the
defaults (p = 13, m = 8192). The sketch-lane numbers below are therefore
p = 13 unless stated; the round-7 arms pass flags before positionals and
rerun at the intended p = 21 / p = 19. The inventory findings are
oracle-based and unaffected. Always put flags first.

## Results

Wall (5 genomes, whole bench process incl. the oracle fill): index lane
fill 4.8 s + greedy ~0 ms; sketch-only lane 0.33 s; peak RSS 722 MB
(oracle-dominated; the sketch lane's own working set is ~260 MB of
sketches).

H floor (mean Hamming error per distinct k-mer, oracle):

| lane | B=0 | B=1 | B=2 | B=4 | B=8 |
|---|---|---|---|---|---|
| index-seeded (exact gains) | 0.4446 | 0.3135 | 0.2592 | 0.2385 | 0.2385 (4 accepted) |
| sketch-only (differential pairs) | 0.4446 | 0.4243 | 0.4135 | 0.4154 | 0.4398 (8 accepted) |

The initial H = 0.44 is the pangenome itself: ~44% of distinct k-mers sit
in multi-strain classes that the species-level tree over-covers (core +
accessory clouds). The index lane removes almost half that error with 4
grafts. The sketch lane, honestly, does NOT improve H on this data: its
accepted grafts' exact marginals are frequently negative (est-vs-exact
pairs like 1,706,130 vs -46,134) because the differential count is
**inflated by core polymorphism** - a core bucket where an outgroup genome's
minimum diverges (a SNP/indel makes a different k-mer the bucket min) is
counted as differential even though the content is core. This is the
real-data divergence the synthetics could not show: the min-hash
differential filter conflates "accessory content" with "core sequence
divergence". Documented as the round-6 limitation; the est-vs-exact table
is in the run log.

## Named-graft inventory (what got named, and does it make sense)

Index lane grafts (tree-leaf masks, best biological reading):
`{JM109, B3008, B766, B1109}` (all but b2207) first - **b2207 is the
outgroup** - then `{JM109, ?}`, `{JM109, B1109}`, `{b2207, ?}`. The
four-strain block is a ~1M-kmer accessory genome shared by everything
except b2207; biologically plausible as a large shared accessory backbone
(virulence/plasmid-derived), expert-light judgement: the four strains
group together against b2207 consistently in every analysis arm.

Sketch-lane grafts include 8 strain pairs; the k-mer inventory
(`BNS_INVENTORY=1`; exclusive content located on the first carrier's
records, spans >= 1 kbp):

| graft (carriers) | exclusive kmers | located spans |
|---|---|---|
| JM109 + B1109 | 51,515 | JM109 chr: 234-235.5 kbp (1.25 kbp); **1619-1635.6 kbp: four clustered islands 1.05-1.44 kbp** |
| JM109 + b2207 | 78,962 | JM109 chr: 1651.4-1652.6 kbp (1.17 kbp) |
| B766 + B3008 | 220,957 | B766 chr: 4138.9-4140 kbp (1.0 kbp) |
| JM109 + B3008, JM109 + B766, B3008 + B1109, B3008 + b2207, b2207 + B1109 | 10-64 k | no >= 1 kbp span on the first carrier (content on the other carrier or fragmented) |

The 1.0-1.4 kbp clustered islands at specific loci (e.g. four islands
within ~17 kbp at ~1.62-1.64 Mbp of JM109 shared with B1109) are the
size and clustering signature of **mobile elements** - prophage remnants,
insertion-sequence clusters, or genomic islands. Judgement, honestly
labeled: no functional annotation was run offline, so this is
expert-light plausibility (size, clustering, strain-restriction), not
annotation. The plasmid records themselves are strain-unique here
(no cross-strain plasmid sharing survived as a graft at k = 17).

## Resolution check (the n/m rule on real data)

At p = 21 (n/m ~ 2.4) the large accessory blocks (10^5-10^6 kmers) are
fully resolved; the 1 kbp islands (10^3 kmers) survive with margin
(10^3 x (m/n)^2 ~ 180 registers). The round-5 resolution floor was not
the binding constraint at this genome size; the binding constraint is the
core-polymorphism estimator bias above.

## Honest limitations

1. The sketch-only lane's differential estimator cannot distinguish
   accessory content from core sequence divergence on real pangenomes;
   it names real pairs (the inventory above) but its gain magnitudes are
   unusable for greedy acceptance (H oscillates). The exact index lane is
   the quality reference on real data today; the sketch lane is the
   cost reference (0.33 s vs 4.8 s, 260 MB vs 700+ MB) with named-candidate
   recall. Closing the bias (e.g. divergence-robust filters, or anchoring
   the differential on k-mer-presence runs rather than per-bucket minima)
   is the next lever.
2. Adaptive tau in differential mode is a multiplicative noise band
   (max(floor, k x sqrt(median))); at pangenome scale the median
   differential is signal, not background - the additive form would
   reject everything.
3. Five strains is a small pangenome; the 64-leaf PlacementTree cap
   bounds any expansion - a LightTaxTree-based PlacementTree would be
   needed for tens of strains.

## Reproduce

```
/usr/bin/python3 data/make_ecoli_pangenome.py
BNS_INVENTORY=1 ./grafting_scale_bench data/ecoli_pangenome/taxmap.txt \
    data/ecoli_pangenome/names.txt data/ecoli_pangenome/genomes.list \
    -k 17 -s 2097152 -p 2 -N 64
```

## Round 7: divergence-robust acceptance on the Zymo pangenome

PI mandate: make sketch-based ACCEPTANCE divergence-robust. Three
mechanisms implemented and measured (each a separate `BNS_MODE` arm; exact
oracle for evaluation only):

1. **Exact-cover regime switch (`covers`).** Per-node FULL k-mer sets
   (full = union, core = intersection; ExactCovers in grafting.h) make the
   differential count exact: |inter parts| - |core(parent)|, error 0 by
   construction. Cost on Zymo: ~415 MB for 18 nodes (5 x 40 MB leaf sets +
   184 MB species union + 30 MB species core; unary-chain nodes alias).
   Memory per mode on this dataset: sketch lane ~260 MB (p=21), exact
   covers ~415 MB, exact index (KmerClasses) ~700 MB. BUT: the measured
   arm accepts 8 grafts whose exact marginals are STILL often negative
   (est 2.08M vs exact 86k; three negatives) - the exact COUNT is not the
   MARGINAL. Round-7 finding: exact covers fix the estimator, not the
   acceptance metric. The small-pangenome operating point is the INDEX
   lane (exact class marginals), which needs no new machinery.
2. **Replication across register budgets (`replication`).** Independent
   builds at p=19/16/13 (register folding is INVALID for this: registers
   hold 64-p value bits, so a fold loses the top bits - unit-tested): all
   10 candidates' differentials agree within ~1% (1.72M/1.71M/1.70M). The
   agreement filter passes everything; accepted grafts reproduce round
   6's sign failures exactly (est 1.85M vs exact +154k; 3 of 8 exact
   marginals negative; H oscillates 0.42-0.44).
3. **Multi-seed CV (`multiseed`, R=3, p=19).** Across-seed CV is 0.00-0.01
   - the inflation is deterministic across hash seeds, so the CV filter
   cannot reject it either.

**The honest round-7 verdict:** the round-6 core-polymorphism inflation is
DETERMINISTIC - stable across register budgets (p=13..19), across hash
seeds (CV ~0.01), and reproduced by the exact-count computation - so no
variance/replication filter can remove it, and the exact-count variant
proves the problem is the acceptance METRIC (differential count != Hamming
marginal), not the estimator. The mandate's (b) outgroup-minima Hamming
guard is infeasible on Wang-hashed registers by design: a SNP changes the
hash wholesale (avalanche), so no Hamming-1/2 neighborhood of the bucket
minimum exists to guard; the guard becomes implementable only under the
position-encoding / masked-reversible hash family (the PI's WangHash
reversibility line), where a SNP shifts the encoded word by a known
amount. That hash swap is the concrete next lever this measurement points
at.

**Acceptance table (Zymo, after = replication arm; before = round 6):**
| accepted grafts | sign agreement | example est / exact marginal |
|---|---|---|
| round 6 sketch (before) | 5/8 | 1,706,130 / -46,134 |
| round 7 replication (after) | 5/8 | 1,850,741 / +154,545 (next: 1,834,720 / +82,341; then -46,134) |
| round 7 exact covers | 5/8 | 2,077,443 / +86,787 (next: -200,007) |

H floors (B=8): index lane 0.4446 -> 0.2385; every sketch-count arm
(round 6/7) oscillates at 0.41-0.44 - no sketch-count acceptance variant
improves H on this data. Recommended operating point, confirmed: small
pangenomes run the exact index lane; sketch-approximate naming serves
candidate RECALL at metagenomic scale with exact verification of the top
candidates (the round-4 hybrid), and the acceptance-metric fix (marginal-
shaped sketch estimators, or the position-encoding hash that makes the
outgroup guard real) is the research next step.

Wall: covers arm greedy 10.5 s (exact set intersections per round);
replication arm 14.1 s (three extra sketch builds at p=19/16/13);
multiseed ~3 x sketch build; all arms' sketch sets ~50-260 MB.
