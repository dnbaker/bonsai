# Sketch A/B: RangeMinHash vs BBMinHash for genome-sketch candidates

Date: 2026-10-10. Machine: x86-64, single thread, `g++-15 -O3 -march=native`.
Bench: `bench/sketch_ab_bench.cpp` (`make sketch_ab_bench`), corpus `test/`
(five real genomes: E. coli K12 ~4.6 Mbp, three ~0.7–1.5 Mbp genomes,
phiX 5.4 kbp), k = 17 canonical k-mers, the same `sketch::WangHash` the
sketchers apply internally; exact Jaccard on the hashed k-mer sets is ground
truth.

## Configurations

| config | registers | bits/register (b) | live bytes | final bytes |
|---|---|---|---|---|
| RMM-2048 | 2048 | 64 (full minima) | 16384 | 16384 |
| BBMH-2048 | 2048 | 16 | 16384 | 4096 |
| BBMH-8192 | 8192 | 16 | 65536 | 16384 |

RMM-2048 and BBMH-8192 hold the same total register content (2048 x 64b =
8192 x 16b = 128 Kbit); BBMH-2048 matches RMM-2048's register count and live
bytes at one quarter of the stored content.

`BBitMinHasher<uint64_t>(p, b=16)`: 64-bit register words, 16 retained
quotient bits per register. The b choice is the precision-critical part of
the instantiation. The bias-corrected b-bit estimator subtracts the
chance-agreement floor 2^-b, and distinct k-mer pairs whose quotients collide
add ~n^2 x 2^-b / 2^p of spurious agreement. At b <= 5 (forced when the
register word itself is 16 bits, since p + b <= 16) both terms dominate at
low Jaccard and destroy the tested thresholds; at b = 16 both are negligible
at genome scale (verified: H0 sd 1.2e-4 below).

## Part 1 — real genomes

Sketch time includes gz parsing, canonical k-mer rolling hash, updates, and
finalize (min of 5 reps; total over the 5 genomes). Jaccard time is per
final-pair comparison (min over 2000 passes of the 10 pairs). Accuracy is
against exact Jaccard on all 10 pairs.

| config | sketch total (ms) | Jaccard (ns/pair) | MAE | RMSE | max abs err |
|---|---|---|---|---|---|
| RMM-2048 | 750.97 | 14700 | 0.000547 | 0.000858 | 0.001927 |
| BBMH-2048 | 729.31 | 3500 | 0.001117 | 0.002122 | 0.005392 |
| BBMH-8192 | 689.09 | 4700 | 0.000149 | 0.000225 | 0.000468 |

Per-pair detail (J_true, J_est, |err|); the informative pairs are the two
same-species pairs at J ~ 0.022 and the E. coli / Shigella pair at J ~ 0.63:

| pair (J_true) | RMM-2048 | BBMH-2048 | BBMH-8192 |
|---|---|---|---|
| 0.022198 | 0.023438 (0.00124) | 0.023423 (0.00122) | 0.021836 (0.00036) |
| 0.022578 | 0.023926 (0.00135) | 0.026352 (0.00378) | 0.022690 (0.00011) |
| 0.633763 | 0.631836 (0.00193) | 0.639155 (0.00539) | 0.633295 (0.00047) |

The cross-genome pairs at J ~ 3e-4 read as 0.0 from all three sketches (and
as at most 3.5e-4 from BBMH-8192), far below candidate thresholds.

## Part 2 — H0 noise floor

32 iid-random 1500 bp genomes (496 pairs), J_true ~ 1e-4. Distribution of
the J estimate under the null:

| config | mean | sd | max | pairs >= 0.005 | pairs >= 0.02 |
|---|---|---|---|---|---|
| RMM-2048 | 0.00000 | 0.00000 | 0.00000 | 0/496 | 0/496 |
| BBMH-2048 | 0.00001 | 0.00010 | 0.00145 | 0/496 | 0/496 |
| BBMH-8192 | 0.00002 | 0.00012 | 0.00145 | 0/496 | 0/496 |

RMM is exact here because 1500 bp genomes underfill its 2048 registers, so
its bottom-k lists are the full sorted hash sets. BBMH's corrected estimator
carries a small residual noise (max 1.5e-3 over 496 draws) and still clears
the tested thresholds with 3x headroom.

## Part 3 — relative error vs J

8 kbp genomes with injected shared blocks; mean |J_est - J_true| / J_true
over 12 replicates per level. At this size all registers saturate
(~8k distinct k-mers), so this is the register-limited regime.

| target J | RMM-2048 | BBMH-2048 | BBMH-8192 |
|---|---|---|---|
| 0.002 | 0.542 | 0.380 | **0.131** |
| 0.008 | 0.158 | 0.193 | **0.070** |
| 0.02 | 0.139 | 0.091 | **0.045** |
| 0.05 | 0.057 | 0.091 | **0.040** |
| 0.15 | 0.035 | 0.040 | **0.021** |
| 0.5 | 0.017 | 0.024 | **0.011** |

## Verdict (matched memory budget)

At the matched-content point (BBMH-8192 vs RMM-2048, both 128 Kbit of
register content), BBMinHash with b = 16 is **faster and more precise on
every measured axis**: 8% faster sketching end to end (689 ms vs 751 ms over
the five-genome corpus), 3.1x faster pairwise Jaccard (4.7 us vs 14.7 us),
and 3.7x lower mean absolute Jaccard error (1.5e-4 vs 5.5e-4). The accuracy
margin grows as J shrinks (4.2x better at J = 0.002): bottom-k throws away
all but its 2048 smallest hashes, while the b-bit sketch routes every k-mer
into its register, so at equal register content the b-bit estimator uses
strictly more of the input.

At matched register count (BBMH-2048 vs RMM-2048) the ordering inverts on
accuracy (2x worse MAE: 4x less register content) while staying faster on
both timings. That is the expected content-scaling, and it brackets the
matched-content result.

Caveats: (1) the sketch-time delta is dominated by gz parsing and rolling
hash; the update-path advantage shows fully only in the small-genome limit
(phiX: 2.98 ms RMM vs 1.02/0.55 ms BBMH) and in the 3-4x jaccard speedup.
(2) BBMH's unfinalized register file is 4x larger at the matched-content
point (64 KiB vs 16 KiB live); the stored finalized sketches match at 16
KiB. (3) All numbers are single-machine, single-thread, k = 17.
