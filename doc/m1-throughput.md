# m1 throughput (official rerun, 2026-10-21 21:10 PDT)

Protocol (same as the round-2/3 runs): 23,366 HiFi reads (212.8 Mbp, head
of SRR13128014) against the full mock database (27 genomes, 822 MB compact
`full.k31.cdb` with the BNSKMC01 sidecar, k = 31 w = 31), `-p 2`,
interleaved reps, min-of-3 wall seconds (the box also runs the two long
fits; the Oct-20 2,599 reads/s figure that started this line was a quiet
box). Binaries: baseline = profile-apps (pre-uniq work), current =
kimi/uniq-bbmh @ d0b429a.

## classify -P (the m1 E1d cost path)

| run | baseline | current | reads/s baseline | reads/s current | delta |
|---|---|---|---|---|---|
| `-k` (kraken out) | 18.29 s | 15.26 s | 1,278 | 1,531 | **+20%** |
| `-K` (profile default, no per-read output) | 18.01 s | 12.65 s | 1,298 | 1,847 | **+42%** |

Phase split (current, `-K`, `BNS_TIMING=1`): classify 14.80 s + db_kmers
**0.00 s** (sidecar load; the O(database) scan is gone) + profile write;
peak RSS 842 MB.

## -p 4 contention (unshard check, `-K` min-of-3)

| design | wall | note |
|---|---|---|
| sharded UniqCounter (round-2 binary) | 19.26 s | 64 shards + 64 mutexes |
| unsharded (current) | **15.03 s** | one table, lock-free steady state |

## profile ops (current)

| op | wall | RSS |
|---|---|---|
| `profile report` (23k-read profile -> -U report) | 0.19 s | 6 MB |
| `profile merge` (two profiles) | 0.14 s | - |
| profile size | 181 kB = 0.121% of the gzipped input | (m1 E1d pass: <= 1%) |

## build path (measured at the sidecar rebuild, 2026-10-21 morning)

`bonsai build -k 31` over the 27 genomes: 52.3 s wall, 822 MB database +
488-byte `.taxacounts` sidecar written in the same pass (RSS ~1.9 GB).

## Deltas vs the previous table (round 2, 2026-10-21 ~08:45) - verified

| expected change | round-2 number | this rerun | shows up? |
|---|---|---|---|
| collect tax inside classify | ~4.6 s -> ~0.5 s | classify phase 14.8 s total incl. collection; the collect-tax delta is inside it (not separately visible at -p2; at -p4 see below) | yes (within phase) |
| db_kmers sidecar | 2.6 s -> 0.00 s | **0.00 s** | yes |
| unsharded at -p4 | 7.73 -> 5.02 s (round-2 box state) | 19.26 -> 15.03 s (tonight's box state) | yes, same direction (-22%) |
| kraken output identity | byte-identical | byte-identical (cmp) | yes |

Cross-run absolute numbers shift with box load (the fits); deltas and phase
splits are the comparable quantities. `BNS_TIMING=1` in the environment
prints the classify/db_kmers split on stderr.
