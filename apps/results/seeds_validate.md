# Spaced seeds on real strain pairs (seeds.md's pending validation)
| strain set | pattern | reads to correct strain leaf | to species | classified | strand check (fwd vs rev leaf frac) |
|---|---|---|---|---|---|
| podar | k21_contig | 0.920 | 1.000 | 1.000 (unc 0.000) | 0.919 vs 0.922 OK |
| podar | k21_spaced (spaced) | 0.923 | 1.000 | 1.000 (unc 0.000) | 0.923 vs 0.924 OK |
| podar | k25_contig | 0.917 | 1.000 | 1.000 (unc 0.000) | 0.918 vs 0.920 OK |
| podar | k25_spaced (spaced) | 0.923 | 1.000 | 1.000 (unc 0.000) | 0.922 vs 0.923 OK |
| ecoli | k21_contig | 0.147 | 1.000 | 1.000 (unc 0.000) | 0.146 vs 0.146 OK |
| ecoli | k21_spaced (spaced) | 0.150 | 1.000 | 1.000 (unc 0.000) | 0.150 vs 0.151 OK |
| ecoli | k25_contig | 0.147 | 1.000 | 1.000 (unc 0.000) | 0.147 vs 0.147 OK |
| ecoli | k25_spaced (spaced) | 0.151 | 1.000 | 1.000 (unc 0.000) | 0.150 vs 0.151 OK |
- 10000 simulated 150 bp reads per run at 1% error from the first strain; the database holds every strain as a leaf.
- overall: PASS (strand symmetry must hold for every configuration)
