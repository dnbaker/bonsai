# Spaced-seed design (Monte Carlo sensitivity; codon rates are an assumption)
| weight | pattern | name | P(>=3 hits) uniform 2% / 5% | P(>=3 hits) codon 2% / 5% | mean hits codon 3% |
|---|---|---|---|---|---|
| 21 | `111111111111111111111` | contiguous | 1.000 / 0.992 | 1.000 / 0.994 | 67.6 |
| 21 | `10011111111000111100101011010111` | uniform_opt | 1.000 / 1.000 | 1.000 / 1.000 | 62.0 |
| 21 | `10011111111000111100101011010111` | codon_opt | 1.000 / 1.000 | 1.000 / 1.000 | 62.0 |
| 25 | `1111111111111111111111111` | contiguous | 0.999 / 0.964 | 1.000 / 0.960 | 57.8 |
| 25 | `11110111011101011111110011110111` | uniform_opt | 1.000 / 0.995 | 1.000 / 0.995 | 54.7 |
| 25 | `11101111011101011111110111100111` | codon_opt | 1.000 / 0.994 | 1.000 / 0.993 | 54.6 |
- Next: ~~validate on real held-out divergence (podar strain pairs, mock E. coli strains), after the spaced-seed strand check passes.~~ Done 2026-10-21 (seeds_validate, apps/results/seeds_validate.md): after the both-strand canonical fix (1d810a6) and the spaced-DB regeneration on kimi/uniq-bbmh, the strand check passes on all four rebuilt spaced DBs (fwd == rev, 99.9% of reads to species), and on real strain pairs the optimized spaced patterns match or slightly beat contiguous k21/k25 at equal weight (podar S. parasanguinis pair: 92.3% vs 92.0% correct strain leaf; D6331 E. coli 5 strains: 15.1% vs 14.7%), 100% classified, strand-symmetric on every configuration.
