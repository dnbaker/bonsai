# Spaced-seed design (Monte Carlo sensitivity; codon rates are an assumption)
| weight | pattern | name | P(>=3 hits) uniform 2% / 5% | P(>=3 hits) codon 2% / 5% | mean hits codon 3% |
|---|---|---|---|---|---|
| 21 | `111111111111111111111` | contiguous | 1.000 / 0.992 | 1.000 / 0.994 | 67.6 |
| 21 | `10011111111000111100101011010111` | uniform_opt | 1.000 / 1.000 | 1.000 / 1.000 | 62.0 |
| 21 | `10011111111000111100101011010111` | codon_opt | 1.000 / 1.000 | 1.000 / 1.000 | 62.0 |
| 25 | `1111111111111111111111111` | contiguous | 0.999 / 0.964 | 1.000 / 0.960 | 57.8 |
| 25 | `11110111011101011111110011110111` | uniform_opt | 1.000 / 0.995 | 1.000 / 0.995 | 54.7 |
| 25 | `11101111011101011111110111100111` | codon_opt | 1.000 / 0.994 | 1.000 / 0.993 | 54.6 |
- Next: validate on real held-out divergence (podar strain pairs, mock E. coli strains), after the spaced-seed strand check passes.
