# Database configuration sweep (k x w; D6331 mock database, 27 genomes)
| k (weight) | w / span | build GB | db MB | sim classified | sim BC to truth | divergent sp. recovery 0/1/2/5% | ill classified | ont classified | hifi classified | containment AUROC true vs host | BC ill vs hifi | strain-leaf frac | reads to absent relatives | strand check fwd / rev to species (k-mer hits) |
|---|---|---|---|---|---|---|---|---|---|---|---|---|---|---|
| 21 | 21 (all) | 1.21 | 487 | 1.000 | 0.062 | 1.00/1.00/1.00/0.99 | 0.991 | 1.000 | 0.999 | 0.990 | 0.118 | 0.237 | 176 | 0.999 / 0.999 (hits 418326 / 418344) |
| 25 | 25 (all) | 1.51 | 621 | 1.000 | 0.063 | 1.00/1.00/1.00/0.95 | 0.991 | 1.000 | 0.996 | 0.970 | 0.117 | 0.238 | 122 | 0.999 / 0.999 (hits 389767 / 389790) |
| 31 | 31 (all) | 1.91 | 822 | 1.000 | 0.064 | 1.00/1.00/0.98/0.82 | 0.990 | 1.000 | 0.996 | 0.933 | 0.117 | 0.237 | 78 | 0.999 / 0.999 (hits 349588 / 349619) |
| 21 | span 27 (spaced) | 1.21 | 487 | 0.670 | 0.105 | 0.52/0.53/0.52/0.52 | 0.617 | 1.000 | 0.952 | 1.000 | 0.148 | 0.511 | 3117 | 0.999 / 0.025 (hits 399115 / 2025) |
| 25 | span 31 (spaced) | 1.51 | 621 | 0.650 | 0.104 | 0.51/0.52/0.51/0.50 | 0.603 | 1.000 | 0.812 | 0.972 | 0.128 | 0.508 | 2192 | 0.999 / 0.016 (hits 371353 / 1628) |
| 25 | span 35 (spaced) | 1.51 | 621 | 0.650 | 0.105 | 0.51/0.51/0.51/0.50 | 0.603 | 1.000 | 0.809 | 1.000 | 0.131 | 0.509 | 2308 | 0.999 / 0.015 (hits 359044 / 1524) |
| 31 | span 37 (spaced) | 1.91 | 822 | 0.642 | 0.100 | 0.51/0.51/0.50/0.48 | 0.599 | 0.999 | 0.737 | 0.972 | 0.114 | 0.506 | 1632 | 0.999 / 0.010 (hits 332318 / 1238) |
- divergent sp. recovery: fraction of the mutated B. fragilis reads assigned to species 817 (truth 8000 reads).
- the read cache in test/sweep_cache is temporary; delete after the sweep.
