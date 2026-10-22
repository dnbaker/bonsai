# Grafting round 5: two-tier per-node sketches + adaptive tau (design note)

PI directives (swarm board 2026-10-21) on top of the round-4 sketch-only
lane (doc/sketch-only-naming.md):

1. **Two-tier per-node sketches.** Small nodes (clade distinct-k-mer count
   n <= EXACT_MAX = 512) carry an EXACT cover (a sorted hash set of the
   clade's k-mers; zero estimator error). Above the exact capacity the node
   falls back to bbminh; the MIDDLE band (512 < n < DENSE_MIN = 25,000)
   keeps a PARTIALLY-FILLED bbminh (expected fill 1 - e^{-n/m}, m = 8192
   registers at p = 13: ~6% fill at n = 512, ~95% at n = 25k). Nodes with
   n >= DENSE_MIN hold a densified sketch. Tier decisions need per-node
   clade cardinalities: production reads them from the persisted BNSKMC01
   sidecar (`CompactDatabase::taxon_counts()`, folded up the tree); the
   lane computes the same numbers sketch-only from register fill
   (linear counting, n_hat = -m ln(V/m)) when no sidecar is at hand.

2. **Adaptive tau.** Round 4's absolute threshold could not fit both the
   8 kb background regime (chance sharing ~2 k-mers) and the 30 kb regime
   (~53): the probe's 104-kmer module was named at 8 kb and lost in the
   noise floor at 30 kb. The acceptance threshold is now per node:
       background_v = median over child pairs (i,j) of est|A_i n A_j|
       sigma_v    = sqrt(background_v)          (matching-register counting noise)
       tau_v      = background_v + k * sigma_v,  k = 4 (floor TAU_FLOOR = 32
                    when background_v ~ 0)
   The background enters every pair, so the median cancels it; what
   survives the margin is pair-private (differential) sharing. Predicted:
   8 kb probe tau ~ 32 + eps (module 104 clears); 30 kb probe tau ~ 53 +
   4*7.3 ~ 82 (104 clears, margin 22); the GTDB-shaped genus with a 484-
   kmer background: tau ~ 484 + 4*22 = 572 (the {sp0,sp1} differential at
   884 clears; the 484-background pairs do not). Both round-4 failure
   cases should flip to named. Edge formation uses tau_v; gain acceptance
   in graft_greedy_sketch uses the same per-node threshold (an optional
   threshold for the index lane exists for parity experiments).

3. **Union rule (test-resolved in the sketch repo).** Merge RAW registers
   (register-wise min; the empty sentinel uint64 max = +inf) and densify
   ONCE, at estimation time, on a COPY, only for nodes whose tier says so.
   NEVER union pre-densified sparse sketches: min(dens(A), dens(B)) !=
   densify(A ∪ B) (0.179 vs 0.235 at 30% fill in the sketch repo's
   regression test). The lane's clade sketches are always raw unions from
   per-genome sketches, so the rule holds by construction; the densifier
   is `BBitMinHasher::densify` / `detail::densifybin` - ALREADY VENDORED
   at the current pin (no pin bump needed; maint-2's 3922e43 formalizes
   the API and stays a follow-up).

**Estimator bias/variance, honestly (bbminh pair estimator, p = 13, b = 16).**
Guarded pair intersection (both-empty registers excluded from the equality
count): disjoint 10k/10k keys -> 0; overlap 8k of 10k -> 8,075 (true 8,000,
+0.9%). The estimate counts matching registers ~ Binomial(m, p), so the
counting std is ~ sqrt(n_hat_shared); the gain proxy takes a MIN over
|P| choose 2 pair estimates (pessimistic tail) and a MAX for contamination
(optimistic tail), biasing gains down - that is the direction that drops
borderline candidates, and exactly why tau must scale with the background
rather than sit absolute. Plain (partially-filled, linear counting on
empty registers) is 2-8x sharper at small n than pre-densified - hence
per-node choice: plain in the middle band, densified copies only for
n >= 25k (where plain and densified agree anyway: 95%+ fill). b = 16 keeps
the chance-agreement floor (m x 2^-b = 0.125 registers) below the noise
for any tau >= TAU_FLOOR. Densification itself copies a neighboring
register value into each empty register (the sketch repo's backward-pass
densifier): it removes the linear-counting sparsity correction at the
price of a small downward bias, and it is IDEMPOTENT only on raw unions -
hence rule 3.

**Exact covers.** Round-5 plumbs the tier but the cover BUILD remains
production wiring: at DB-build time the builder streams every genome's
k-mers anyway, so materializing <= 512 hashes per small clade is free
there, and the sidecar already persists the per-taxon counts that decide
the tier. The lane accepts an optional cover provider; unit tests
exercise the exact tier with hand-built covers. Without a cover, an
exact-tier node falls back to its sketch (documented, not silent).

## Stage-H measurement (bench, stage H; probes data/make_syn_sketch_probe.py)

| probe | p | abs16 | adaptive (per-node tau) |
|---|---|---|---|
| 8 kb bg, 120 bp module | 13 | named, H 0.00252 | named, H 0.00252 (tau = floor 32) |
| 8 kb bg, 120 bp module | 16 | named, H 0.00252 | named, H 0.00252 |
| 30 kb bg, 120 bp module | 13 | NOT named (H 0.00556) | NOT named (H 0.00556) |
| 30 kb bg, 120 bp module | 16 | NOT named (see below) | **named, H 0.00307 (round-4 failure flipped)** |

Two mechanisms, both measured. (1) At 30 kb and p = 13 the module's gain
estimate is 0 even for the correct component - register resolution (n/m ~
3.7), not tau: the fix is the register budget rule, p chosen so n/m <= ~1
per sketch (8 kb -> p = 13; 30 kb -> p = 16). The design note's tier
thresholds are absolute k-mer counts; the register budget scales with them.
(2) At 30 kb and p = 16 the ABSOLUTE tau = 16 arm fails differently: noise
edges (pair intersections ~16-30) chain extra species into the component
(mask grows 24 -> 30 genomes), the min-over-pairs core collapses, and the
diluted candidate's estimate drops to 0. The adaptive arm's floor (tau = 32)
rejects those edges, keeps the component at the 12 carrier species, and
names the graft. The floor is doing exactly the work the absolute grid
cannot: regime-appropriate edge admission without a hand-tuned tau.

Densifier: `BBitMinHasher::densify` / `detail::densifybin` already vendored
at the current submodule pin (used on copies for DENSE-tier estimates only);
the maint-2 pin bump (3922e43) formalizing the API remains a follow-up.
Exact covers: tier plumbed via SketchNodePolicy::exact_max; the cover
provider is production wiring (built at DB-build time; the BNSKMC01
sidecar's per-taxon counts folded up the tree give the tier decisions for
free); without a provider, exact-tier nodes use their plain sketch.
