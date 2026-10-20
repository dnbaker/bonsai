# bonsai profile applications

This directory holds analyses built on bonsai sample profiles (`classify -P`, `bonsai profile ...`). Mock
communities come first, because their truth is known. Real corpora follow, and they are **streamed** (no raw
reads on disk). Each analysis states its question, design, endpoints with pass thresholds, and what a pass or a
fail would mean. Each one writes `apps/results/<id>.json` (all numbers, the seed and the parameters) and
`apps/results/<id>.md` (verdicts).

```
PY=/Users/dnb13/miniconda3/envs/jtopo-311/bin/python3.11    # stdlib only; any python >= 3.9 works
$PY apps/run.py setup          # mock taxonomy + databases from test/zymo, test/entero (local data)
$PY apps/run.py m1 m2 m3 m4 m5 # or: all
```

Reads are produced by a process (a simulator, `cat` of a local run, `curl` of an ENA URL) and piped straight into
`bonsai classify`. bonsai reads gzip, zstd and plain input transparently. Mocks and real corpora go through the
same code path (`bnsapps/stream.py`).

## What a profile is, and why it should matter

A profile holds, for each taxon, the sample's own reads and k-mer hits, plus a mergeable counter of the distinct
database k-mers its reads hit. That counter is exact up to 512 distinct k-mers and a MinHash or HyperLogLog
sketch beyond that. Three properties are testable claims, not assumptions:

1. **Sufficiency and mergeability.** Profiles of read shards merge to the profile of the whole run, byte for byte.
   The classify -U report, clade sketches and cross-sample statistics then come from profiles alone, with no
   need to revisit reads.
2. **Taxon-resolved k-mer content carries signal that abundance does not.** Examples are false-positive taxa
   (many reads, few distinct k-mers), strain sharing, and strain replacement at equal abundance.
3. **Streaming scale.** A corpus can be profiled while it is downloaded, keeping kilobytes per sample.

**Limitation, in every analysis:** counters see *database* k-mers only. A strain absent from the database
appears only through the database k-mers it shares, that is, as a pattern of k-mer loss. Its novel k-mers are
invisible.

## Mock analyses (local data: Zymo D6331 references, 8 Enterobacterales relatives absent from the mock)

The databases are 19 D6331 genomes (the two yeasts are left out to bound memory; E. coli strains B1109, B2207,
B3008, B766 and JM109 are leaves under species 562), plus `entero/` relatives (Shigella, Klebsiella, Citrobacter,
Enterobacter, Proteus, Morganella, Yersinia, Serratia) that the mock does not contain.

### M1: sufficiency, streaming, cost
- **Data:** the local SRR13128014 HiFi head (about 23k reads), and simulated 150 bp reads.
- **E1a:** profiles of K = 2, 4 and 8 random read shards, merged, equal the whole-run profile byte for byte.
  *Pass:* all identical.
- **E1b:** `profile report` equals `classify -U` byte for byte, for the HyperLogLog and the MinHash profiles.
- **E1c:** the profile of a piped stream equals the profile from the file.
- **E1d (reported):** reads/s, Mbp/s, peak RSS, and profile bytes as a fraction of input bytes.
  *Pass:* profile size ≤ 1% of the gzipped input.
- **Meaning:** a fail in E1a–c means the "profiles are sufficient statistics" claim, and everything built on
  merging, is false.

### M2: false-positive taxa (KrakenUniq's claim, on our counters)
- **Data:** the real HiFi reads, and simulated reads at 1% error. Present taxa are the 15 mock species. Absent
  taxa are the 8 relatives.
- **Scores per species:** clade reads, distinct k-mers, and coverage (distinct k-mers / database k-mers of
  the clade).
- **E2:** AUROC of present vs absent for each score, with bootstrap 95% CI. Also the false positives at the
  threshold that keeps every present species.
- **Expected:** coverage ≥ reads.
- **Meaning:** if coverage does not beat reads, the counters add nothing for presence calls on this data. With
  15 vs 8 taxa the CI is wide; this is a sanity check, not a paper result.
- **Learned from the first run:** a database that holds every organism in the sample produces no false-positive
  taxa. On the real HiFi reads, 0/8 relatives got a read, and every score's AUROC was set by the rarest true
  members. False-positive taxa come from reads of organisms *absent* from the database, hitting conserved k-mers.
  So M2 also runs the mock plus abundant Vibrio and Pseudomonas reads (on disk, not in the database). It runs a
  uniform community and a staggered one, with present species log-uniform in 50–20000 reads. The staggered one
  is the discriminating case.

### M3: strain resolution (5 E. coli strains)
- **E3a (strains in the database):** samples carry 1–2 E. coli strains, at E. coli coverage
  c ∈ {0.1, 0.3, 1, 3}x, on a mock background. The endpoint is the AUROC of "the two samples share ≥ 1 strain"
  against k-mer Jaccard at clade 562. It is compared with the Bray-Curtis of strain-leaf read counts.
- **E3b (strain absent from the database):** leave one strain out of the database (5 folds) and repeat on pairs
  that involve the missing strain. This is the case abundance at strain leaves cannot handle.
- **E3 detection limit:** the lowest coverage with AUROC ≥ 0.95.
- **Meaning:** E3b > 0.9 at ≥ 1x supports strain tracking against incomplete databases. E3b ≈ 0.5 means k-mer
  sharing needs the strain in the database, which is then the limitation to state.

### M4: differential abundance, calibration and power
Replicate mock samples are simulated as reads ~ Multinomial(N, p_i), with p_i ~ Dirichlet(α·p), α = 50, so each
replicate has biological dispersion. Each experiment is 6 vs 6 samples, repeated R = 20 times with seeds logged.
- **E4a null:** both groups come from one distribution. The endpoints are the per-clade FPR at α = 0.05 and a KS
  test of p-value uniformity, for the read test and the k-mer test. A Wilcoxon rank-sum on proportions is the
  baseline.
  *Pass:* FPR ≤ 0.05 + 2 SE.
- **E4b abundance shift:** one species ×f, f ∈ {1.5, 2, 4}. The endpoint is power, the fraction of experiments
  with q < 0.1 for that species.
- **E4c strain swap at equal abundance:** group A carries JM109 and group B carries B766. The endpoints, at rank
  species, are the read-test rejection rate (expected ≈ α) and the k-mer-test power (expected ≥ 0.8).
- **Meaning:** E4c is the headline. A k-mer test that is calibrated (E4a) and powerful where abundance tests
  are blind by construction (E4c) is the strongest argument for keeping taxon-resolved counters.

### M5: nearest-profile search
- **Data:** 20 synthetic "subjects", each with a composition and an E. coli strain set, and 3 samples each
  (resampled, with 20% composition drift).
- **E5:** top-1 same-subject retrieval, and the AUROC of same- vs different-subject pairs, for tree EMD,
  Bray-Curtis, k-mer Jaccard and containment.
- **Expected:** k-mer metrics beat abundance metrics when subjects differ by strain.
- **Meaning:** this is synthetic. The real endpoint is R1.

## Real corpora (streamed; each needs the PI's OK on volume before running)
- **R0 mock runs on ENA:** other D6331/D6300 runs (Illumina, ONT), using a stream head of about 200 MB each.
  The endpoints are M2 on real short reads, and cross-platform profile consistency (Bray-Curtis between
  platforms at species level).
- **R1 longitudinal cohort:** repeated stool samples per subject. The endpoint is same-subject retrieval (M5 on
  real data) by k-mer vs abundance metrics.
- **R2 mother–infant pairs:** the endpoint is strain-sharing calls compared with published transmission.
- **R3 case–control:** the endpoints are DA by reads and by k-mers, null calibration by label permutation, and
  agreement with published markers.
- **Database for R1–R3:** gut species beyond the mock are needed. That is a reference download (size to be
  stated when proposed). Building a database from streamed references needs two passes today: a cardinality
  estimate, then the build.

## What it would take to be a compelling paper
- **Claim:** "Mergeable, taxon-resolved k-mer profiles make metagenomic corpora analyzable without raw reads,
  and expose strain-level differences that abundance cannot."
- **Evidence needed:**
  1. M1 holds on real data at corpus scale (≥ 1000 samples streamed; cost per sample reported).
  2. Calibrated tests (M4a on real data, by label permutation), and power where abundance is blind (M4c, then
     R3 with an orthogonal check such as strain calls from inStrain/StrainPhlAn on a subset).
  3. Retrieval/identity (R1) and strain sharing (R2), each against the obvious tool.
- **Comparisons a reviewer will ask for:**
  - sourmash (FracMinHash sample sketches, containment search, gather). This is the closest prior art: our
    difference is per-clade counters tied to a taxonomy, and tests over them.
  - KrakenUniq and Kraken2/Bracken (distinct k-mers, abundance).
  - StrainPhlAn4 and inStrain (strain sharing).
  - ALDEx2, ANCOM-BC and MaAsLin2 (DA).
- **Ablations:** counter kind and precision (HLL vs MinHash, p), k, and database completeness (holdouts, as in
  M3b).
- **Honest scope:** database k-mers only; small groups limit the attainable p (2 per group: p ≥ 1/3); MinHash
  Jaccard error ≈ 1/sqrt(2^p).
- **Kill criteria:** M4c k-mer power < 0.5, or M4a FPR > 0.1, means no DA paper. R1 showing k-mer ≤ abundance
  means the search claim goes. If both fail, the contribution is engineering (bonsai profiles as a format),
  not a methods paper.

## Status (2026-10-20, first full mock pass; numbers in results/)
- **M1 PASS on real HiFi reads.** Merging 2, 4 or 8 shards reproduces the whole run; `report` equals `-U`; a
  pipe equals the file. The profile is 0.12% of the gzipped input.
- **M2.** With reads of organisms absent from the database, in a staggered community: reads AUROC 0.89
  (4/8 relatives above the rarest true species) vs distinct k-mers / coverage 1.00 (0/8). n = 15 vs 8.
- **M3 partial.**
  - With strains in the database, strain-leaf abundance wins at every coverage.
  - With the strain absent, k-mer Jaccard beats leaf abundance only at 3x (0.96 vs 0.85). At <= 1x the two are
    about level (0.79-0.87 vs 0.81-0.87).
  - A post-hoc coverage-normalized statistic did not help.
  - Low-coverage strain sharing is limited by sampling, not by the statistic.
- **M4 holds after the calibration follow-up.**
  - k-mer test FPR: 0.080 over the first 20 null experiments, then 0.048 +/- 0.006 over 60 more (m4cal; KS
    p 0.93). Pooled 0.056.
  - Strain swap at equal abundance: k-mer power 0.80, read test 0.10.
  - Read test calibrated (0.043 and 0.063 +/- 0.008); power 0.80 at x2.
- **M5 negative for whole-sample k-mer search.** Abundance metrics beat k-mer metrics in both synthetic
  regimes (beta 5: Bray-Curtis top-1 0.92 vs k-mer 0.68). At these depths, whole-sample counters reflect
  resampling. If search is to use k-mers at all, it has to be per clade.

### Counts and sets together (PI: "both info together is stronger than one alone")
- **M2 read saturation.** Score = distinct k-mers / E[distinct | reads], where E = g (1 - exp(-r w / g)).
  - g is the clade's database k-mer count; w is the sample's mean hits per classified read.
  - Staggered + novel: AUROC 1.00, with a margin. Every present species scores >= 0.34 (one at 209 reads); every
    absent relative scores <= 0.20 (0.06-0.09 at 230-270 spurious reads).
  - The score does not grow with abundance, so one threshold should transfer across depths, unlike raw
    distinct counts.
  - Hits are the wrong sampling unit (each read is a run of ~120 distinct k-mers); a hits-based saturation
    score failed (AUROC 0.16-0.98).
- **M3 capture-recapture.** Each sample's clade universe comes from distinct = g (1 - exp(-hits / g)), giving
  the coverage-corrected genome Jaccard I / (gA + gB - I), with I = shared gA gB / (uA uB).
  - Strain absent from the database: best k-mer statistic at 1x (0.888 vs Jaccard 0.861, leaf 0.874).
  - Worse at <= 0.3x, where few duplicates leave g unidentifiable and MinHash error (~3%) dominates.
  - Below Jaccard at 3x.
  - Post hoc, same samples.

### R0: real mock runs streamed from ENA (42 runs, 60 MB R1 heads, ~2.2 GB transferred, ~0 disk)
- **Host background.** Human-only controls classify 1.30e-3 of reads (6-9 species with reads). The D6331
  dilution series sits on that floor from 1:10^3 on.
- **Dilution slope.** After subtracting the host background, the microbial fraction falls ~10x per 10x dilution
  through 1:10^3; the slope over the clear-excess range is in results/r0.md.
- **Coverage >= 1e-3** removes every host artefact (0 species in the human-only controls), but loses true species
  as depth falls (12, then 9, then 2-3 of 15 species from 1:1 to 1:10^3).
- **Blanks.** Template-free extraction blanks have 21-23 species with reads, 8-9 passing coverage: kit
  contamination (the "kitome"), which needs a negative-control model rather than a threshold.
- **Platform bias.** Species Bray-Curtis is 0.07-0.14 within a platform, 0.16 Illumina vs HiFi, and 0.43 ONT
  vs either.
- **D6300** (absent members abundant and distant): AUROC 1.0 for every score.
- **Read saturation with M2's pre-set threshold (0.3) does NOT transfer.**
  - 1:1: 1-3/15 species pass. 1:10^2-10^3: 9-11/15 pass. Absent relatives: 0 throughout.
  - Cause: the normalizer w (sample-wide mean hits per classified read) is dragged down by host-artefact reads,
    and a true read's hits split between its species and conserved ancestor-level k-mers.
  - Fix to try: a per-clade expected hits/read from the database (the fraction of a genome's k-mers that fall
    within its clade), recorded at build time.

### R1-R2 (2026-10-20)
- **R1** (R0 profiles, no new data):
  - Protocol bias is multiplicative and consistent (within/across SD 0.15). Cross-validated calibration cuts the
    ONT-vs-Illumina Bray-Curtis from 0.42 to 0.07.
  - Host-background subtraction restores composition at 1:10^2 (0.42 -> 0.11) and 1:10^3 (0.76 -> 0.21).
  - A decontam-style frequency test separates mock members from contaminants and host artefacts (AUROC 0.98).
- **Report columns** (0c92af9): per-clade lambda, present k-mers and containment. True species sit at
  containment 2e-3 to 1.5e-2; host artefacts at ~1e-6 (thousands of reads on ~10 k-mers). The host-spiked
  libraries are heavily duplicated (lambda 3-19 on few distinct k-mers), so lambda carries duplication as well as
  coverage.
- **R2: gather paper samples (podar, Zymo D6300), database = gather's matched genomes (circular for presence).**
  - Agreement with the paper's read mapping: bonsai per-species reads give Spearman 0.80 (podar) and 0.97 (Zymo).
    Genome level is 0.79 / 0.34 against all mapping, and 0.73 / 0.26 against unique mapping.
  - 91% (podar) and 70% (Zymo) of reads resolve to a genome leaf.
  - Gather false positives (genomes mapping covers < 50%): bonsai coverage AUROC 0.96 / 0.94, gather 0.97 / 0.88.
  - **NOT a fair head-to-head with gather.** Its f_unique_weighted is an abundance-weighted fraction of the
    sample's known k-mers, not a read fraction. It gives Cryptococcus 23% of Zymo D6300, where the paper's own
    mapping and bonsai give < 1%.
  - **Fair version (to do):** run sourmash gather ourselves on the same 300 MB heads and genomes, and score both
    against mapping.
  - hu-s1 and iHMP are held: the podar build peaked at 7.18 GB, and hu-s1 projects to ~10 GB.

## Extensions beyond pairwise (PI 2026-10-20): design and pre-specified endpoints
Every node carries two dimensions, count (reads, hits) and cardinality (distinct k-mers), joined through lambda
(coverage, containment). Every extension below therefore yields one statistic per node, per dimension.

### X1. Query vs panel of normals (PoN)
- **Panel:** per node and dimension, the sufficient moments of a beta-binomial over N normal profiles
  (rate mu_v and node-specific overdispersion rho_v; containment mean and dispersion). It is mergeable like
  profiles: adding normals adds moments, and leave-one-out is a subtraction.
- **Query:** the predictive tail p at the query's own depth replaces taxdiff's single global rho.
- **Endpoints:**
  - (a) Leave-one-normal-out calibration: p uniform (KS) and FPR <= alpha + 2 SE, per dimension.
  - (b) Spiked deviations (merge a synthetic profile into a held-out normal) are detected, and localized to the
    spiked node, per dimension.
  - (c) Power against panel size N.
- **Kill criterion:** if (a) fails for every reasonable N, the beta-binomial does not describe between-sample
  variation and needs replacing (e.g. a Dirichlet-tree).

### X2. e-values
- **Testing:** per node, e = predictive likelihood ratio of query vs panel (alternative: a mixture over effect
  sizes, so e is a valid betting score). Nested clades make per-node tests dependent, so FDR uses e-BH (Wang &
  Ramdas), which is valid under arbitrary dependence, unlike BH. Merging is valid by averaging across nodes or sets,
  and by product across independent samples, e.g. repeated samples of one subject.
- **Search:** E-value = the expected number of database samples scoring >= s under the null, from a tail fit to
  null scores across the database.
- **Endpoints:** e-BH FDR on the X1 spike-ins (realized FDR <= alpha); power relative to BH on the same p-values;
  and E-value calibration on held-out queries (count of null hits with E <= x is about x).

### X3. Interpretation
- **Localization:** the minimal set of deepest nodes whose conditional tests explain the difference, each with
  effect sizes and CIs (as in the ONT-vs-Illumina call within Enterobacteriaceae).
- **Set enrichment (GSEA-like):** nodes ranked by the signed statistic. Sets come from traits (BugSigDB
  signatures, butyrate producers, oxygen tolerance, Gram stain...).
  - The null must be tree-aware: traits are phylogenetically clustered, so naive label permutation is
    anti-conservative. Use within-clade permutation or tree-rotation nulls.
  - **Endpoint:** FPR of the naive vs the tree-aware null on random, phylogenetically clustered sets (the naive
    null should fail; the tree-aware one should not).

### X4. Function: a second tree
- The bonsai "taxonomy" can be any tree. A functional tree (pathway -> gene family -> gene, with k-mers as leaves)
  reuses profiles, count + cardinality, taxdiff, the PoN and e-values unchanged.
- Cardinality then gives direct evidence of gene or pathway content, not function inferred from which taxa are
  present.
- Links to metabolism go through annotations (KEGG / GMM modules, AGORA2 reconstructions).
- Needs: genome annotation of the database (a download / compute decision for the PI).
- **Endpoint (mock):** pathway presence and absence in D6331 and podar, from the functional tree, against the
  reference genomes' annotations.

Priority: X1 + X2 (testing), then X3 enrichment, then X4. X1 code (panel file and fitting) can be written and
unit-tested now. Real evaluation needs a normals corpus (to be decided with the PI).
