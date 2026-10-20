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
