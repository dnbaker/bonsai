# Literature review for bonsai profile applications (2026-10-20)

This review is scoped to what bears on the profile applications and their results so far. Each entry gives what
exists, what it implies for us, and what to do. It is not exhaustive. It comes from one round of web searches plus
prior knowledge; items marked *(memory)* were not re-checked in this pass.

## 1. Coverage correction of k-mer containment: the M3 / M6 / R0 low-coverage problem
- **sylph** (Shaw & Yu; bioRxiv 2023, Nat Biotech 2024 per press) models the coverage of a genome's k-mers in a sample
  as a zero-inflated Poisson. Divergent k-mers make the zero inflation. It infers an effective coverage, then a
  coverage-adjusted containment ANI. It was the most accurate of 7 profilers on CAMI2 Marine, and its ANI gave a
  signal orthogonal to abundance in a Parkinson's MWAS.
  [PMC](https://www.ncbi.nlm.nih.gov/pmc/articles/PMC12339375/), [bioRxiv](https://www.biorxiv.org/content/10.1101/2023.11.20.567879.full.pdf)
  - **For us:** this is the principled version of "counts + sets together". Our per-clade counters give
    multiplicity = hits / distinct. Under a zero-truncated Poisson that solves for an effective coverage λ.
    - The expected distinct count for a fully present clade is g (1 - e^-λ), where g is the clade's database k-mers.
    - u / that is containment π, and π^(1/k) is an ANI.
  - **It replaces two things that failed:**
    - **R0's sample-wide w normalizer:** λ is per clade, so host reads do not move it. False-positive taxa pile hits
      on a few k-mers, so they get high λ and tiny π.
    - **M6's internal-node bias:** estimate per taxon and sum up the tree.
  - **Caution:** at very low coverage (multiplicity ≈ 1), λ is unidentifiable from (hits, distinct) alone. sylph
    uses per-k-mer counts. A small per-clade multiplicity histogram (counts of k-mers seen 1, 2, 3+ times) would be
    a cheap profile extension that pins λ down.
  - **For the paper:** sylph is a required comparison for species detection and ANI. Our claimed difference is
    mergeable, taxonomy-resolved counters, and tests over them.
- **sourmash** FracMinHash containment and containment ANI (Irber et al. 2022; Hera et al. *(memory)*) estimate the
  same quantities on whole sketches.

## 2. False-positive taxa and distinct k-mers: M2, R0
- **KrakenUniq** (Breitwieser, Baker & Salzberg, Genome Biol 2018) counts the unique k-mers each taxon's reads cover,
  with HyperLogLog (16 KB, ~1% error). It catches false positives especially at low abundance.
  [Genome Biol](https://genomebiology.biomedcentral.com/track/pdf/10.1186/s13059-018-1568-0)
  - **For us:** our `-U` report is this. Profiles add mergeability, MinHash for set comparisons, and cross-sample
    tests. Claim novelty only for those.
- **Host reads classified as microbes:** Gihawi, Salzberg et al. (mBio 2023) showed that much of the "cancer
  microbiome" bacterial signal was human reads, misclassified through database contamination. A data transformation
  then created artificial per-tumor signatures.
  [PMC](https://www.ncbi.nlm.nih.gov/pmc/articles/PMC10418105/)
  - **For us:** R0's human-only controls classify 1.3e-3 of reads, across 6-9 species, which is exactly this failure
    mode. Any host-rich application (tissue, low biomass) needs three things:
    - (a) a host-aware database (human/T2T in it, or a host filter);
    - (b) the coverage or ZIP containment filter, which removed every host artefact in R0;
    - (c) reporting of the host background.
  - **Paper endpoint:** R0's human-only false-positive count, before and after each safeguard.

## 3. Contamination and negative controls: the kitome in R0's blanks
- **decontam** (Davis et al., Microbiome 2018) uses two signals. Frequency: contaminant frequency rises as sample
  DNA concentration falls. Prevalence: contaminants are more prevalent in negative controls. It needs several
  controls and does not handle cross-contamination.
  [Microbiome](https://link.springer.com/doi/10.1186/S40168-018-0605-2)
- **SCRuB** (Nat Biotechnol 2023) models each sample as a mixture of contamination and true sources, and models
  well-to-well leakage using plate positions.
  [PubMed](https://pubmed.ncbi.nlm.nih.gov/36928429/)
  - **For us:** R0's dilution series is a direct test of decontam's frequency principle. True D6331 taxa track
    microbial DNA input (slope -0.98); kit contaminants should stay flat or rise in relative terms as input falls.
  - **taxdiff fits here:** "taxa over-represented in a sample relative to its extraction blanks, beyond
    overdispersion" is a principled contamination screen at every node.
  - **To do:** R0b, kitome detection on the dilution series and blanks, scored against D6331 truth.

## 4. Protocol bias: R0's ONT vs Illumina gap
- McLaren, Willis & Callahan (eLife 2019): metagenomic measurement bias is multiplicative. True relative abundances
  are multiplied by taxon- and protocol-specific efficiencies that are constant across samples, so the bias is
  correctable with mock controls.
  [eLife](https://elifesciences.org/articles/46923)
  - **For us:**
    - R0 shows species Bray-Curtis of 0.43 between ONT and Illumina, vs 0.07-0.14 within a platform.
    - Under the multiplicative model, ratios between taxa within a protocol cancel the bias. taxdiff's
      sample-vs-sample log rate ratios on one protocol are unaffected, but cross-protocol comparisons need
      calibration.
    - **Test now, no new data:** estimate per-taxon efficiencies per platform from the D6331 runs, check they are
      constant across replicates and dilutions (the model's prediction), then apply them and re-measure the
      cross-platform Bray-Curtis.
    - **Profile op:** `profile calibrate` (efficiencies from a mock profile, applied to other profiles).

## 5. Differential abundance: `profile diff`, M4
- **Nearing et al.** (Nat Commun 2022): 14 DA methods on 38 datasets give strikingly different results. Counts track
  sample size, depth and effect size. ALDEx2 and ANCOM-II are the most consistent, and the authors recommend
  consensus. A 2025 F1000 follow-up re-tests this with synthetic data (metaSPARSim, sparseDOSSA2).
  [DOAJ](https://doaj.org/article/10cb599d20bf4079bcf69b1849ac9695)
- **corncob** (Martin, Witten & Willis, AoAS 2020): beta-binomial regression with covariate-dependent overdispersion.
  It tests differential abundance and, separately, differential variability (dysbiosis); the two find different
  taxa. [arXiv](https://arxiv.org/pdf/1902.02776)
- **Compositional methods** *(memory)*: ANCOM-BC/BC2, LinDA, MaAsLin2/3 and ALDEx2's CLR.
- **For us:**
  - **Benchmark:** `diff` reads/k-mers on real cohorts against ALDEx2, ANCOM-BC2, corncob and MaAsLin, in Nearing's
    design: concordance, plus FPR by within-group label shuffles. Synthetic power uses sparseDOSSA2-like simulators,
    not only our Dirichlet model.
  - **Differential variability:** a natural addition on counters too. Is the within-group k-mer dispersion
    different?
  - **taxdiff's beta-binomial ρ:** this is corncob's model for two samples. For groups, fit ρ per taxon and test
    with covariates, which is corncob-like.

## 6. Tree-structured testing: taxdiff
- **Dirichlet-tree multinomial (DTM):** a cascade of local Dirichlet-multinomials at the internal nodes of the tree.
  **PhyloScan** (Tang, Ma & Nicolae, AoAS 2018) scans the tree so nodes borrow strength from parents and children.
  **phyloMDA** (R) implements DTM and tree-guided estimation. **PhAAT** tests branches and refines them.
  A paired-multinomial taxonomic-tree test identifies differing subtrees.
  [phyloMDA](https://www.ncbi.nlm.nih.gov/pmc/articles/PMC9169257/), [paired tree test](https://ar5iv.labs.arxiv.org/html/1702.04808)
  - **For us:** taxdiff's test conditional on the branching ancestor is a local DTM binomial, so cite DTM and
    PhyloScan as the frame.
  - **Gaps to close:**
    - **Hierarchical FDR over the tree** (Yekutieli-style top-down testing *(memory)*), instead of flat BH over nodes.
    - **Siblings:** a sequential sibling decomposition (child i against the remaining siblings) makes the node G
      statistics additive. That gives an exact LLR budget per subtree.
    - **Effect sizes with CIs** at every node (done).

## 7. Strain sharing and transmission: M3
- **Valles-Colomer et al.** (Nature 2023): more than 9,700 metagenomes and over 10M strain-sharing events.
  Mother-to-infant gut transmission is substantial; cohabiting adults share a median 12% (gut) and 32% (oral) of
  strains. [PMC commentary](https://www.ncbi.nlm.nih.gov/pmc/articles/PMC10989920/)
- **inStrain** popANI vs consensus ANI: in the inStrain benchmarks, StrainPhlAn's minimum reported threshold is
  99.97%, while inStrain's popANI reaches 99.99996%. [inStrain benchmarks](https://instrain.readthedocs.io/en/latest/benchmarks.html)
  - **For us:** 31-mer containment cannot reach SNP-level popANI resolution. So "same strain" claims from counters
    should be framed as "candidate sharing at containment ANI ≥ x" (sylph-style). Confirm a subset with
    inStrain/StrainPhlAn; do not claim transmission.
  - **M3's measured limit:** k-mer sharing tracks strains absent from the database only at ≥ 3x. ZIP coverage
    correction (§1) is the principled attempt to push that down.

## 8. Corpus-scale baselines and search: the PI's "large collections of corpora"
- **sourmash branchwater** (Irber, Pierce-Ward & Brown; bioRxiv 2022, JOSS): FracMinHash containment search over
  ~1M SRA metagenomes. It recommends queries ≥ 50 kb.
  [bioRxiv](https://www.biorxiv.org/content/10.1101/2022.11.02.514947.full.pdf), [JGI service](https://branchwater.jgi.doe.gov/about)
- **Logan** (2024): assemblies (unitigs and contigs) of 26.8M SRA accessions (48 Pbp, 88% of the SRA by bases),
  compressed more than 100x, on a public AWS bucket. Logan-Search queries the whole SRA in under 5 min.
  [PMC](https://www.ncbi.nlm.nih.gov/pmc/articles/PMC12424806/), [GitHub](https://github.com/IndexThePlanet/Logan)
- **curatedMetagenomicData** *(memory)*: ~20k uniformly processed human metagenomes with curated metadata. This is
  the natural baseline corpus.
- **For us:**
  - **Prior art:** corpus search exists (branchwater, Logan-Search). Our contribution is taxonomy-resolved,
    mergeable profiles with tests, so baselines become statistical objects:
    - a per-node reference distribution (mean, overdispersion ρ, prevalence) built by merging and summarizing a
      corpus;
    - taxdiff of a sample against its baseline;
    - "how much does X matter" as variance partitioning of per-node rates across corpus metadata (study, platform,
      body site, disease).
  - **Streaming cost:** ENA runs ~0.9 MB/s here, so 1000 samples x 60 MB heads ≈ 18 h.
  - **Logan unitigs are ~100x smaller.** But profiles need read multiplicities (counts), and unitigs keep only k-mer
    abundance annotations *(to check: Logan/BCALM headers carry mean k-mer abundance)*. An abundance-weighted
    unitig classifier could approximate read counts. Worth a feasibility test before any large stream.
  - **Baseline needs:** a broad gut database (GTDB or UHGG representatives). That is a new download, to be sized and
    asked.

## 9. Chaining-based classification (PI's tangent)
- No 2024-25 metagenomic classifier built on colinear chaining turned up in this search. Nearby work:
  - MetaTrinity (arXiv 2023) skips chaining for speed, using seed-hit sorting and edit-distance approximation.
  - mapquik (Ekim, Sahlin, Berger *(memory for authors)*) uses linear-time pseudochaining on k-min-mers.
  - Das & Schatz (2022) use sampling and sketching for long-read classification.
  - The exact-chaining theory of Rizzo, Cáceres & Mäkinen (2025) shows the ChainX heuristic can be suboptimal.
  - MetaMaps (Dilthey et al. 2019) uses approximate mapping plus EM *(memory)*.
  [MetaTrinity](https://arxiv.org/pdf/2311.02029), [chaining theory](https://arxiv.org/abs/2506.11750v2)
  - **For us:** the gap looks real, though this is not exhaustive. The cheap test is chain bonsai's own hits by
    reference position per clade, then re-run M2's staggered+novel condition and R0's host controls. Scattered
    conserved-k-mer hits from absent organisms or host reads should fail to chain. The Centrifuger backbone (r-index)
    is the scalable version.

## What this changes in the plan (priority order)
1. **ZIP effective-coverage containment per clade (§1).** It replaces the read-saturation normalizer and fixes M6's
   per-node estimate. If possible, extend profiles with a small per-clade multiplicity histogram. Re-run M2, M3, M6
   and R0. Compare against sylph on R0.
2. **Protocol-bias calibration (§4)** on R0's D6331 platform data. Test the multiplicative model, then add
   `profile calibrate`.
3. **Contamination screen (§3)** on R0's dilution series and blanks: taxdiff against blanks, and a decontam-style
   frequency test.
4. **Tree FDR and sibling decomposition** for taxdiff (§6). Benchmark `diff` against ALDEx2, ANCOM-BC2 and corncob on
   a real cohort (§5).
5. **Corpus baselines (§8).** Size the gut database and a baseline corpus (curatedMetagenomicData accessions, or a
   Logan-unitig feasibility test). Needs a PI OK on volume.
6. **Chaining prototype (§9)** on bonsai hits.
