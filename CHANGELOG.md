# Changelog

## Unreleased

### Added: sample profiles and cross-sample analysis
- `classify -P <out.bprof>`: the sample's classification result as sufficient statistics
  (per taxon: own reads and k-mer hits, parent, rank, clade database k-mers, and the raw
  distinct-k-mer counter), written before clade folding. Self-contained (no taxonomy needed
  to read it) and mergeable: merging profiles of read subsets gives byte for byte the
  profile of all the reads. `-m` counts with b-bit MinHash (tight shared-k-mer estimates).
- `bonsai profile report`: the `classify -U` report (identical) and `-B` sketches from a profile.
- `bonsai profile merge`: lanes/replicates/pools.
- `bonsai profile compare`: pairwise Bray-Curtis, presence Jaccard and taxonomy earth mover's
  distance (unit edges); `-K` per-clade distinct-k-mer Jaccard and containment (same vs
  different strains). `bonsai profile table`: taxon x sample matrix (reads or distinct k-mers).
- `bonsai profile index` / `search`: nearest samples to a query profile by tree EMD,
  Bray-Curtis, k-mer Jaccard or containment.
- `bonsai profile diff`: differential abundance between two groups: per clade, Welch t on
  log proportions with permutation p and BH q; `-K` PERMANOVA on per-clade k-mer Jaccard
  distances (k-mer content shifts at equal abundance) and B's pooled k-mers unseen in A;
  community PERMANOVA on tree EMD and Bray-Curtis. Seeded; exact enumeration when small.
  Counters hold database k-mers only, so k-mer content is measured over the database.
- `make check`: unit tests and the CLI smokes (`test/cli_profile_smoke.sh` is new).

### Changed: compact k-mer -> taxon databases by default
- `bonsai build` (lex LCA, the default mode) writes the compact `cht_map_t` database, as `-Z`
  did; `-L` writes the legacy khash database. Non-lex modes still build khash (logged). A gzip
  output path or `-z` needs `-L` (compact databases are uncompressed); `-z` is now accepted
  (it was documented but missing from the option string). classify reads both formats.

### Changed: classify's parent map on the sketch cht table
- classify reads the taxonomy into `CompactTaxMap` (sketch cht; same parsing as
  build_parent_map) instead of khash. The classifier, confidence walk-up, identity
  placement, lineage evidence and tree resolution are generic over the parent map, and
  regress_classify_modes checks khash and CompactTaxMap parents give identical output.
- resolve_tree stops at a taxon missing from the taxonomy (it read past kh_end before).

### Changed (results change, -U)
- A distinct-k-mer counter is exact iff it holds at most 512 distinct k-mers, else a sketch,
  whatever order hits arrived in (before, a taxon with 513 to ~1023 distinct k-mers could be
  reported exactly or by its sketch depending on hit order). The -U report now depends only on
  the k-mer sets.
- The -U report and profiles look parents up in the sketch cht table (`CompactTaxMap`).

## v0.6.0 (2026-10-20)

### Fixed
- Chunked FASTQ reading: a 1-byte heap overflow when headers carry a comment (e.g. SRA HiFi
  `@SRR... 1/1`), and uninitialized record-slot reuse when a chunk had more reads than any
  earlier one (the record array also could not grow past 4096). Long-read classify output was
  nondeterministic (on Zymo D6331 HiFi ~1 run in 4 aborted; survivors had a few corrupted or
  empty-named records per 23k reads). Now ASan-clean and byte-identical across runs/threads.
- Errors exit 1 instead of aborting.

### Added (all opt-in; default classify output unchanged)
- Compact LCA databases on sketch `cht_map_t`: `bonsai build -Z`, `bonsai compact`.
- Exact 64/128-bit k-mer keys for k <= 64; spaced seeds and minimizer windows for 32 < k <= 64.
- Six-frame translated (protein) k-mer databases. Builds with C++20.
- classify: Kraken2-style confidence walk-up and minimum hit groups; lineage evidence (`-E`);
  hits-only (`-R`) and quality-aware (`-Q`) confidence; identity-capped placement (`-I`);
  KrakenUniq-style distinct-k-mer report (`-U`, shared per-taxon counters, identical at any
  thread count); per-clade b-bit MinHash sketches (`-U ... -B dir`).

### Changed
- `ankerl::unordered_dense` replaces `ska::flat_hash_map`; hll -> sketch v0.20.0
  (concurrent sketch adds call `add_threadsafe`).
