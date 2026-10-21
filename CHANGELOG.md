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

### Fixed: spaced seeds were single-stranded
- Spaced seeds were never canonicalized (the encoder switched canonicalization off for a spaced
  spacer), so a database built from forward-strand references missed reverse-strand reads: on the
  D6331 mock, forward-only reads classified 0.999 and reverse-only 0.010-0.025, i.e. about half of
  all reads lost.
- Spaced keys are now canonical under a both-strand contract: min(pattern on the window, pattern
  on its reverse complement), valid for any pattern (asymmetric included).
- Contiguous seeds are unchanged (classification and -U byte-identical).
- **Spaced databases built before this change must be rebuilt.**

### Fixed: feature-count (-f) and taxonomic-depth (-t) minimized builds
These builds could not have worked; each fix below was needed for `prebuild -f|-t`
then `build -f|-t` to give a usable database.
- `build -f|-t` never reached the minimized build: the mode test `LEX == mode || ENTROPY` was
  always true, so it built an entropy LCA map, wrote it over the phase-1 input (the first
  positional) and read the intended output path as a genome.
- `prebuild` decremented the stored spacing a second time (Database already stores it as
  differences), writing -1s; build then encoded different k-mers ("Missing kmer from database").
- `kh_set` (khash64.h) dropped the slot kh_put returned for a new key: the value was written
  past the end of the values array and the slot kept uninitialized memory (the minimized
  database's taxa were garbage). Its CAS retry also never reread the old value.
- `hash_score` (the phase-1 minimizer score) had its lookup test inverted, reading vals[kh_end].
- `update_feature_counter` stored a depth as the taxon, never incremented the genome count, and
  took LCAs of packed values (now as update_frequency_counter does).
- Minimizer ties under the phase-1 score are broken by a k-mer hash (the low bits were the
  taxid, so every window kept a different k-mer: no minimization at w > k).
- cli_classify_smoke: -f and -t end to end (phase-1 input untouched, reads classified), and
  with w > k: density as lex, and k-mers shared between genomes avoided (lex keeps ~450 of a
  test set, -f/-t keep none).

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
