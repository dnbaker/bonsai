# Changelog

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
