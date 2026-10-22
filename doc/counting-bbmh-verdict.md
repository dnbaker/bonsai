# Counting-minhash verdict for UniqCounter (PI directive, 2026-10-21)

Question: does the sketch submodule offer a *counting* b-bit MinHash that
could back per-taxon counters in `include/bonsai/uniq.h`, and would it let us
drop the exact SharedTally phase (or the exact/atomic counters)?

## What the submodule (hll @ 13e0c0d, sketch 0.20.1) actually has

No counting variant of `BBitMinHasher` exists. `BBitMinHasher<T>` (sketch/bbmh.h:912)
registers hold only the per-bucket minimum hash; `add_threadsafe` is an
atomic-min CAS (bbmh.h:987, already load-gated). Nothing counts multiplicities.

Counting structures that do exist, and why none is a drop-in:

| structure | header | what it counts | merge semantics | lock-free adds |
|---|---|---|---|---|
| `bcqf_t<KEY_BITS,VAL_BITS>` | sketch/cqf.h:1437 | exact per-fingerprint counts to 2^VAL_BITS-1 | count SUM (fp collisions merge counts, documented overcount) | no: slot shifting under a table spin/mutex |
| `bcf_t<FP_BITS,VAL_BITS>` | sketch/cqf.h:2510 | fixed-capacity fingerprint counting filter | sum | no |
| `cbfbase_t` | sketch/cbf.h:33 | counting Bloom (saturating, probabilistic) | sum of noisy counts | atomic increments possible |
| `ccmbase_t` | sketch/ccm.h:189 | count-min sketch (heavy hitters) | sum | yes |
| `bmh_t` / BagMinHash | sketch/bmh.h:219 | multiset cardinality / weighted Jaccard | n/a (finalized estimator) | no |

## Verdict

**Keep bbmh + the exact phase; do not adopt a counting sketch.** Three reasons:

1. **Merge semantics.** The profile contract is set union with byte-exact
   merge: counter merge is register-max, state depends only on the hash set
   (uniq.h canonicalize). Every counting structure in the submodule merges by
   count *addition* — that is multiplicity semantics, not set-union semantics.
   "Max per fingerprint" (the union-preserving way to merge counts) is not
   offered by any of them.
2. **Locks.** Quotient-filter inserts shift slots and take a table lock; the
   PI's directive #2 was to get *out* of the lock business. bbmh post-promotion
   adds are already lock-free CAS.
3. **Nothing to remove.** The atomic SharedTally reads/kmers counters are not
   per-k-mer structures — two words per taxon, batch fetch_add — and the
   report's duplication metric (hits / distinct) comes from those totals, not
   from per-k-mer counts. A counting sketch would not replace them.

The natural "counting bbmh" — register = (min, multiplicity-of-the-min) — is
order-independent and union-mergeable (max of mins; sum multiplicities on
ties), so it *would* preserve byte-exactness. But the multiplicity it tracks
is that of a single hash value per register, i.e. a weighted-minhash quantity;
it is not per-taxon duplication and removes no existing counter. Not worth the
register width.

Follow-up if per-k-mer multiplicity ever matters (e.g. abundance-aware
containment): `bcqf_t<24,8>` per hot taxon, merged by explicit
max-with-collision-overcount, serialized separately — a new counter *kind*,
not a replacement.
