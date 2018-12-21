Bonsai: Flexible Taxonomic Analysis and Extension [![Build Status](https://travis-ci.com/dnbaker/bonsai.svg?branch=master)](https://travis-ci.com/dnbaker/bonsai)
===============

Bonsai contains varied utilities for taxonomic analysis and classification using exact subsequence matches. These include:
* A high-performance, generic taxonomic classifier
  * Efficient classification
    * 20x as fast, single-threaded, as Kraken in our benchmarks, while demonstrating significantly better threadscaling.
  * Arbitrary, user-defined spaced-seed encoding.
    * *Reference compression* by windowing/minimization schemes.
    * *Generic minimization* including by taxonomic depth, lexicographic value, subsequence specificity, or Shannon entropy.
  * Parallelized pairwise Jaccard Distance estimation using HyperLogLog sketches, which has recently migrated to [dashing](https://github.com/dnbaker/dashing).
* An unsupervised method for taxonomic structure discovery and correction. (metatree)
* A threadsafe, SIMD-accelerated HyperLogLog implementation, which has migrated to [hll](https://github.com/dnbaker/hll).
* Scripts for downloading reference genomes from new (post-2014) and old RefSeq.

Tools have been compiled using both zlib and zstd, which means that they can transparently consume zlib-, zstd-, and uncompressed files.

All of these tools are experimental. Use at your own risk.


Build Instructions
=================

`cd bonsai && make bonsai`

Unit Tests
=================
We use the Catch testing framework. You can build and run the tests by:

`cd bonsai && make unit && ./unit`


Usage
================

Usage instructions are available in each executable by executing it with no options or providing the `-h` flag.


For classification purposes, the commands involved are `bonsai prebuild`, `bonsai build`, and `bonsai classify`.
prebuild is only required for taxonomic or feature minimization strategies, for which case database building requires double the memory requirements.
Unless you're very sure you know what you're doing, we recommend simply `bonsai build` with either Entropy or Lexicographic minimization.

To build a database with k = 31, window size = 50, minimized by entropy, from a taxonomy in `ref/nodes.dmp` and a nameidmap in `ref/nameidmap.txt` and store it in in `bns.db`
```
bonsai build -e -w50 -k31 -p20 -T ref/nodes.dmp -M ref/nameidmap.txt bns.db `find ref/ -name '*.fna.gz'`
```

To prepare the above, the script in `python/download_genomes.py` can be used. The default of downloading all available genomes can be run by `python python/download_genomes.py --threads 20 all`.
This places downloaded genomes by default into the paths listed above in the `bonsai build` command. These paths can be altered; see `python/download_genomes.py -h/--help` for details.
