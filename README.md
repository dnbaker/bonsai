Bonsai: Flexible Taxonomic Analysis and Extension
===============

Bonsai contains varied utilities for taxonomic analysis and classification using exact subsequence matches. These include:
* A high-performance, generic taxonomic classifier
  * Efficient classification
    * 20x as fast, single-threaded, as Kraken in our benchmarks, while demonstrating significantly better threadscaling.
  * Arbitrary, user-defined spaced-seed encoding.
    * *Reference compression* by windowing/minimization schemes.
    * *Generic minimization* including by taxonomic depth, lexicographic value, subsequence specificity, or Shannon entropy.
  * Parallelized pairwise Jaccard Distance estimation using HyperLogLog sketches and is dramatically more accurate than comparable tools while also significantly outperforming them in speed.
* An unsupervised method for taxonomic structure discovery and correction.
* A threadsafe, SIMD-accelerated HyperLogLog implementation.
* Scripts for downloading reference genomes from new (post-2014) and old RefSeq.

Tools can be built to work with zstd instead of gzip by being built with a '_z' suffix. (e.g., bonsai_z).

All of these tools are experimental. Use at your own risk.


Build Instructions
=================

`make bonsai`

Unit Tests
=================
We use the Catch testing framework. You can build and run the tests by:

`cd bonsai && make unit && ./unit`

