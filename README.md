Bonsai: Flexible Taxonomic Analysis and Extension
===============

Bonsai contains varied utilities for taxonomic analysis and classification using exact subsequence matches. These include:
* A high-performance, generic taxonomic classifier
  * Efficient classification
    * 20x as fast as Kraken in our benchmarks.
  * Arbitrary, user-defined spaced-seed encoding.
    * *Reference compression* by windowing/minimization schemes.
    * *Generic minimization* including by taxonomic depth, lexicographic value, or subsequence specificity.
  * Parallelized pairwise Jaccard Distance estimation using HyperLogLog sketches.
* An unsupervised method for taxonomic structure discovery and correction.
* A threadsafe, SIMD-accelerated HyperLogLog implementation.
* A HyperLogLog-based sequence distance estimator which is more accurate and orders of magnitude faster than comparable tools.
* Scripts for downloading reference genomes from new (post-2014) and old RefSeq.

Tools can be built to work with zstd instead of gzip by being built with a '_z' suffix. (e.g., bonsai_z).

All of these tools are experimental. Use at your own risk.


Build Instructions
=================

`make bonsai`

Unit Tests
=================
We use the Catch testing framework. You can build and run the tests by:

`cd make && make unit && ./unit`

