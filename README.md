[EMP]: Exact-Match Playground
===============

EMP contains varied utilities for taxonomic analysis and classification using exact subsequence matches. These include:
* A high-performance, generic taxonomic classifier
    0. Efficient classification
        1. *20x as fast as Kraken* in our benchmarks.
    1. Arbitrary, user-defined spaced-seed encoding.
    2. *Reference compression* by windowing/minimization schemes.
    3. *Generic minimization* including by taxonomic depth, lexicographic value, or subsequence specificity.
* An unsupervised method for taxonomic structure discovery and correction.
* A threadsafe, SIMD-accelerated HyperLogLog implementation.
* Scripts for downloading reference genomes from new (post-2014) and old RefSeq.

All of these tools are experimental. Use at your own risk.


Build Instructions
=================

`make`

Unit Tests
=================
We use the Catch testing framework. You can build and run the tests by:

`make unit`


Prior Work
================

We modified [khash](https://github.com/attractivechaos/klib) for our primary database,
borrowed and modified code from [Kraken](https://github.com/DerrickWood/kraken) for
kmer encoding, taxonomy tree building and querying, and classification.
