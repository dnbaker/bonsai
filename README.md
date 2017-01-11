[EM]: Exact-Match Playground
===============

Contains utilities for exploring the tradeoff space between various ways of encoding exact sequence databases.

Don't ask, don't touch, don't use, Don't Panic.


Build Instructions
=================

`make`

Unit Tests
=================
We use the Catch testing framework. You can build and run the tests by:

`make unit_tests`


Prior Work
================

We modified [khash](https://github.com/attractivechaos/klib) for our primary database,
borrowed and modified code from [Kraken](https://github.com/DerrickWood/kraken) for XOR masking,
kmer encoding, taxonomy tree building and querying, and classification.

We used a [custom allocator from StackOverflow](http://stackoverflow.com/questions/12942548/making-stdvector-allocate-aligned-memory)
for allocating aligned memory for our HyperLogLog implementation.

