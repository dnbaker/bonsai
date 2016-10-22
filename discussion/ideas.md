Issues:
  1. Minimizer
    1. Frequency
      1. Total number of observations
      2. Total number of features sharing observation
    2. Lexicographic
      1. Correct for frequent AAA and ACA lmers. (lmer referring to the minimizer sequence) [Discussed in some of our papers]
  2. Windows
    1. Sliding [like most do]
      1. With sufficiently-sized windows for long read mapping, 
    2. Fixed windows in sequence.
      1. This reduces the number of potential locations but reduces advantages due to minimizer presence.
      2. Particularly useful for long sequence to long sequence alignment, such as [minimap](https://github.com/lh3/minimap).
  3. Redundancy
    1. Difficult to test between kmers of complex spaced seeds, easier to do with unspaced seeds.
    2. Systematically skip positions. [consider every i positions]
      1. Complex to coordinate with spaced seeds.
    3. Randomly select positions to use/skip.
  4. Seed generation
    1. Rasbhari
    2. Simply select ones which have been generated already.
  5. Paired-end seed use. (Large gaps in spaced seed to simulate kmer co-occurrence.)

Applications:
  1. Single genome
    1. Transcript assignment and subsequent de novo assembly using [fermi-lite](https://github.com/lh3/fermi-lite). (RNA/DNA)
    2. Bin genome by exon and intermediate regions. Treat coding exons as separate due to their improved uniqueness.
  2. Many genomes
    1. Use taxonomy to constrain

TODO:
  1. Templated kmer encoding with minimizers, window-size, minimizer criterion.
    2. File I/O (easiest would be using kseq to "eat" both fastq and fasta.
    3. To multithread, doing multiple smaller hash tables and performing unions would be a hacky way to get there.
    4. For our hash map, I think I want to use libcuckoo. For lookups, there are up to two cache misses,
       and it is easy to get over 90% occupancy with high performance, and it scales better than even
       the Intel concurrent hashmap.
  2. Establish classification strategy
    1. Ambiguity resolution ?
    2. Use NCBI taxonomy?
      1. In addition to Kraken, Centrifuge (just published) uses it. I think this strategy is entirely appropriate.

Analysis:
  1. 
  2. 
