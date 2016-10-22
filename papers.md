# Papers

## Still in draft form


### Spacing seeds -- selecting several spaced seeds to use
[Multiple spaced seeds for homology search](http://bioinformatics.oxfordjournals.org/content/23/22/2969.full.pdf)
[rasbhari](https://arxiv.org/pdf/1511.04001v2.pdf)

This results in a separate table and separate encodings of the query sequence for each spaced seed. I would likely consider using
the libcuckoo hash table, which incorporates the advances in
[this paper](http://www.cs.princeton.edu/~xl/Publications_files/cuckoo-eurosys14.pdf).

### Storing signatures rather than kmers themselves
[KMC 2: Fast and resource-frugal k-mer counting.](https://arxiv.org/abs/1407.1507) # Stores only signatures, not kmers, during counting.
  This results in imperfect collision resolution, but only a few (5-7) bits need to be stored per key. Unfortunately, you can not get the kmer itself out of the database.
  This makes this a sketch-like data structure.
  Because building a concurrent, performant hash table is nontrivial, I think I want to prototype using existing implementations.

### Constraining value space with taxonomy
[Centrifuge](http://genome.cshlp.org/content/early/2016/10/17/gr.210641.116.abstract)


### How many seeds to use
[Multiple seeds sensitivity using a single seed with threshold.](https://www.ncbi.nlm.nih.gov/pubmed/25747382)





### Practical Examples

|Title | Tool | Authors | Link|
|------|------|--------|------|
|Kraken| Kraken| Wood, Salzberg| https://ccb.jhu.edu/software/kraken/|
|Kallisto| | https://arxiv.org/pdf/1510.07371|

