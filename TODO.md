Updated 1/3/18
For a few points to remember:

1. Jaccard index is itself just an approximation of ANI, and a multiset Jaccard index would be closer to ANI.
2. Metannot/wavelet trie compression and "memoized" compression (Florian project) as points of comparison.
3. The Counting Quotient Filter is space efficient and the RSQF is quite fast. Is there any way to modify it to perform either multiset cardinality estimation (e.g. for multiset Jaccard)
   or an improvement over standard bloom filters for search?

TODO:

1. Test where nodes are being added by greedy methods. Be willing to consider adding a surplus of nodes and trim them down later.
2. Modify tree construction to give a given tax id to the first genome that contains it, but be willing to modify it and move it to be a child
if more genomes labeled with that taxid occur.
  1. This function should take as an argument a set of genomes, annotations for taxids for each, and the NCBI taxonomy, and then construct as concise
     a taxonomy as possible, providing each end genome with a node, even if there have to be a significant number of children
     for the 'stem' nodes.
3. Come up with a metric for how much is gained/lost vs Kraken and vs full CDBG.
4. Is HLL or bloom filter sketching an option for getting at kmer set overlaps with smaller memory requirements?
