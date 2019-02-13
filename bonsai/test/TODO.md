Tests:
  1. Encoder
    1. Kmers
      1. Unspaced kmers
      2. Particular spaced kmers
      3. Windowing
    2. Minimizing (Test the functions themselves. We'll test the hash tables we build elsewhere)
      1. Lex [Dura lex, sed lex]
      2. Tax
      3. Feature
    3. Counting
      1. Estimate Cardinality
  2. Feature Min
    1. LCA Map
    2. Feature count map
    3. lca2depth
    4. fill_set_*
  3. HLL
    1. Small
    2. Normal value
    3. Addition
  4. NCBI
  5. Spacer
    1. Encoding
    2. Decoding
  6. Util
     1. Read/write khash maps.
     2. LCA
     3. Node depth
     4. Taxonomic Tree Map building
  7. Qmap
    1. Values are always BF until the window is filled.
    2. The element is the lowest element.
