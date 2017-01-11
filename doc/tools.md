### Tools in the code base.

#### Khash
khash is our core hash table, originally written for klib and used as the core hash table
in htslib. Being written in C, one manually has to insert keys and get locations in the table.
Implemented with macros, the functions require a "name" argument to specify which khash table we're
using.

"all" is our 64-bit integer hash set.
"64" is our 64-bit integer hash map with 64-bit integer keys.
"c" is our 64-bit integer hash map with 32-bit integer keys for counting, which "c" stands for.

One gets a value from a khash map using a key being searched for, the name of the hash table,
and a pointer to the hash table.

If it is present in the table, it returns an index into the hash->vals and hash->keys arrays,
and the value is present at kh_val(hash_ptr, index). Otherwise, it returns kh_end(hash).
If it is a set, hash->vals is null.

An example using "64":
```c++
ki = kh_get(64, hash_ptr, key);
if(ki == kh_end(hash_ptr)) {
    fprintf(stderr, "Key is missing from table.\n");
} else {
    fprintf(stderr, "Key is present and available at index %" PRIu64 ".\n", ki);
    fprintf(stderr, "Value is %" PRIu64 ".\n", kh_val(hash_ptr, ki));
}
```

For our multithreading, we use std::async to create futures. When the spawned thread is complete,
(which we test with is_ready(future)), we get the value using future.get().

#### Internal typedefs/structs

**spvec_t** is a vector (dynamically sized array) of spaces between positions in the string which we
are encoding. It is std::vector<uint8_t>.
When it is provided to spacer, we increment all of its values so that we can directly
increment our index variable by the spaces.
(IE, if the argument value is 0 at a position, we
increment our position in the string by 1.)


**Spacer** is our struct containing kmer size,
window size, and our spacing. Because it contains the necessary information for decoding,
we can decode to a string using it. When we encode or decode, we XOR with XOR_MASK to help
distribute the keys in our table and avoid selecting poly-A kmers for our lexicographic minimizers.

**Encoder** is the core encoding structure. It requires a spacer, and it works by
using a qmap_t (see below) to keep track of best-scoring kmers in a window.
Between using different strings,
we use the "assign" method using the string size and string length to tell the Encoder to use that string.
next_minizer() provides the next minimizer. Before doing so, you need to call has_next_kmer() to
find out whether or not the string has been exhausted.

**qmap_t** contains a queue and a tree map of kmers and scores. This keeps track of the best-scoring
kmers in a given window. You should get its next value by using next_value(kmer, score), which inserts
the kmer and its score into the queue and map and removes the element we're dropping from our map.

**hll_t** is a HyperLogLog implementation which takes hash values and estimates the number
that were inserted after insertion is complete.
We access the value using report(), and we add hash values using the add() function.
We need to call wang_hash on a kmer before adding it.
est_err() returns the estimated error given the theoretical performance of the data structure.
