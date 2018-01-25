from __future__ import division
import fa
import sys
from fa import chunker

if __name__ == "__main__":
    from sys import stderr
    import argparse
    parser = argparse.ArgumentParser(description=(
        "Create a set of synthetic genomes consisting "
        "of subgroups per tax level. Some kmers are unique, "
        "some are shared, and this provides a case where we can test"
        " the efficacy and behavior of our bitmap method."))
    parser.add_argument("-n", "--num-nucleotides-per-leaf",
                        type=int, default=13000)
    parser.add_argument("-N", "--num-nucs-shared-per-subgroup",
                        type=int, default=2000)
    parser.add_argument("-l", "--num-nucs-shared-per-level",
                        type=int, default=8000)
    parser.add_argument("-d", "--tree-depth",
                        type=int, default=4)
    parser.add_argument("-s", "--split-size", type=int,
                        default=3,
                        help=("Number of subgroups for "
                              "each parent node."))
    parser.add_argument("--parent-map", "-p",
                        help="Path to which to write synthetic taxonomy.",
                        default="nodes.dmp")
    parser.add_argument("-S", "--subgroup-size", type=int,
                        default=3,
                        help="Number of genomes for each subgroup")
    parser.add_argument("-o", "--outdir", default="./", type=str)
    parser.add_argument("--name-id-map", "-m", default="synth_nameidmap.txt")
    args = parser.parse_args()

    # Variables/settings for constructing synthetic genome
    # and accessory files.
    mult_per_layer = args.split_size * args.subgroup_size
    depth = args.tree_depth
    nleaves = mult_per_layer ** (depth - 1)
    leaf_seqs = [fa.SeqId(fa.gen_seq(args.num_nucleotides_per_leaf), i) for
                 i in range(nleaves)]
    nleaf_seq = len(leaf_seqs)
    outdir = args.outdir
    name_id_map = outdir + args.name_id_map

    # Variables for constructing the parent_map dictionary.
    pcmap = {}
    used_seqids = set(i.taxid() for i in leaf_seqs)
    ctax = max(used_seqids) + 1
    last_layer = []

    for i in range(1, depth):
        nchunks = nleaf_seq // (mult_per_layer ** i)
        chunk_size = nleaf_seq // nchunks
        assert nleaf_seq % chunk_size == 0
        stderr.write(f"chunk size: {chunk_size}. nchunks: {nchunks}.\n")
        for seqsetid, seqset in enumerate(chunker(leaf_seqs, chunk_size)):
            print("seqset len: %i" % len(seqset), file=stderr)
            add = fa.gen_seq(args.num_nucs_shared_per_level)
            for seq in seqset:
                seq.seq += add
                seq.subsets[i] = seqsetid
            for sssid, seqsubset in enumerate(chunker(seqset,
                                                      args.subgroup_size)):
                # print("seqsubset len: %i" % len(seqsubset), file=stderr)
                add = fa.gen_seq(args.num_nucs_shared_per_subgroup)
                for seq in seqset:
                    seq.seq += add
                    seq.subgroups[i] = seqsetid
            if i == 1:  # or it not last_layer
                # Add leaf node to parent connections
                for seq in seqset:
                    pcmap[seq.taxid()] = ctax + seqsetid
        if i > 1:
            # Add higher nodes to parent connections
            if i == depth - 1:
                pcmap.update((el, 1) for el in last_layer)
                break
                # This leaves the loop on the last layer in the tree
                # because the root is 1 by construction
            else:
                # pcmap.update((tax, i + ctax) for tax in
                #              last_layer[i:i+mult_per_layer] for
                #              i in range(mult_per_layer))
                for i in range(mult_per_layer):
                    for tax in last_layer[i:i + mult_per_layer]:
                        pcmap[tax] = i + ctax
        last_layer = [ctax + i for i in range(nchunks)]
        used_seqids.update(last_layer)
        ctax = max(used_seqids) + 1
    del used_seqids
    del ctax
    del last_layer
    {seq.write(outdir + seq.filename()) for seq in leaf_seqs}
    print("[1/3] Successfully created synthetic genomes.", file=stderr)

    filenames = [outdir + seq.filename() for seq in leaf_seqs]
    fa.write_nameid_map(name_id_map, filenames)
    print("[2/3] Successfully wrote nameidmap to %s." % name_id_map,
          file=stderr)

    fa.write_parent_map(args.parent_map, pcmap)
    print("[3/3] Successfully wrote child->parent map.", file=stderr)
    stderr.write("Genomes: %s\n" % ', '.join(filenames))
    stderr.write("Nameidmap: %s\n" % name_id_map)
    stderr.write("Taxonomy: %s\n" % args.parent_map)
