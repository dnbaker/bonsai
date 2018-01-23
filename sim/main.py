from __future__ import division
import fa
import sys
from fa import chunker

if __name__ == "__main__":
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
    parser.add_argument("-S", "--subgroup-size", type=int,
                        default=3,
                        help="Number of genomes for each subgroup")
    parser.add_argument("-o", "--outdir", default=None, type=str)
    parser.add_argument("--name-id-map", "-m", default=None, type=str)
    args = parser.parse_args()
    npl = args.num_nucleotides_per_leaf
    mult_per_layer = args.split_size * args.subgroup_size
    depth = args.tree_depth
    nleaves = mult_per_layer ** (depth - 1)
    leaf_seqs = [fa.SeqId(fa.gen_seq(npl), i) for i in range(nleaves)]
    nleaf_seq = len(leaf_seqs)
    outdir = args.outdir if args.outdir else "./"
    name_id_map = args.name_id_map
    if not name_id_map:
        name_id_map = "./synth_nameidmap.txt"
    for i in range(1, depth):
        nchunks = nleaf_seq // (mult_per_layer ** i)
        chunk_size = nleaf_seq // nchunks
        assert nleaf_seq % chunk_size == 0
        sys.stderr.write(f"chunk size: {chunk_size}. nchunks: {nchunks}.\n")
        for seqsetid, seqset in enumerate(chunker(leaf_seqs, chunk_size)):
            print("seqset len: %i" % len(seqset), file=sys.stderr)
            add = fa.gen_seq(args.num_nucs_shared_per_level)
            for seq in seqset:
                seq.seq += add
                seq.subsets[i] = seqsetid
            for sssid, seqsubset in enumerate(chunker(seqset,
                                                      args.subgroup_size)):
                print("seqsubset len: %i" % len(seqsubset), file=sys.stderr)
                add = fa.gen_seq(args.num_nucs_shared_per_subgroup)
                for seq in seqset:
                    seq.seq += add
                    seq.subgroups[i] = seqsetid
    used_seqids = set(i.id for i in leaf_seqs)
    filenames = []
    for seq in leaf_seqs:
        fn = outdir + seq.filename()
        seq.write(fn)
        filenames.append(fn)
    print("Successfully created synthetic genomes.", file=sys.stderr)
    print("I normallly would continue, but I'm killing "
          "this so I can examine logging.")
    sys.exit(1)
    fa.write_nameid_map(name_id_map, filenames)
    with open(args.parent_map, "w") as f:
        raise NotImplementedError("I still need to create the fake taxonomy"
                                  " file. At that point, I'll have a full "
                                  "synthetic taxonomy to use.")
    for i in range(1, depth):
        nnodes = mult_per_layer ** i
        cmax = max(used_seqids)
        print("For tree at height %i..." % i)
    sys.stderr.write("Genomes: %s\n" % ', '.join(filenames))
