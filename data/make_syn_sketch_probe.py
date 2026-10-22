#!/usr/bin/env python3
"""Generate the round-4 sketch-only-naming probe: one genus, `fanout`
species x 2 genomes, random background plus a cross-species module block
carried by all genomes of the first half of the species (the seeded-hull
scenario: the exact class exists, but the sketch-only lane must NAME it
from bbminhash intersection cardinalities alone).

Writes <outdir>/{taxmap.txt, names.txt, genomes.list, genomes/*.fa}.
"""
import argparse
import os
import random


def main() -> None:
    ap = argparse.ArgumentParser()
    ap.add_argument("-o", "--outdir", default="data/syn_sketch_probe")
    ap.add_argument("--fanout", type=int, default=24)
    ap.add_argument("--per", type=int, default=2)
    ap.add_argument("--len", type=int, default=30000, help="random background bases")
    ap.add_argument("--module-len", type=int, default=120)
    ap.add_argument("--seed", type=int, default=424242)
    args = ap.parse_args()
    rng = random.Random(args.seed)

    gdir = os.path.join(args.outdir, "genomes")
    os.makedirs(gdir, exist_ok=True)
    taxf = open(os.path.join(args.outdir, "taxmap.txt"), "w")
    namef = open(os.path.join(args.outdir, "names.txt"), "w")
    listf = open(os.path.join(args.outdir, "genomes.list"), "w")
    taxf.write("1\t|\t0\t|\n2\t|\t1\t|\n50\t|\t2\t|\n")

    module = "".join(rng.choice("ACGT") for _ in range(args.module_len))

    def genome(name: str, carry: bool) -> None:
        seq = "".join(rng.choice("ACGT") for _ in range(args.len))
        if carry:
            seq += module
        path = os.path.join(gdir, name + ".fa")
        with open(path, "w") as fh:
            fh.write(">%s\n%s\n" % (name, seq))
        listf.write(path + "\n")

    half = args.fanout // 2
    for s in range(args.fanout):
        sp = 100 + s
        taxf.write("%d\t|\t50\t|\n" % sp)
        for g in range(args.per):
            name = "sp%02dg%d" % (s, g)
            namef.write("%s\t%d\n" % (name, sp))
            genome(name, carry=s < half)

    taxf.close(); namef.close(); listf.close()
    print("wrote %d genomes to %s (module %d bp in species 0..%d)" %
          (args.fanout * args.per, args.outdir, args.module_len, half - 1))


if __name__ == "__main__":
    main()
