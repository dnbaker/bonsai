#!/usr/bin/env python3
"""Generate a module-structured synthetic dataset for grafting_scale_bench.

Hierarchy: DOM x PHYLA x CLASSES x ORDERS x FAMILIES x GENERA x SPECIES x
GENOMES_PER_SPECIES, rooted at taxid 1 (depth-6 nodes are genera, matching
the bench's stage-D scan). Each genome is random background plus:
  - a species block shared by all genomes of a species,
  - a genus module shared by all genomes of a genus,
  - a family module shared by all genomes of alternating families,
so the observed k-mer classes are exactly the species/genus/family sets and
the taxonomy tree's implied masks over-cover (the grafting error signal).

Writes <outdir>/{taxmap.txt, names.txt, genomes.list, genomes/*.fa}.
"""
import argparse
import os
import random


def main() -> None:
    ap = argparse.ArgumentParser()
    ap.add_argument("-o", "--outdir", default="data/syn_scale")
    ap.add_argument("--dom", type=int, default=2)
    ap.add_argument("--phy", type=int, default=4)
    ap.add_argument("--cls", type=int, default=2)
    ap.add_argument("--ord", type=int, default=2, dest="ords")
    ap.add_argument("--fam", type=int, default=2)
    ap.add_argument("--gen", type=int, default=2)
    ap.add_argument("--spe", type=int, default=4)
    ap.add_argument("--per", type=int, default=2, help="genomes per species")
    ap.add_argument("--len", type=int, default=60000, help="bases per genome")
    ap.add_argument("--seed", type=int, default=1234567)
    args = ap.parse_args()
    rng = random.Random(args.seed)

    gdir = os.path.join(args.outdir, "genomes")
    os.makedirs(gdir, exist_ok=True)
    taxf = open(os.path.join(args.outdir, "taxmap.txt"), "w")
    namef = open(os.path.join(args.outdir, "names.txt"), "w")
    listf = open(os.path.join(args.outdir, "genomes.list"), "w")

    next_tax = [10]

    def node(parent: int) -> int:
        t = next_tax[0]
        next_tax[0] += 1
        taxf.write(f"{t}\t|\t{parent}\t|\n")
        return t

    def chunk(seq: str, n: int) -> str:
        return "\n".join(seq[i:i + n] for i in range(0, len(seq), n))

    root = 1
    taxf.write(f"{root}\t|\t0\t|\n")
    # One rogue module per order, shared by the even-species genomes: a
    # union of whole species-pair subtrees that no taxonomy node represents
    # (species nodes are singleton pairs), so grafting adopts them cleanly.
    rogue_modules = {}
    for d in range(args.dom):
        for p in range(args.phy):
            for c in range(args.cls):
                for o in range(args.ords):
                    rogue_modules[(d, p, c, o)] = "".join(rng.choice("ACGT")
                                                          for _ in range(400))
    for d in range(args.dom):
        dtax = node(root)
        for p in range(args.phy):
            ptax = node(dtax)
            for c in range(args.cls):
                ctax = node(ptax)
                for o in range(args.ords):
                    otax = node(ctax)
                    for f in range(args.fam):
                        ftax = node(otax)
                        # Cross-genus module: shared by all genera of this
                        # family on alternating families.
                        fam_module = ("".join(rng.choice("ACGT")
                                               for _ in range(300))
                                      if f % 2 == 0 else None)
                        for g in range(args.gen):
                            gtax = node(ftax)
                            genus_module = "".join(rng.choice("ACGT")
                                                   for _ in range(200))
                            for s in range(args.spe):
                                stax = node(gtax)
                                species_block = "".join(rng.choice("ACGT")
                                                        for _ in range(150))
                                for h in range(args.per):
                                    seqid = f"g{d}_{p}_{c}_{o}_{f}_{g}_{s}_{h}"
                                    # Distinct background per genome; shared
                                    # structure comes only from the blocks, so
                                    # observed k-mer classes are exactly the
                                    # species/genus/family sets and the
                                    # taxonomy tree over-covers. The rogue
                                    # module cuts across families within an
                                    # order (a non-clade set the tree cannot
                                    # represent without grafting).
                                    rogue = rogue_modules[(d, p, c, o)] if s % 2 == 0 else ""
                                    seq = ("".join(rng.choice("ACGT")
                                                   for _ in range(args.len))
                                           + species_block + genus_module
                                           + (fam_module or "") + rogue)
                                    path = os.path.join(gdir, seqid + ".fa")
                                    with open(path, "w") as fh:
                                        fh.write(f">{seqid}\n{chunk(seq, 80)}\n")
                                    namef.write(f"{seqid}\t{stax}\n")
                                    listf.write(path + "\n")
    for fh in (taxf, namef, listf):
        fh.close()
    n = (args.dom * args.phy * args.cls * args.ords * args.fam
         * args.gen * args.spe * args.per)
    print(f"wrote {n} genomes to {args.outdir}")


if __name__ == "__main__":
    main()
