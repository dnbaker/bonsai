#!/usr/bin/env python3
"""Fixture for the real-data grafting test (round 6): the five D6331
Escherichia coli strain references (chromosome + plasmid records) as a
mini-pangenome under species 562. Writes <outdir>/{taxmap.txt,names.txt,
genomes.list}; genomes.list points at the repo's Zymo references."""
import os

REPO = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
GEN = os.path.join(REPO, "test", "zymo", "refs", "D6331.refseq", "genomes")

STRAINS = [
    ("Escherichia_coli_B1109.fasta", 910001, "B1109"),
    ("Escherichia_coli_b2207.fasta", 910002, "b2207"),
    ("Escherichia_coli_B3008.fasta", 910003, "B3008"),
    ("Escherichia_coli_B766.fasta", 910004, "B766"),
    ("Escherichia_coli_JM109.fasta", 910005, "JM109"),
]
CHAIN = [(1, 0), (2, 1), (1224, 2), (1236, 1224), (91347, 1236), (543, 91347),
         (561, 543), (562, 561)]


def main():
    out = os.path.join(REPO, "data", "ecoli_pangenome")
    os.makedirs(out, exist_ok=True)
    with open(os.path.join(out, "taxmap.txt"), "w") as fh:
        for t, p in CHAIN:
            fh.write("%d\t|\t%d\t|\n" % (t, p))
        for _, t, _ in STRAINS:
            fh.write("%d\t|\t562\t|\n" % t)
    with open(os.path.join(out, "names.txt"), "w") as fh:
        with open(os.path.join(out, "genomes.list"), "w") as gl:
            for f, t, name in STRAINS:
                path = os.path.join(GEN, f)
                with open(path) as gf:
                    tok = gf.readline()[1:].split()[0]
                fh.write("%s\t%d\n" % (tok, t))
                fh.write("%s\t%d\n" % (name, t))
                gl.write(path + "\n")
    print("wrote %s (5 strains, %d taxonomy nodes)" % (out, len(CHAIN) + len(STRAINS) + 1))


if __name__ == "__main__":
    main()
