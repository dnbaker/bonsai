"""Spaced-seed validation on real strain pairs (seeds.md's pending item).

Two strain-pair sets, reads simulated from ONE strain at 1% error, classified
against a database holding every strain of the pair/set as its own leaf under
a shared species: which fraction of reads lands on the correct strain leaf?

  podar:  the S. parasanguinis genome pair gather matched (SRR606249 set)
  mock:   the five D6331 E. coli strains (910001-910005 under species 562)

Patterns (seeds.md): contiguous at the seed weight, and the optimized spaced
patterns at a larger span. The both-strand canonical fix (1d810a6) is
required for the reverse-strand half of the simulated reads to be counted at
all, so each run also checks forward- vs reverse-only reads classify alike.

Pass criterion (PI 2026-10-20): spaced seeds must not lose correct-leaf
fraction vs contiguous at the same weight on real divergence; the point of
the spacing is sensitivity at 2-5% divergence (the sweep's DIVERGENCE sims
cover mutated genomes; here divergence is the strains' real SNPs).
"""
import os
import random
import subprocess

from . import common, mock, simreads, stream

SEED = 20261021
N_READS = 10000
ZYMO = os.path.join(common.REPO, "test", "zymo", "refs", "D6331.refseq", "genomes")
REFS = os.path.join(common.REPO, "test", "gather_refs")
def mask_gaps(mask):
    """0/1 pattern -> bonsai -S gap string (k-1 between-base gap counts)."""
    pos = [i for i, c in enumerate(mask) if c == "1"]
    return ",".join(str(b - a - 1) for a, b in zip(pos, pos[1:]))


PATTERNS = {   # name -> (k weight, span, -S spacing string or None)
    "k21_contig": (21, 21, None),
    "k21_spaced": (21, 27, mask_gaps("10011111111000111100101011010111")),
    "k25_contig": (25, 25, None),
    "k25_spaced": (25, 31, mask_gaps("11110111011101011111110011110111")),
}

SETS = {
    # name -> (genomes [(path, leaf taxid)], species taxid, reads-from leaf, label)
    "podar": ([(os.path.join(REFS, "GCA_009717815.1.fna.gz"), 10011),
               (os.path.join(REFS, "GCA_003470855.1.fna.gz"), 10012)], 1001, 10011,
              "S. parasanguinis BIOML-A3 vs CWD strain pair (podar gather set)"),
    "ecoli": ([(os.path.join(ZYMO, "Escherichia_coli_B1109.fasta"), 910001),
               (os.path.join(ZYMO, "Escherichia_coli_b2207.fasta"), 910002),
               (os.path.join(ZYMO, "Escherichia_coli_B3008.fasta"), 910003),
               (os.path.join(ZYMO, "Escherichia_coli_B766.fasta"), 910004),
               (os.path.join(ZYMO, "Escherichia_coli_JM109.fasta"), 910005)], 562, 910001,
              "D6331 E. coli strains (B1109 reads)"),
}


def write_taxonomy(path, species, leaves):
    rows = [(1, 1, "no rank"), (2, 1, "superkingdom"), (species, 2, "species")] + \
           [(leaf, species, "strain") for leaf in leaves]
    with open(path, "w") as fh:
        for t, p, r in rows:
            fh.write(f"{t}\t|\t{p}\t|\t{r}\t|\n")


def build_db(tag, genomes, species, k, w, spacing):
    nodes = common.work("seedsval", "nodes.dmp")
    write_taxonomy(nodes, species, [t for _, t in genomes])
    s2t = common.work("seedsval", f"{tag}.seq2tax")
    with open(s2t, "w") as fh:
        for path, t in genomes:
            fh.write(f"{mock.first_token(path)}\t{t}\n")
    db = common.work("seedsval", f"{tag}.cdb")
    if not (common.reuse(db) and os.path.exists(db + ".rss")):
        common.memory_gate(f"build {tag}")
        args = [common.BONSAI, "build", "-k", str(k), "-w", str(w)] + (["-S", spacing] if spacing else []) + \
               ["-p", "2", "-T", nodes, "-M", s2t, db] + [g for g, _ in genomes]
        r = subprocess.run(["/usr/bin/time", "-l", *args], capture_output=True)
        if r.returncode:
            raise RuntimeError(r.stderr.decode()[-1500:])
    return db, nodes


def leaf_fraction(report, leaf, species):
    rows = {int(r["taxid"]): r for r in common.read_tsv(report)[2]}
    tot = sum(int(r["taxon_reads"]) for r in rows.values())   # per-read calls, no clade double-count
    correct = int(rows[leaf]["clade_reads"]) if leaf in rows else 0
    spec = int(rows[species]["clade_reads"]) if species in rows else 0
    unc = int(rows[0]["clade_reads"]) if 0 in rows else 0
    return {"correct_leaf": correct / N_READS, "to_species": spec / N_READS,
            "classified": tot / N_READS, "unclassified": unc / N_READS, "n": tot}


def run():
    res = {"seed": SEED, "n_reads": N_READS, "sets": {}}
    lines = ["# Spaced seeds on real strain pairs (seeds.md's pending validation)",
             "| strain set | pattern | reads to correct strain leaf | to species | classified | strand check (fwd vs rev leaf frac) |",
             "|---|---|---|---|---|---|"]
    ok_all = True
    for set_name, (genomes, species, leaf, label) in SETS.items():
        truth = [(p, t, N_READS) for p, t in genomes if t == leaf]
        res["sets"][set_name] = {"label": label, "configs": {}}
        for cfg, (k, w, spacing) in PATTERNS.items():
            tag = f"{set_name}.{cfg}"
            db, nodes = build_db(tag, genomes, species, k, w, spacing)
            row = {}
            for strands, key in ((None, "both"), ("+", "fwd"), ("-", "rev")):
                prof = common.work("seedsval", tag, f"{key}.bprof")
                rep = prof[:-6] + ".report.tsv"
                if not common.reuse(rep):
                    stream.profile_stream(simreads.chunks(simreads.reads(truth, 150, 0.01, SEED, strand=strands)),
                                          prof, db, nodes, uniq=rep)
                row[key] = leaf_fraction(rep, leaf, species)
            strand_ok = abs(row["fwd"]["correct_leaf"] - row["rev"]["correct_leaf"]) < 0.02
            ok_all &= strand_ok
            res["sets"][set_name]["configs"][cfg] = dict(row, strand_ok=strand_ok)
            lines.append(f"| {set_name} | {cfg}{' (spaced)' if spacing else ''} | "
                         f"{row['both']['correct_leaf']:.3f} | {row['both']['to_species']:.3f} | "
                         f"{row['both']['classified']:.3f} (unc {row['both']['unclassified']:.3f}) | "
                         f"{row['fwd']['correct_leaf']:.3f} vs {row['rev']['correct_leaf']:.3f} {'OK' if strand_ok else 'FAIL'} |")
            common.write_results("seeds_validate", res, lines)
    lines.append(f"- {N_READS} simulated 150 bp reads per run at 1% error from the first strain; the database holds every strain as a leaf.")
    lines.append(f"- overall: {'PASS' if ok_all else 'FAIL'} (strand symmetry must hold for every configuration)")
    common.write_results("seeds_validate", res, lines)
    return res
