"""ZymoBIOMICS D6331 mock (local test/zymo) and absent Enterobacterales
relatives (test/entero): taxonomy, seq2tax and cached databases.

Lineages and the genome table are bench/zymo_holdout.cpp's (NCBI main ranks;
E. coli strains as local leaves 910001-910005 under species 562; the two yeasts
left out to bound memory)."""
import gzip
import json
import os
import resource
import subprocess
import time

from . import common

ZYMO = os.path.join(common.REPO, "test", "zymo", "refs", "D6331.refseq", "genomes")
ENTERO = os.path.join(common.REPO, "test", "entero")

# file, taxid
MOCK = [
    ("Akkermansia_muciniphila.fasta", 239935), ("Bacteroides_fragilis.fasta", 817),
    ("Prevotella_corporis.fasta", 28128), ("Bifidobacterium_adolescentis.fasta", 1680),
    ("Clostridioides_difficile.fasta", 1496), ("Clostridium_perfringens.fasta", 1502),
    ("Faecalibacterium_prausnitzii.fasta", 853), ("Roseburia_hominis.fasta", 301301),
    ("Veillonella_rogosae.fasta", 423477), ("Enterococcus_faecalis.fasta", 1351),
    ("Lactobacillus_fermentum.fasta", 1613), ("Fusobacterium_nucleatum.fasta", 851),
    ("Methanobrevibacter_smithii.fasta", 2173), ("Salmonella_enterica.fasta", 28901),
    ("Escherichia_coli_B1109.fasta", 910001), ("Escherichia_coli_b2207.fasta", 910002),
    ("Escherichia_coli_B3008.fasta", 910003), ("Escherichia_coli_B766.fasta", 910004),
    ("Escherichia_coli_JM109.fasta", 910005),
]
STRAINS = {910001: "B1109", 910002: "B2207", 910003: "B3008", 910004: "B766", 910005: "JM109"}
# absent from the mock (entero/*.fna.gz), species taxid
ABSENT = [("shigella.fna.gz", 623), ("klebsiella.fna.gz", 573), ("citrobacter.fna.gz", 545),
          ("enterobacter.fna.gz", 550), ("proteus.fna.gz", 584), ("morganella.fna.gz", 582),
          ("yersinia.fna.gz", 632), ("serratia.fna.gz", 615)]
MOCK_SPECIES = [239935, 817, 28128, 1680, 1496, 1502, 853, 301301, 423477, 1351, 1613, 851, 2173, 28901, 562]
ABSENT_SPECIES = [t for _, t in ABSENT]

PARENTS = [
    (1, 0), (2, 1), (2157, 1),
    (74201, 2), (203494, 74201), (48461, 203494), (1647988, 48461), (239934, 1647988), (239935, 239934),
    (976, 2), (200643, 976), (171549, 200643), (815, 171549), (816, 815), (817, 816),
    (171552, 171549), (838, 171552), (28128, 838),
    (201174, 2), (1760, 201174), (85004, 1760), (31953, 85004), (1678, 31953), (1680, 1678),
    (1239, 2), (186801, 1239), (186802, 186801),
    (186804, 186802), (1870884, 186804), (1496, 1870884),
    (31979, 186802), (1485, 31979), (1502, 1485),
    (216572, 186802), (216851, 216572), (853, 216851),
    (186803, 186802), (841, 186803), (301301, 841),
    (909932, 1239), (1843489, 909932), (31977, 1843489), (29465, 31977), (423477, 29465),
    (91061, 1239), (186826, 91061), (81852, 186826), (1350, 81852), (1351, 1350),
    (33958, 186826), (2742598, 33958), (1613, 2742598),
    (32066, 2), (203490, 32066), (203491, 203490), (203492, 203491), (848, 203492), (851, 848),
    (28890, 2157), (183925, 28890), (2158, 183925), (2159, 2158), (2172, 2159), (2173, 2172),
    (1224, 2), (1236, 1224), (91347, 1236), (543, 91347),
    (561, 543), (562, 561), (910001, 562), (910002, 562), (910003, 562), (910004, 562), (910005, 562),
    (590, 543), (28901, 590), (620, 543), (623, 620), (570, 543), (573, 570),
    (544, 543), (545, 544), (547, 543), (550, 547),
    (1903414, 91347), (583, 1903414), (584, 583), (581, 1903414), (582, 581),
    (1903411, 91347), (629, 1903411), (632, 629), (613, 1903411), (615, 613),
]
_RANK = {"domain": [2, 2157], "phylum": [74201, 976, 201174, 1239, 32066, 28890, 1224],
         "class": [203494, 200643, 1760, 186801, 909932, 91061, 203490, 183925, 1236],
         "order": [48461, 171549, 85004, 186802, 1843489, 186826, 203491, 2158, 91347],
         "family": [1647988, 815, 171552, 31953, 186804, 31979, 216572, 186803, 31977, 81852, 33958, 203492,
                    2159, 543, 1903414, 1903411],
         "species": MOCK_SPECIES + ABSENT_SPECIES, "strain": list(STRAINS)}
RANK = {t: r for r, ts in _RANK.items() for t in ts}
PARENT = dict(PARENTS)


def rank_of(t):
    if t in RANK:
        return RANK[t]
    return "genus" if any(PARENT.get(c) == t and RANK.get(c) == "species" for c in PARENT) else "no rank"


def mock_path(f):
    return os.path.join(ZYMO, f)


def absent_path(f):
    return os.path.join(ENTERO, f)


def first_token(path):
    op = gzip.open if path.endswith(".gz") else open
    with op(path, "rt") as fh:
        line = fh.readline()
    if not line.startswith(">"):
        raise ValueError(f"{path}: not FASTA")
    tok = line[1:].split()[0]
    if "|" in tok:
        raise ValueError(f"{path}: '|' in the first header; bonsai's get_taxid parses those differently")
    return tok


def write_taxonomy(path):
    with open(path, "w") as fh:
        for c, p in PARENTS:
            fh.write(f"{c}\t|\t{p if c != 1 else 1}\t|\t{'no rank' if c == 1 else rank_of(c).replace('domain', 'superkingdom')}\t|\n")
    return path


def write_seq2tax(path, genomes):
    """genomes: [(path, taxid)]; bonsai maps each file by its first header's first token
    (an unknown name silently becomes taxid 1, so every name is checked to be unique)."""
    seen = {}
    with open(path, "w") as fh:
        for g, t in genomes:
            tok = first_token(g)
            if tok in seen and seen[tok] != t:
                raise ValueError(f"first-header name {tok} of {g} is also another genome's")
            seen[tok] = t
            fh.write(f"{tok}\t{t}\n")
    return path


def genomes(which="full", drop=()):
    g = [(mock_path(f), t) for f, t in MOCK if t not in drop]
    if which == "full":
        g += [(absent_path(f), t) for f, t in ABSENT]
    return g


def build_db(name, gens, k=31, threads=4):
    """Compact (default) database; cached by name. Returns (db, nodes, info)."""
    nodes = write_taxonomy(common.work("tax", "nodes.dmp"))
    db = common.work("db", f"{name}.k{k}.cdb")
    info_path = db + ".info"
    if os.path.exists(db) and os.path.exists(info_path):
        with open(info_path) as fh:
            return db, nodes, json.load(fh)
    s2t = write_seq2tax(common.work("db", f"{name}.seq2tax"), gens)
    common.memory_gate("build " + name)
    t0 = time.time()
    r0 = resource.getrusage(resource.RUSAGE_CHILDREN).ru_maxrss
    common.bonsai("build", "-k", k, "-p", threads, "-T", nodes, "-M", s2t, db, *[g for g, _ in gens])
    info = {"seconds": round(time.time() - t0, 1), "bytes": os.path.getsize(db), "genomes": len(gens), "k": k,
            "max_rss_children": max(r0, resource.getrusage(resource.RUSAGE_CHILDREN).ru_maxrss)}
    with open(info_path, "w") as fh:
        json.dump(info, fh)
    return db, nodes, info


def read_fasta(path):
    """Concatenated sequence of all records (uppercase)."""
    op = gzip.open if path.endswith(".gz") else open
    seq = []
    with op(path, "rt") as fh:
        for line in fh:
            if not line.startswith(">"):
                seq.append(line.strip().upper())
    return "".join(seq)
