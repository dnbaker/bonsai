"""The sourmash gather paper's samples (dib-lab/2021-paper-sourmash-gather-pipeline):
podar SRR606249 (Shakya et al. 64-genome synthetic community), Zymo D6300
SRR12324253, Alaska oil reservoir SRR1976948 (hu-s1), iHMP p8808mo11. For each,
the genomes gather matched (GenBank 2020.07) with lineages, and per-genome read
mapping summaries (coverage), fetched from the paper's repository.

A per-sample database of exactly those genomes lets us compare abundance
estimates against gather and mapping on the same reference set. It is circular
for presence (the database is gather's matches), so it evaluates abundance and
statistics, not database search. Taxonomy is built from the lineage strings
(local ids; each genome is its own strain-level leaf), with no NCBI taxdump."""
import csv
import gzip
import io
import os
import re
import subprocess
import urllib.request

from . import common

REPO = "https://raw.githubusercontent.com/dib-lab/2021-paper-sourmash-gather-pipeline/HEAD/gather/inputs"
REFS = os.path.join(common.REPO, "test", "gather_refs")   # untracked, like test/zymo
CURL = "/usr/bin/curl"
RANKS = ["superkingdom", "phylum", "class", "order", "family", "genus", "species"]
# repo sample name -> (ENA run accession, label); p8808mo11 is an HMP2 library name (SRR5650070)
SAMPLES = {"SRR606249": ("SRR606249", "podar"), "SRR12324253": ("SRR12324253", "zymo_d6300"),
           "SRR1976948": ("SRR1976948", "hu_s1"), "p8808mo11": ("SRR5650070", "ihmp_p8808mo11")}


def fetch(url):
    return subprocess.run([CURL, "-sSfL", url], capture_output=True, check=True).stdout


def gather_rows(sample):
    path = common.work("gather", f"{sample}.gather.with-lineages.csv")
    if not common.reuse(path):
        with open(path, "wb") as fh:
            fh.write(fetch(f"{REPO}/gather/{sample}.gather.with-lineages.csv"))
    with open(path) as fh:
        return list(csv.DictReader(fh))


def mapping_rows(sample):
    path = common.work("gather", f"{sample}.mapping.csv")
    if not common.reuse(path):
        with open(path, "wb") as fh:
            fh.write(fetch(f"{REPO}/mapping/{sample}.summary.csv"))
    with open(path) as fh:
        return {r["genome_id"]: r for r in csv.DictReader(fh)}


def accession(row):
    return row["name"].split()[0]


def genome_path(acc):
    return os.path.join(REFS, f"{acc}.fna.gz")


def download_genome(acc):
    """NCBI FTP by accession: list the GCA directory, take the folder for this
    version, fetch *_genomic.fna.gz. Returns bytes downloaded (0 if cached)."""
    out = genome_path(acc)
    if os.path.exists(out) and os.path.getsize(out) > 0:
        return 0
    os.makedirs(REFS, exist_ok=True)
    pre, num = acc.split("_")
    digits = num.split(".")[0]
    base = f"https://ftp.ncbi.nlm.nih.gov/genomes/all/{pre}/{digits[0:3]}/{digits[3:6]}/{digits[6:9]}/"
    listing = fetch(base).decode(errors="replace")
    m = re.search(r'href="(' + re.escape(acc) + r'_[^"/]+)/?"', listing)
    if not m:
        raise RuntimeError(f"{acc}: no assembly folder under {base}")
    folder = m.group(1)
    try:
        data = fetch(f"{base}{folder}/{folder}_genomic.fna.gz")
    except subprocess.CalledProcessError:
        # suppressed assemblies keep their reports but not sequence: fetch the
        # GenBank sequences listed in the assembly report from ENA instead
        report = fetch(f"{base}{folder}/{folder}_assembly_report.txt").decode(errors="replace")
        accs = [l.split("\t")[4] for l in report.splitlines() if l and not l.startswith("#") and len(l.split("\t")) > 4
                and l.split("\t")[4] not in ("na", "")]
        if not accs:
            raise RuntimeError(f"{acc}: no sequence accessions in the assembly report")
        fasta = b"".join(fetch("https://www.ebi.ac.uk/ena/browser/api/fasta/" + ",".join(accs[i:i + 50]))
                         for i in range(0, len(accs), 50))
        if not fasta.startswith(b">"):
            raise RuntimeError(f"{acc}: ENA returned no FASTA")
        data = gzip.compress(fasta)
        print(f"[gather] {acc} suppressed at NCBI: {len(accs)} sequences from ENA", flush=True)
    tmp = out + ".part"
    with open(tmp, "wb") as fh:
        fh.write(data)
    os.replace(tmp, out)
    return len(data)


def taxonomy(sample, rows):
    """nodes.dmp and seq2tax for the sample's genomes from gather's lineages.
    Returns (nodes, seq2tax, genome taxid by accession, species taxid by accession, rank by taxid)."""
    ids, parent, rank = {(): 1}, {1: 1}, {1: "no rank"}
    nxt = [2000000]
    def node(path, r):
        if path not in ids:
            ids[path] = nxt[0]
            nxt[0] += 1
            parent[ids[path]] = ids[path[:-1]]
            rank[ids[path]] = r
        return ids[path]
    genome_id, species_id = {}, {}
    gid = 3000000
    for row in rows:
        lin = [x.strip() or "unassigned" for x in row["lineage"].split(";")][:7]
        lin += ["unassigned"] * (7 - len(lin))
        t = 1
        for i, r in enumerate(RANKS):
            t = node(tuple(lin[:i + 1]), r)
        acc = accession(row)
        genome_id[acc], species_id[acc] = gid, t
        parent[gid], rank[gid] = t, "strain"
        gid += 1
    nodes = common.work("gather", sample, "nodes.dmp")
    with open(nodes, "w") as fh:
        for c, p in parent.items():
            fh.write(f"{c}\t|\t{p}\t|\t{rank[c]}\t|\n")
    s2t = common.work("gather", sample, "seq2tax")
    with open(s2t, "w") as fh:
        for acc, g in genome_id.items():
            with gzip.open(genome_path(acc), "rt") as gh:
                tok = gh.readline()[1:].split()[0]
            fh.write(f"{tok}\t{g}\n")
    return nodes, s2t, genome_id, species_id, rank


def build(sample, threads=4):
    rows, got, missing = [], 0, []
    for r in gather_rows(sample):
        try:
            got += download_genome(accession(r))
            rows.append(r)
        except (RuntimeError, subprocess.CalledProcessError) as e:
            missing.append(accession(r))
            print(f"[gather] {sample}: {accession(r)} unavailable, excluded ({e})", flush=True)
    nodes, s2t, gid, sid, rank = taxonomy(sample, rows)
    db = common.work("gather", sample, "db.k31.cdb")
    if not common.reuse(db):
        common.memory_gate("build " + sample)
        r = subprocess.run(["/usr/bin/time", "-l", common.BONSAI, "build", "-k", "31", "-p", str(threads), "-T", nodes, "-M", s2t, db,
                            *[genome_path(accession(r)) for r in rows]], capture_output=True)
        if r.returncode:
            raise RuntimeError(f"build {sample} failed: {r.stderr.decode()[-1500:]}")
        rss = next((int(l.split()[0]) for l in r.stderr.decode().splitlines() if "maximum resident set size" in l), None)
        with open(db + ".rss", "w") as fh:
            fh.write(str(rss))
        print(f"[gather] built {sample}: peak RSS {rss / 1e9 if rss else float('nan'):.2f} GB", flush=True)
    return {"rows": rows, "db": db, "nodes": nodes, "genome_id": gid, "species_id": sid, "rank": rank, "downloaded_bytes": got,
            "missing": missing}
