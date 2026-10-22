"""Gut reference database from UHGG v2.0.2 (MGnify unified human gut genome
catalogue): species ranked by prevalence proxy = number of genomes recovered for
the species and the number of countries they came from (independent samples
across populations); the top N_SPECIES representatives (GTDB lineages) form the
database. Laptop constraint: build memory ~25 B per genome bp for all k-mers, so
the database stores minimizers (window W over k = K; density ~ 2 / (W - K + 2)).

Light steps (overnight-safe): metadata, selection, representative downloads, and a
measured calibration build on CAL_GENOMES species before any full build."""
import csv
import gzip
import os
import subprocess

from . import common

BASE = "https://ftp.ebi.ac.uk/pub/databases/metagenomics/mgnify_genomes/human-gut/v2.0.2"
REFS = os.path.join(common.REPO, "test", "uhgg_refs")   # untracked, like test/zymo
CURL = "/usr/bin/curl"
N_SPECIES = int(os.environ.get("GUTDB_SPECIES", 300))
K, W = 31, int(os.environ.get("GUTDB_WINDOW", 64))
CAL_GENOMES = 50
RANKS = ["superkingdom", "phylum", "class", "order", "family", "genus", "species"]


def metadata():
    path = common.work("gutdb", "genomes-all_metadata.tsv")
    if not common.reuse(path):
        subprocess.run([CURL, "-sSfL", "-o", path + ".part", f"{BASE}/genomes-all_metadata.tsv"], check=True)
        os.replace(path + ".part", path)
    with open(path, newline="") as fh:
        return list(csv.DictReader(fh, delimiter="\t"))


def species_table(rows):
    sp = {}
    for r in rows:
        rep = r["Species_rep"]
        d = sp.setdefault(rep, {"rep": rep, "n_genomes": 0, "countries": set(), "lineage": None, "length": None})
        d["n_genomes"] += 1
        if r.get("Country"):
            d["countries"].add(r["Country"])
        if r["Genome"] == rep:
            d["lineage"], d["length"] = r["Lineage"], int(r["Length"])
    out = [d for d in sp.values() if d["lineage"] and d["lineage"].startswith("d__Bacteria")]
    out.sort(key=lambda d: (-d["n_genomes"], -len(d["countries"]), d["rep"]))
    return out


def rep_path(rep):
    return os.path.join(REFS, f"{rep}.fna.gz")


def download_rep(rep):
    """species_catalogue/<rep[:11]>/<rep>/genome/<rep>.fna, stored gzipped. (The same directory also has
    <rep>.gff and <rep>_InterProScan.tsv: annotations for the taxon x gene-kind tree, X5.)"""
    out = rep_path(rep)
    if os.path.exists(out) and os.path.getsize(out) > 0:
        return 0
    os.makedirs(REFS, exist_ok=True)
    url = f"{BASE}/species_catalogue/{rep[:11]}/{rep}/genome/{rep}.fna"   # e.g. MGYG0000000/MGYG000000001
    data = subprocess.run([CURL, "-sSfL", url], capture_output=True, check=True).stdout
    with gzip.open(out + ".part", "wb", compresslevel=6) as fh:
        fh.write(data)
    os.replace(out + ".part", out)
    return len(data)


def taxonomy(species, name):
    """nodes.dmp / seq2tax from GTDB lineages (d__;p__;...;s__), each representative its own species leaf."""
    ids, parent, rank = {(): 1}, {1: 1}, {1: "no rank"}
    nxt = [4000000]
    def node(path, r):
        if path not in ids:
            ids[path] = nxt[0]; nxt[0] += 1
            parent[ids[path]], rank[ids[path]] = ids[path[:-1]], r
        return ids[path]
    sp_id = {}
    for d in species:
        lin = [x.split("__", 1)[-1] or "unassigned" for x in d["lineage"].split(";")][:7]
        lin += ["unassigned"] * (7 - len(lin))
        t = 1
        for i, r in enumerate(RANKS):
            t = node(tuple(lin[:i]) + ((lin[i] + "|" + d["rep"]) if i == 6 and lin[i] == "unassigned" else lin[i],), r)
        sp_id[d["rep"]] = t
    nodes = common.work("gutdb", name, "nodes.dmp")
    with open(nodes, "w") as fh:
        for c, p in parent.items():
            fh.write(f"{c}\t|\t{p}\t|\t{rank[c]}\t|\n")
    s2t = common.work("gutdb", name, "seq2tax")
    with open(s2t, "w") as fh:
        for d in species:
            with gzip.open(rep_path(d["rep"]), "rt") as gh:
                fh.write(f"{gh.readline()[1:].split()[0]}\t{sp_id[d['rep']]}\n")
    return nodes, s2t, sp_id


def build(name, species, threads=2):
    """Measured minimizer build (/usr/bin/time). Returns (db, nodes, peak RSS bytes)."""
    nodes, s2t, _ = taxonomy(species, name)
    db = common.work("gutdb", name, f"db.k{K}w{W}.cdb")
    rss_path = db + ".rss"
    if common.reuse(db) and os.path.exists(rss_path):
        return db, nodes, int(open(rss_path).read())
    common.memory_gate("build " + name)
    r = subprocess.run(["/usr/bin/time", "-l", common.BONSAI, "build", "-k", str(K), "-w", str(W), "-p", str(threads), "-T", nodes,
                        "-M", s2t, db, *[rep_path(d["rep"]) for d in species]], capture_output=True)
    if r.returncode:
        raise RuntimeError(f"build {name} failed: {r.stderr.decode()[-1500:]}")
    rss = next(int(l.split()[0]) for l in r.stderr.decode().splitlines() if "maximum resident set size" in l)
    with open(rss_path, "w") as fh:
        fh.write(str(rss))
    return db, nodes, rss


def run():
    sp = species_table(metadata())
    top = sp[:N_SPECIES]
    got = sum(download_rep(d["rep"]) for d in top)
    mbp = sum(d["length"] for d in top) / 1e6
    cal = top[:CAL_GENOMES]
    db, nodes, rss = build("cal", cal)
    cal_mbp = sum(d["length"] for d in cal) / 1e6
    per_bp = rss / (cal_mbp * 1e6)
    res = {"n_species_catalogue": len(sp), "n_selected": len(top), "selected_mbp": mbp, "downloaded_bytes": got,
           "k": K, "w": W, "cal_genomes": len(cal), "cal_mbp": cal_mbp, "cal_peak_rss": rss, "rss_per_bp": per_bp,
           "projected_full_rss": per_bp * mbp * 1e6,
           "top": [{"rep": d["rep"], "n_genomes": d["n_genomes"], "countries": len(d["countries"]),
                    "species": d["lineage"].split(";")[-1]} for d in top]}
    lines = ["# Gut reference database (UHGG v2.0.2 representatives, minimizer windows)",
             f"- {len(sp)} bacterial species in the catalogue; top {len(top)} by recovered genomes / countries: {mbp:.0f} Mbp",
             f"- most prevalent: " + ", ".join(f"{d['lineage'].split(';')[-1][3:]} ({d['n_genomes']})" for d in top[:8]),
             f"- calibration build ({len(cal)} species, {cal_mbp:.0f} Mbp, k={K}, w={W}): peak RSS {rss / 1e9:.2f} GB = "
             f"{per_bp:.2f} B/bp -> projected full build {per_bp * mbp / 1e3:.2f} GB"]
    common.write_results("gutdb", res, lines)
