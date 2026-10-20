"""Database configuration sweep: k x minimizer window w (PI 2026-10-20: default
should be w = k, no minimization; and k = 31 was itself never chosen from data).
Mock database (D6331 + Enterobacterales relatives, 27 genomes) so every build is
light; one configuration at a time.

configs: k in {21, 25, 31} x w in {k (all k-mers), 48, 64}
data: (a) simulated D6331 community, 150 bp at 1% error, with B. fragilis replaced
          by a mutated copy at divergence d in {0, 1, 2, 5}% (exact truth);
      (b) real heads cached once on disk (temporary; deleted after): 2 Illumina D6331,
          2 human-only controls, 1 ONT D6331, plus the local HiFi head.
metrics per config: build peak RSS, database bytes; per input: classified fraction,
classify peak RSS and time; (a) Bray-Curtis of species composition to truth and the
fraction of the divergent species' reads assigned to it; (b) containment AUROC of
true mock species vs taxa seen in the human-only controls, and species Bray-Curtis
Illumina vs HiFi."""
import math
import os
import random
import subprocess

from . import common, m1, mock, simreads, stream

SEED = 20261029
CONFIGS = [(k, w) for k in (21, 25, 31) for w in (k, 48, 64)]
DIVERGENCE = (0.0, 0.01, 0.02, 0.05)
CACHE = os.path.join(common.REPO, "test", "sweep_cache")   # temporary real-read heads
REAL = {"ill1": "SRR33710519", "ill2": "SRR38158100", "hum1": "SRR38158102", "hum2": "SRR38158103", "ont": "ERR14714116"}
HEAD = 60_000_000


def cache_heads():
    from . import r0
    os.makedirs(CACHE, exist_ok=True)
    got = 0
    for tag, acc in REAL.items():
        out = os.path.join(CACHE, f"{acc}.head.fastq.gz")
        if os.path.exists(out) and os.path.getsize(out) > 0:
            continue
        url = r0.ena_search(f"run_accession={acc}")[0]["r1"]
        subprocess.run(["sh", "-c", f"{stream.CURL} -sSfL '{url}' | head -c {HEAD} > '{out}.part'"], check=False)
        os.replace(out + ".part", out)
        got += os.path.getsize(out)
    return got


def mutated_genome(path, d, seed):
    """A copy of a genome with substitutions at rate d (cached by d)."""
    key = f"{os.path.basename(path)}.mut{d}"
    if key not in simreads._GENOMES:
        g = list(simreads.genome(path))
        rng = random.Random(seed)
        alt = {"A": "CGT", "C": "AGT", "G": "ACT", "T": "ACG"}
        if d > 0:
            log1m = math.log1p(-d)
            pos = int(math.log(1 - rng.random()) / log1m)
            while pos < len(g):
                c = g[pos]
                g[pos] = alt[c][rng.randrange(3)] if c in alt else c
                pos += 1 + int(math.log(1 - rng.random()) / log1m)
        simreads._GENOMES[key] = "".join(g)
    return key


def build(k, w):
    nodes = mock.write_taxonomy(common.work("tax", "nodes.dmp"))
    gens = mock.genomes("full")
    s2t = mock.write_seq2tax(common.work("sweep", "full.seq2tax"), gens)
    db = common.work("sweep", f"full.k{k}w{w}.cdb")
    rssf = db + ".rss"
    if not (common.reuse(db) and os.path.exists(rssf)):
        common.memory_gate(f"build k{k} w{w}")
        r = subprocess.run(["/usr/bin/time", "-l", common.BONSAI, "build", "-k", str(k), "-w", str(w), "-p", "2", "-T", nodes,
                            "-M", s2t, db, *[g for g, _ in gens]], capture_output=True)
        if r.returncode:
            raise RuntimeError(r.stderr.decode()[-1500:])
        rss = next(int(l.split()[0]) for l in r.stderr.decode().splitlines() if "maximum resident set size" in l)
        open(rssf, "w").write(str(rss))
    return db, nodes, int(open(rssf).read())


def species_counts(rep):
    _, _, rows = common.read_tsv(rep)
    tot = sum(float(r["taxon_reads"]) for r in rows)
    unc = float(next(r for r in rows if r["taxid"] == "0")["clade_reads"])
    sp = {int(r["taxid"]): r for r in rows if r["rank"] == "species"}
    return tot, unc, sp


def bc(a, b):
    sa, sb = sum(a.values()), sum(b.values())
    if not sa or not sb:
        return float("nan")
    return 1 - sum(min(a.get(t, 0) / sa, b.get(t, 0) / sb) for t in set(a) | set(b))


def run():
    got = cache_heads()
    rng = random.Random(SEED)
    res = {"seed": SEED, "cache_bytes_downloaded": got, "configs": {}}
    lines = ["# Database configuration sweep (k x w; D6331 mock database, 27 genomes)",
             "| k | w | build GB | db MB | sim classified | sim BC to truth | divergent sp. recovery 0/1/2/5% | ill classified | ont classified | hifi classified | containment AUROC true vs host | BC ill vs hifi |",
             "|---|---|---|---|---|---|---|---|---|---|---|---|"]
    base = [(mock.mock_path(f), t, 8000) for f, t in mock.MOCK]
    for k, w in CONFIGS:
        db, nodes, rss = build(k, w)
        c = {"build_rss": rss, "db_bytes": os.path.getsize(db), "sim": {}, "real": {}}
        # (a) simulated, with a divergent B. fragilis
        for d in DIVERGENCE:
            key = mutated_genome(mock.mock_path("Bacteroides_fragilis.fasta"), d, SEED)
            comp = [(p, t, n) if t != 817 else (key, 817, n) for p, t, n in base]
            prof = common.work("sweep", f"k{k}w{w}", f"sim_d{d}.bprof")
            rep = prof[:-6] + ".report.tsv"
            if not common.reuse(rep):
                stream.profile_stream(simreads.chunks(simreads.reads(comp, 150, 0.01, SEED + int(d * 1000))), prof, db, nodes, uniq=rep)
            tot, unc, sp = species_counts(rep)
            truth = {}
            for _, t, n in comp:
                s = 562 if t in mock.STRAINS else t
                truth[s] = truth.get(s, 0) + n
            obs = {t: float(r["clade_reads"]) for t, r in sp.items()}
            c["sim"][d] = {"classified": 1 - unc / tot, "bc_truth": bc(obs, truth), "recovery_817": obs.get(817, 0) / 8000}
        # (b) real heads
        inputs = {tag: os.path.join(CACHE, f"{acc}.head.fastq.gz") for tag, acc in REAL.items()}
        inputs["hifi"] = m1.HIFI
        for tag, path in inputs.items():
            prof = common.work("sweep", f"k{k}w{w}", f"{tag}.bprof")
            rep = prof[:-6] + ".report.tsv"
            if not common.reuse(rep):
                cost = stream.profile_stream(stream.file_source(path), prof, db, nodes, uniq=rep)
                c["real"].setdefault(tag, {})["classify_rss"] = cost["peak_rss_bytes"]
                c["real"][tag]["seconds"] = cost["seconds"]
            tot, unc, sp = species_counts(rep)
            c["real"].setdefault(tag, {})["classified"] = 1 - unc / tot if tot else float("nan")
            c["real"][tag]["species"] = {str(t): {"reads": float(r["clade_reads"]), "contain": float(r["containment"]) if r.get("containment", "NA") != "NA" else 0.0}
                                         for t, r in sp.items()}
        true_c = [c["real"][t]["species"].get(str(s), {}).get("contain", 0.0) for t in ("ill1", "ill2") for s in mock.MOCK_SPECIES]
        host_c = [v["contain"] for t in ("hum1", "hum2") for v in c["real"][t]["species"].values() if v["reads"] > 0]
        c["contain_auroc"] = common.auroc(true_c, host_c) if host_c else float("nan")
        ill = {int(t): v["reads"] for t, v in c["real"]["ill1"]["species"].items()}
        hifi = {int(t): v["reads"] for t, v in c["real"]["hifi"]["species"].items()}
        c["bc_ill_hifi"] = bc(ill, hifi)
        res["configs"][f"k{k}w{w}"] = c
        sim = c["sim"]
        lines.append(f"| {k} | {w}{' (all)' if w == k else ''} | {rss / 1e9:.2f} | {c['db_bytes'] / 1e6:.0f} | "
                     f"{sim[0.0]['classified']:.3f} | {sim[0.0]['bc_truth']:.3f} | "
                     + "/".join(f"{sim[d]['recovery_817']:.2f}" for d in DIVERGENCE) + " | "
                     f"{c['real']['ill1']['classified']:.3f} | {c['real']['ont']['classified']:.3f} | {c['real']['hifi']['classified']:.3f} | "
                     f"{common.fmt(c['contain_auroc'])} | {common.fmt(c['bc_ill_hifi'])} |")
        common.write_results("sweep", res, lines)   # partial results after each config
    lines.append("- divergent sp. recovery: fraction of the mutated B. fragilis reads assigned to species 817 (truth 8000 reads).")
    lines.append("- the read cache in test/sweep_cache is temporary; delete after the sweep.")
    common.write_results("sweep", res, lines)
