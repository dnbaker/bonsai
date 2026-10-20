"""R0: real mock-community runs streamed from ENA (PI OK 2026-10-20). Each run's
R1 is streamed (head of MAX_BYTES; Illumina files are ordered by flowcell tile,
not organism) into classify -m -P -U against the full database; only profiles
and reports reach disk. Bytes are tallied per run (results r0.json).

E0a dilution series (PRJNA1454942): D6331 in human HCT116 DNA at microbe:host
    1:1 .. 1:10^8 (3 replicates), D6331 DNA-only, human-only and extraction
    blanks. Endpoints: classified read fraction vs nominal ratio (log-log slope
    over the measurable range), mock species detected (>= 1 read; coverage >=
    COV_MIN), Bray-Curtis to the DNA-only profile, and taxa called in the
    human-only and blank controls (reads vs coverage).
E0b platform consistency: species Bray-Curtis among Illumina (PRJNA1268178,
    PRJNA1454942 DNA-only), ONT (PRJEB82315) and the local HiFi head.
E0c D6300 (PRJNA1268178): of its 8 bacteria, E. coli, Salmonella, L. fermentum
    and E. faecalis are in the database; Pseudomonas, Staphylococcus, Listeria
    and Bacillus are not. Present = those 4 species; absent = every other
    database species. AUROC by reads, distinct k-mers and coverage (M2 on real
    reads with real absent organisms)."""
import itertools
import math
import os
import re
import urllib.parse
import urllib.request

from . import common, m1, mock, stream

MAX_BYTES = int(float(os.environ.get("R0_MAX_BYTES", 60e6)))   # ~0.9 MB/s from ENA here: ~1 min per run
ONT_RUNS = int(os.environ.get("R0_ONT_RUNS", 3))
COV_MIN = 1e-3
D6300_IN_DB = [562, 28901, 1613, 1351]


def ena_search(query):
    q = urllib.parse.urlencode({"result": "read_run", "query": query, "format": "tsv", "limit": 500,
                                "fields": "run_accession,instrument_platform,library_layout,read_count,fastq_bytes,fastq_ftp,sample_title"})
    with urllib.request.urlopen("https://www.ebi.ac.uk/ena/portal/api/search?" + q, timeout=120) as r:
        lines = r.read().decode().splitlines()
    head = lines[0].split("\t")
    out = []
    for line in lines[1:]:
        d = dict(zip(head, line.split("\t")))
        urls, sizes = d["fastq_ftp"].split(";"), d["fastq_bytes"].split(";")
        i = 1 if len(urls) == 3 else 0    # paired with an unpaired file first: take _1
        d["r1"], d["r1_bytes"] = "https://" + urls[i], int(sizes[i])
        out.append(d)
    return out


def runs():
    sets = []
    for r in ena_search("study_accession=PRJNA1454942"):
        t = r["sample_title"]
        m = re.search(r"ratio 1:(?:10\^)?(\d+)", t)
        rep = re.search(r"replicate (\d)", t)
        if "DNA only, replicate" in t and "D6331" in t:
            kind, ratio = "d6331_only", 0
        elif m:
            kind, ratio = "dilution", 0 if m.group(1) == "1" and "1:1," in t else int(m.group(1))
        elif "Human DNA only" in t:
            kind, ratio = "human_only", None
        elif "negative control" in t:
            kind, ratio = "blank", None
        else:
            continue
        if kind == "dilution":   # "1:1" -> log10 ratio 0; "1:10^k" -> k
            ratio = 0 if "1:1," in t else int(re.search(r"1:10\^(\d)", t).group(1))
        sets.append(dict(r, kind=kind, log10_host=ratio, rep=int(rep.group(1)) if rep else 0))
    for r in ena_search("study_accession=PRJNA1268178"):
        if "D6331" in r["sample_title"]:
            sets.append(dict(r, kind="illumina_d6331", log10_host=None, rep=int(r["sample_title"][-1])))
        elif "D6300" in r["sample_title"]:
            sets.append(dict(r, kind="d6300", log10_host=None, rep=int(r["sample_title"][-1])))
    ont = sorted((r for r in ena_search("study_accession=PRJEB82315") if r["sample_title"] == "D6331_Zymo"),
                 key=lambda r: r["run_accession"])[:ONT_RUNS]
    sets += [dict(r, kind="ont_d6331", log10_host=None, rep=i + 1) for i, r in enumerate(ont)]
    return sets


def profile_run(r, db, nodes):
    acc = r["run_accession"]
    prof = common.work("r0", "profiles", f"{acc}.bprof")
    rep = common.work("r0", "reports", f"{acc}.tsv")
    if common.reuse(prof) and common.reuse(rep):
        return prof, rep, None
    cost = stream.profile_stream(stream.url_source(r["r1"], MAX_BYTES), prof, db, nodes, uniq=rep, threads=4)
    return prof, rep, cost


def species_reads(rep):
    _, _, rows = common.read_tsv(rep)
    tot = sum(int(r["taxon_reads"]) for r in rows)
    out = {"_total": tot, "_unclassified": int(next(r for r in rows if r["taxid"] == "0")["clade_reads"])}
    for r in rows:
        if r["rank"] == "species":
            out[int(r["taxid"])] = {"reads": int(r["clade_reads"]), "distinct": float(r["clade_uniq_kmers"]),
                                    "coverage": float(r["cov"]) if r["cov"] != "NA" else 0.0}
    return out


def bc(a, b):
    ka = {t: v["reads"] for t, v in a.items() if isinstance(t, int)}
    kb = {t: v["reads"] for t, v in b.items() if isinstance(t, int)}
    sa, sb = sum(ka.values()), sum(kb.values())
    if not sa or not sb:
        return float("nan")
    return 1 - sum(min(ka.get(t, 0) / sa, kb.get(t, 0) / sb) for t in set(ka) | set(kb))


def run():
    db, nodes, _ = mock.build_db("full", mock.genomes("full"))
    sets = runs()
    res = {"max_bytes": MAX_BYTES, "runs": {}, "bytes_streamed": 0}
    sp = {}
    for r in sets:
        prof, rep, cost = profile_run(r, db, nodes)
        sp[r["run_accession"]] = species_reads(rep)
        res["runs"][r["run_accession"]] = {k: r[k] for k in ("kind", "log10_host", "rep", "sample_title", "r1_bytes")}
        if cost:
            res["runs"][r["run_accession"]]["cost"] = cost
            res["bytes_streamed"] += cost["bytes"]
        print(f"[r0] {r['run_accession']} {r['kind']} {r['log10_host']} rep{r['rep']}: "
              f"{sp[r['run_accession']]['_total']} reads", flush=True)
    lines = [f"# R0: real mock runs streamed from ENA (R1 head <= {MAX_BYTES / 1e6:.0f} MB each)",
             f"- {len(sets)} runs; streamed this invocation: {res['bytes_streamed'] / 1e9:.2f} GB (resumed runs not re-streamed)"]
    by = lambda kind: [r for r in sets if r["kind"] == kind]
    # E0a
    ref = [sp[r["run_accession"]] for r in by("d6331_only")]
    e0a = []
    for r in sorted(by("dilution"), key=lambda r: (r["log10_host"], r["rep"])):
        s = sp[r["run_accession"]]
        frac = 1 - s["_unclassified"] / s["_total"] if s["_total"] else float("nan")
        det = sum(1 for t in mock.MOCK_SPECIES if t in s and s[t]["reads"] > 0)
        detc = sum(1 for t in mock.MOCK_SPECIES if t in s and s[t]["coverage"] >= COV_MIN)
        fp = sum(1 for t in mock.ABSENT_SPECIES if t in s and s[t]["reads"] > 0)
        b = sum(bc(s, x) for x in ref) / len(ref) if ref else float("nan")
        e0a.append({"run": r["run_accession"], "log10_host": r["log10_host"], "rep": r["rep"], "classified_frac": frac,
                    "detected_reads": det, "detected_cov": detc, "absent_with_reads": fp, "bc_to_dna_only": b})
    res["E0a"] = e0a
    pts = [(e["log10_host"], math.log10(e["classified_frac"])) for e in e0a if e["classified_frac"] > 0]
    measurable = [p for p in pts if p[1] > -5.5]
    if len(measurable) >= 3:
        mx = sum(x for x, _ in measurable) / len(measurable)
        my = sum(y for _, y in measurable) / len(measurable)
        slope = sum((x - mx) * (y - my) for x, y in measurable) / sum((x - mx) ** 2 for x, _ in measurable)
    else:
        slope = float("nan")
    res["E0a_slope"] = slope
    lines.append(f"- E0a dilution: log10(classified fraction) vs log10(host:microbe) slope {common.fmt(slope, 2)} "
                 f"over {len(measurable)} measurable runs (ideal -1)")
    for lg in sorted({e["log10_host"] for e in e0a}):
        g = [e for e in e0a if e["log10_host"] == lg]
        lines.append(f"  - 1:10^{lg}: classified {', '.join('%.2e' % e['classified_frac'] for e in g)}; species "
                     f"detected (reads) {[e['detected_reads'] for e in g]}, (cov >= {COV_MIN}) {[e['detected_cov'] for e in g]}; "
                     f"absent relatives with reads {[e['absent_with_reads'] for e in g]}; BC to DNA-only "
                     f"{[round(e['bc_to_dna_only'], 3) for e in g]}")
    for kind in ("human_only", "blank"):
        for r in by(kind):
            s = sp[r["run_accession"]]
            called = {t: v for t, v in s.items() if isinstance(t, int) and v["reads"] > 0}
            lines.append(f"  - {kind} {r['run_accession']}: {s['_total']} reads, {s['_total'] - s['_unclassified']} classified; "
                         f"species with reads {len(called)}, with coverage >= {COV_MIN} "
                         f"{sum(1 for v in called.values() if v['coverage'] >= COV_MIN)}")
            res.setdefault("E0a_controls", {})[r["run_accession"]] = {str(t): v for t, v in called.items()}
    # E0b
    groups = {"illumina": by("illumina_d6331") + by("d6331_only"), "ont": by("ont_d6331")}
    hifi_rep = common.work("r0", "reports", "hifi_local.tsv")
    if not common.reuse(hifi_rep):
        stream.profile_stream(stream.file_source(m1.HIFI), common.work("r0", "profiles", "hifi_local.bprof"), db, nodes, uniq=hifi_rep)
    hifi = species_reads(hifi_rep)
    e0b = {}
    names = {g: [r["run_accession"] for r in rs] for g, rs in groups.items()}
    for g, ns in names.items():
        e0b[f"within_{g}"] = [bc(sp[a], sp[b]) for a, b in itertools.combinations(ns, 2)]
    e0b["illumina_vs_ont"] = [bc(sp[a], sp[b]) for a in names["illumina"] for b in names["ont"]]
    e0b["illumina_vs_hifi"] = [bc(sp[a], hifi) for a in names["illumina"]]
    e0b["ont_vs_hifi"] = [bc(sp[a], hifi) for a in names["ont"]]
    res["E0b"] = e0b
    mean = lambda v: sum(v) / len(v) if v else float("nan")
    lines.append("- E0b species Bray-Curtis (mean): " + "; ".join(f"{k} {common.fmt(mean(v))} (n={len(v)})" for k, v in e0b.items()))
    # E0c
    e0c = {}
    for r in by("d6300"):
        s = sp[r["run_accession"]]
        scores = {t: s.get(t, {"reads": 0, "distinct": 0.0, "coverage": 0.0}) for t in mock.MOCK_SPECIES + mock.ABSENT_SPECIES}
        ev = {}
        for sc in ("reads", "distinct", "coverage"):
            pos = [scores[t][sc] for t in D6300_IN_DB]
            neg = [scores[t][sc] for t in scores if t not in D6300_IN_DB]
            ev[sc] = {"auroc": common.auroc(pos, neg), "absent_above_min_present": sum(1 for v in neg if v >= min(pos))}
        e0c[r["run_accession"]] = {"eval": ev, "scores": {str(t): v for t, v in scores.items()}}
        lines.append(f"- E0c D6300 {r['run_accession']}: " + "; ".join(
            f"{sc} AUROC {common.fmt(v['auroc'])} ({v['absent_above_min_present']} absent >= weakest present)" for sc, v in ev.items()))
    res["E0c"] = e0c
    common.write_results("r0", res, lines)
