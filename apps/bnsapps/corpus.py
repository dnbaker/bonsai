"""Corpus selection from curatedMetagenomicData's curated metadata
(waldronlab/curatedMetagenomicDataCuration, inst/curated/<study>/<study>_sample.tsv).
Metadata only: no reads are fetched here.

Curated schema (2024+): body_site 'feces' (or legacy 'stool'), disease 'Healthy',
age_group 'Adult', ncbi_accession; CRC studies have 'Colorectal Carcinoma' in
target_condition; cases: disease mentions colorectal carcinoma / cancer
(adenomas excluded); controls: disease 'Healthy'.
normals: healthy adult stool, with run accessions, at most PER_STUDY per study
         (seeded choice), from >= MIN_STUDIES studies: the panel of normals.
crc:     cases and controls from studies with a CRC arm, for the no-training utility
         test (each sample queried alone against the panel; AUROC CRC vs control).
Output: apps/results/corpus_<set>.tsv with study, sample, condition, country, runs
(';'), and the R1 URL and bytes of the first run (ENA)."""
import csv
import io
import json
import os
import random
import subprocess
import urllib.parse

from . import common

TREE = "https://api.github.com/repos/waldronlab/curatedMetagenomicDataCuration/git/trees/HEAD?recursive=1"
RAW = "https://raw.githubusercontent.com/waldronlab/curatedMetagenomicDataCuration/HEAD/"
CURL = "/usr/bin/curl"
SEED = 20261027
PER_STUDY = int(os.environ.get("CORPUS_PER_STUDY", 40))
N_NORMALS = int(os.environ.get("CORPUS_NORMALS", 300))
MIN_STUDIES = 5


def fetch(url):
    return subprocess.run([CURL, "-sSfL", url], capture_output=True, check=True).stdout


def study_tables():
    tree = json.loads(fetch(TREE))["tree"]
    paths = [x["path"] for x in tree if x["type"] == "blob" and x["path"].startswith("inst/curated/")
             and x["path"].endswith("_sample.tsv")]
    out = {}
    for p in paths:
        cache = common.work("corpus", "meta", os.path.basename(p))
        if not common.reuse(cache):
            with open(cache, "wb") as fh:
                fh.write(fetch(RAW + p))
        with open(cache, newline="") as fh:
            out[p.split("/")[2]] = list(csv.DictReader(fh, delimiter="\t"))
    return out


def runs_of(row):
    acc = row.get("ncbi_accession") or row.get("NCBI_accession") or ""
    return [a for a in acc.replace(",", ";").split(";") if a.strip().startswith(("SRR", "ERR", "DRR"))]


def ena_r1(runs):
    """R1 https URL and bytes for each run (one ENA search per 100 runs)."""
    out = {}
    for i in range(0, len(runs), 100):
        q = " OR ".join(f"run_accession={r}" for r in runs[i:i + 100])
        url = ("https://www.ebi.ac.uk/ena/portal/api/search?" +
               urllib.parse.urlencode({"result": "read_run", "query": q, "format": "tsv", "limit": 0,
                                       "fields": "run_accession,fastq_ftp,fastq_bytes,instrument_platform"}))
        lines = fetch(url).decode().splitlines()
        head = lines[0].split("\t")
        for line in lines[1:]:
            d = dict(zip(head, line.split("\t")))
            u, b = d.get("fastq_ftp", "").split(";"), d.get("fastq_bytes", "").split(";")
            if not u or not u[0]:
                continue
            j = 1 if len(u) == 3 else 0
            out[d["run_accession"]] = ("https://" + u[j], int(b[j]) if b[j] else 0, d.get("instrument_platform", ""))
    return out


def write(name, rows):
    path = os.path.join(common.RESULTS, f"corpus_{name}.tsv")
    os.makedirs(common.RESULTS, exist_ok=True)
    cols = ["study", "sample_id", "condition", "country", "age_category", "runs", "r1_url", "r1_bytes", "platform"]
    with open(path, "w", newline="") as fh:
        w = csv.DictWriter(fh, cols, delimiter="\t", extrasaction="ignore")
        w.writeheader()
        w.writerows(rows)
    return path


def run():
    rng = random.Random(SEED)
    tabs = study_tables()
    normals_by_study, crc = {}, []
    for study, rows in tabs.items():
        has_crc = any("colorectal carcinoma" in (r.get("target_condition") or "").lower() for r in rows)
        is_crc = lambda d: "colorectal carcinoma" in d.lower() or "colorectal cancer" in d.lower()
        for r in rows:
            runs = runs_of(r)
            if not runs or (r.get("body_site") or "") not in ("feces", "stool"):
                continue
            dis = (r.get("disease") or "").strip()
            rec = {"study": study, "sample_id": r.get("sample_id", ""), "condition": dis, "country": r.get("country", ""),
                   "age_category": r.get("age_group", ""), "runs": ";".join(runs)}
            if dis == "Healthy" and (r.get("age_group") or "") == "Adult":
                normals_by_study.setdefault(study, []).append(rec)
            if has_crc and (dis == "Healthy" or is_crc(dis)):
                crc.append(dict(rec, condition="CRC" if is_crc(dis) else "control"))
    # normals: cap per study, then fill to N_NORMALS across studies (seeded)
    picked = []
    for study in sorted(normals_by_study):
        rows = normals_by_study[study][:]
        rng.shuffle(rows)
        picked += rows[:PER_STUDY]
    rng.shuffle(picked)
    picked = picked[:N_NORMALS]
    crc_ids = {(r["study"], r["sample_id"]) for r in crc}
    normals = [r for r in picked if (r["study"], r["sample_id"]) not in crc_ids]   # keep CRC controls out of the panel
    meta = ena_r1(sorted({r["runs"].split(";")[0] for r in normals + crc}))
    def attach(rows):
        out = []
        for r in rows:
            m = meta.get(r["runs"].split(";")[0])
            if m:
                out.append(dict(r, r1_url=m[0], r1_bytes=m[1], platform=m[2]))
        return out
    normals, crc = attach(normals), attach(crc)
    pn, pc = write("normals", normals), write("crc", crc)
    studies = sorted({r["study"] for r in normals})
    res = {"seed": SEED, "per_study": PER_STUDY, "n_normals": len(normals), "normal_studies": studies,
           "n_crc_cases": sum(1 for r in crc if r["condition"] == "CRC"), "n_crc_controls": sum(1 for r in crc if r["condition"] == "control"),
           "crc_studies": sorted({r["study"] for r in crc}), "platforms": sorted({r["platform"] for r in normals + crc})}
    lines = ["# Corpus selection (curatedMetagenomicData curated metadata; metadata only)",
             f"- normals: {len(normals)} healthy adult stool samples from {len(studies)} studies (<= {PER_STUDY} per study; "
             f"CRC-study controls excluded): {pn}",
             f"- CRC test set: {res['n_crc_cases']} cases, {res['n_crc_controls']} controls from {len(res['crc_studies'])} studies: {pc}",
             f"- platforms: {res['platforms']}",
             f"- {'PASS' if len(studies) >= MIN_STUDIES else 'FAIL'}: >= {MIN_STUDIES} studies in the panel"]
    common.write_results("corpus", res, lines)
