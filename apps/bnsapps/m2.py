"""M2: false-positive taxa (see apps/README.md).

Present = the 15 mock species; absent = 8 Enterobacterales relatives in the
database but not in the mock. Scores per species: clade reads, clade distinct
k-mers, coverage (distinct / database k-mers of the clade). Endpoint: AUROC
present vs absent (bootstrap 95% CI) and false positives at the threshold that
keeps every detected present species.
Conditions:
  hifi        real HiFi reads (D6331 is staggered: its rarest members may be
              undetectable in a 23k-read head, costing every score alike)
  sim150      simulated 150 bp reads at 1% error, mock genomes uniform
  sim150+novel  as sim150 plus abundant reads of organisms absent from the
              database (Vibrio, Pseudomonas: test/entero, not in it), whose
              hits on conserved k-mers are what makes false-positive taxa; this
              is the condition the scores are meant for. A database holding every
              organism in the sample gives no false-positive taxa to separate.
  staggered+novel  as sim150+novel, present species at log-uniform read counts
              in [50, 20000] (seeded), so the rarest true members have read
              counts like the false positives' (the real D6331 regime)."""
import math
import random

from . import common, m1, mock, simreads, stream

SEED = 20261021


def species_scores(report_tsv):
    _, _, rows = common.read_tsv(report_tsv)
    by = {int(r["taxid"]): r for r in rows if r["taxid"] != "0"}
    out = {}
    for t in mock.MOCK_SPECIES + mock.ABSENT_SPECIES:
        r = by.get(t)
        if r is None:
            out[t] = {"reads": 0.0, "distinct": 0.0, "coverage": 0.0}
        else:
            out[t] = {"reads": float(r["clade_reads"]), "distinct": float(r["clade_uniq_kmers"]),
                      "coverage": float(r["cov"]) if r["cov"] != "NA" else 0.0}
    return out


def evaluate(scores):
    res = {}
    for s in ("reads", "distinct", "coverage"):
        pos = [scores[t][s] for t in mock.MOCK_SPECIES]
        neg = [scores[t][s] for t in mock.ABSENT_SPECIES]
        a = common.auroc(pos, neg)
        lo, hi = common.bootstrap_ci(common.auroc, [pos, neg], seed=SEED)
        det = [v for v in pos if v > 0]
        thr = min(det) if det else float("inf")
        res[s] = {"auroc": a, "ci95": [lo, hi], "fp_at_detected_recall": sum(1 for v in neg if v >= thr),
                  "absent_with_reads": sum(1 for t in mock.ABSENT_SPECIES if scores[t]["reads"] > 0),
                  "present_detected": sum(1 for v in pos if v > 0)}
    return res


def run():
    db, nodes, _ = mock.build_db("full", mock.genomes("full"))
    res = {"seed": SEED}
    lines = ["# M2: false-positive taxa (15 present mock species vs 8 absent relatives)"]
    data = {}
    recs = m1.fastq_records(m1.HIFI)
    data["hifi"] = simreads.chunks(iter(recs))
    comp = [(mock.mock_path(f), t, 20000) for f, t in mock.MOCK]
    data["sim150"] = simreads.chunks(simreads.reads(comp, 150, 0.01, SEED))
    novel = [(mock.absent_path("vibrio.fna.gz"), 0, 100000), (mock.absent_path("pseudomonas.fna.gz"), 0, 100000)]
    data["sim150+novel"] = simreads.chunks(simreads.reads(comp + novel, 150, 0.01, SEED + 1))
    rng = random.Random(SEED + 2)
    stag = [(g, t, int(math.exp(rng.uniform(math.log(50), math.log(20000))))) for g, t, _ in comp]
    res["staggered_reads"] = {str(t): n for _, t, n in stag}
    data["staggered+novel"] = simreads.chunks(simreads.reads(stag + novel, 150, 0.01, SEED + 3))
    for name, src in data.items():
        prof = common.work("m2", name, f"{name}.bprof")
        rep = common.work("m2", name, "report.tsv")
        stream.profile_stream(src, prof, db, nodes, uniq=rep)
        sc = species_scores(rep)
        ev = evaluate(sc)
        res[name] = {"scores": {str(t): v for t, v in sc.items()}, "eval": ev}
        lines.append(f"- {name}:")
        for s, e in ev.items():
            lines.append(f"  - {s}: AUROC {common.fmt(e['auroc'])} (95% CI {common.fmt(e['ci95'][0])}-"
                         f"{common.fmt(e['ci95'][1])}); absent relatives above the lowest detected present species "
                         f"{e['fp_at_detected_recall']}/8; present detected {e['present_detected']}/15")
        lines.append(f"  - absent relatives with >= 1 read: {ev['reads']['absent_with_reads']}/8")
        best = max(ev, key=lambda s: ev[s]["auroc"])
        lines.append(f"  - verdict: coverage {'>=' if ev['coverage']['auroc'] >= ev['reads']['auroc'] else '<'} reads "
                     f"(best: {best}); 15 vs 8 taxa, so a sanity check, not a paper result")
    common.write_results("m2", res, lines)
