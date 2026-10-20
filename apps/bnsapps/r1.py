"""R1: easy wins on the R0 profiles (no new data): protocol-bias calibration
and background / contamination handling (apps/LITERATURE.md sections 3-4).

E8 multiplicative bias (McLaren, Willis & Callahan 2019): observed species
   proportions are true proportions times taxon- and protocol-specific
   efficiencies, constant across samples.
   E8a per-species log(ONT / Illumina) bias (centered log ratio) per replicate
       pair: the between-pair SD per species vs the spread of biases across
       species (consistent bias: small within / across ratio).
   E8b calibration, leave-one-ONT-run-out: efficiencies from the other ONT runs
       vs all Illumina D6331 runs, applied to the held-out run; species
       Bray-Curtis to the mean Illumina profile before vs after.
E9 background:
   E9a host background subtraction: expected host-artefact reads per species =
       human-only control rate x the run's read count; subtract (floor 0); species
       Bray-Curtis to the DNA-only profile before vs after, per dilution.
   E9b decontam-style frequency test: per species, slope of log(species reads /
       excess microbial reads) against log10 microbial input over the dilution
       series. True members keep a constant share (slope ~ 0); contaminants and host
       artefacts arrive at a constant absolute level (slope ~ -1)."""
import json
import math

from . import common, mock, r0

PSEUDO = 0.5


def species_counts(acc):
    s = r0.species_reads(common.work("r0", "reports", f"{acc}.tsv"))
    return s, {t: v["reads"] for t, v in s.items() if isinstance(t, int)}


def clr(counts, taxa):
    logs = [math.log(counts.get(t, 0) + PSEUDO) for t in taxa]
    m = sum(logs) / len(logs)
    return [x - m for x in logs]


def bc(a, b, taxa):
    sa, sb = sum(a.get(t, 0) for t in taxa), sum(b.get(t, 0) for t in taxa)
    if sa <= 0 or sb <= 0:
        return float("nan")
    return 1 - sum(min(a.get(t, 0) / sa, b.get(t, 0) / sb) for t in taxa)


def run():
    with open(common.RESULTS + "/r0.json") as fh:
        meta = json.load(fh)["runs"]
    kind = lambda k: sorted(a for a, m in meta.items() if m["kind"] == k)
    taxa = list(mock.MOCK_SPECIES)
    res, lines = {}, ["# R1: bias calibration and background handling on the R0 profiles"]
    # E8
    ill = [species_counts(a)[1] for a in kind("illumina_d6331") + kind("d6331_only")]
    ont = [species_counts(a)[1] for a in kind("ont_d6331")]
    biases = [[x - y for x, y in zip(clr(o, taxa), clr(i, taxa))] for o in ont for i in ill]
    per_taxon_sd = [math.sqrt(sum((b[j] - sum(bb[j] for bb in biases) / len(biases)) ** 2 for b in biases) / (len(biases) - 1))
                    for j in range(len(taxa))]
    mean_bias = [sum(b[j] for b in biases) / len(biases) for j in range(len(taxa))]
    across = math.sqrt(sum((m - sum(mean_bias) / len(mean_bias)) ** 2 for m in mean_bias) / (len(mean_bias) - 1))
    within = sum(per_taxon_sd) / len(per_taxon_sd)
    res["E8a"] = {"taxa": taxa, "mean_log_bias": mean_bias, "per_taxon_sd": per_taxon_sd, "within_sd_mean": within,
                  "across_sd": across, "ratio": within / across}
    lines.append(f"- E8a ONT/Illumina log-efficiency: spread across species SD {common.fmt(across, 2)} (natural log), "
                 f"between replicate pairs mean SD {common.fmt(within, 2)} -> ratio {common.fmt(within / across, 2)} "
                 f"(small = consistent multiplicative bias)")
    before, after = [], []
    ill_mean = {t: sum(i.get(t, 0) / sum(i.values()) for i in ill) / len(ill) for t in taxa}
    for k in range(len(ont)):
        train = [o for j, o in enumerate(ont) if j != k]
        eff = [sum(clr(o, taxa)[j] - clr(i, taxa)[j] for o in train for i in ill) / (len(train) * len(ill)) for j in range(len(taxa))]
        held = ont[k]
        corrected = {t: (held.get(t, 0) + PSEUDO) / math.exp(e) for t, e in zip(taxa, eff)}
        before.append(bc(held, ill_mean, taxa))
        after.append(bc(corrected, ill_mean, taxa))
    res["E8b"] = {"bc_before": before, "bc_after": after}
    mean = lambda v: sum(v) / len(v) if v else float("nan")
    lines.append(f"- E8b leave-one-ONT-run-out calibration: species Bray-Curtis to mean Illumina {common.fmt(mean(before))} "
                 f"-> {common.fmt(mean(after))} (n={len(ont)} held-out runs)")
    # E9a
    hum = [species_counts(a) for a in kind("human_only")]
    rate = {}
    for s, c in hum:
        for t, n in c.items():
            rate[t] = rate.get(t, 0) + n / s["_total"] / len(hum)
    ref = [species_counts(a)[1] for a in kind("d6331_only")]
    ref_mean = {t: sum(r.get(t, 0) / sum(r.values()) for r in ref) / len(ref) for t in taxa}
    e9a = {}
    for acc in kind("dilution"):
        m = meta[acc]
        s, c = species_counts(acc)
        sub = {t: max(0.0, c.get(t, 0) - rate.get(t, 0) * s["_total"]) for t in taxa}
        e9a.setdefault(m["log10_host"], []).append((bc(c, ref_mean, taxa), bc(sub, ref_mean, taxa)))
    res["E9a"] = {str(k): v for k, v in e9a.items()}
    for lg in sorted(e9a):
        v = e9a[lg]
        lines.append(f"- E9a 1:10^{lg}: species Bray-Curtis to DNA-only {common.fmt(mean([x for x, _ in v]))} -> "
                     f"{common.fmt(mean([y for _, y in v]))} after host-background subtraction")
    # E9b frequency test over dilutions 1:1 .. 1:10^3 (microbial excess measurable)
    bg = sum(1 - s["_unclassified"] / s["_total"] for s, _ in hum) / len(hum)
    pts = {}
    for acc in kind("dilution"):
        lg = meta[acc]["log10_host"]
        if lg > 3:
            continue
        s, c = species_counts(acc)
        excess = (1 - s["_unclassified"] / s["_total"] - bg) * s["_total"]
        if excess <= 0:
            continue
        for t, n in c.items():
            pts.setdefault(t, []).append((-lg, math.log((n + PSEUDO) / excess)))
    slopes = {}
    for t, p in pts.items():
        if len({x for x, _ in p}) < 3:
            continue
        mx, my = mean([x for x, _ in p]), mean([y for _, y in p])
        sxx = sum((x - mx) ** 2 for x, _ in p)
        slopes[t] = sum((x - mx) * (y - my) for x, y in p) / sxx / math.log(10)   # per log10 of input
    mock_sl = [v for t, v in slopes.items() if t in mock.MOCK_SPECIES]
    other_sl = [v for t, v in slopes.items() if t not in mock.MOCK_SPECIES]
    res["E9b"] = {"slopes": {str(t): v for t, v in slopes.items()}, "auroc_mock_vs_other": common.auroc(mock_sl, other_sl) if other_sl else float("nan")}
    lines.append(f"- E9b frequency test (slope of log share vs log10 microbial input, 1:1..1:10^3): mock species median "
                 f"{common.fmt(sorted(mock_sl)[len(mock_sl) // 2] if mock_sl else float('nan'), 2)} (n={len(mock_sl)}); other taxa "
                 f"median {common.fmt(sorted(other_sl)[len(other_sl) // 2] if other_sl else float('nan'), 2)} (n={len(other_sl)}); "
                 f"AUROC mock (high slope) vs other {common.fmt(res['E9b']['auroc_mock_vs_other'])}")
    common.write_results("r1", res, lines)
