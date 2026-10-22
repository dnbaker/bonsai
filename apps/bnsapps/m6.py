"""M6: tree-wide comparison of two samples, every node: abundance (read
fractions) and content (shared / unshared k-mer content, coverage-corrected),
from read-level counts and the counters together (PI 2026-10-20).

Per node v and sample s (profile report + `profile compare -K` over all clades):
  r reads, u distinct k-mers (counter), h k-mer hits, g database k-mers of the clade
  frac  = r / classified reads                              abundance
  P_hat = solve u = P (1 - exp(-h / P)), <= g; g if no duplicates are seen
                                                            content present at v
  f     = u / P_hat                                         fraction of it observed
  C_hat = shared / (fA fB), <= min(P_A, P_B)                content in common
  A-only = P_A - C_hat, B-only = P_B - C_hat;  J_content = C_hat / (P_A + P_B - C_hat)
Endpoints (mock; M4 profiles reused, one new depth pair set):
E6a depth invariance: the same community at 1x vs 4x read depth: J_content ~ 1 at
    nodes with >= MIN_READS in both, while raw Jaccard drifts with the depth ratio.
E6b localization: M4 strain-swap pairs (JM109 vs B766): the lowest-J_content nodes
    are the E. coli lineage; AUROC over nodes (E. coli lineage below family vs the rest).
E6c separation: M4 B. fragilis x4 pairs: the largest |log2 frac ratio| is on 817's
    lineage, while J_content there stays near the null (content unchanged)."""
import math
import os
import random

from . import common, m4, mock, simreads, stream

SEED = 20261026
MIN_READS = 20


def report(profile):
    rep = profile[:-6] + ".report.tsv"
    if not common.reuse(rep):
        common.bonsai("profile", "report", profile, rep)
    _, _, rows = common.read_tsv(rep)
    out = {}
    classified = 0
    for r in rows:
        t = int(r["taxid"])
        if t == 0:
            continue
        classified += int(r["taxon_reads"])
        cov = float(r["cov"]) if r["cov"] != "NA" else 0.0
        u = float(r["clade_uniq_kmers"])
        out[t] = {"r": int(r["clade_reads"]), "u": u, "h": float(r["clade_kmers"]), "g": u / cov if cov > 0 else float("nan"),
                  "parent": int(r["parent"])}
    for v in out.values():
        v["frac"] = v["r"] / classified if classified else 0.0
    return out


def present_content(u, h, g):
    p = common.saturation_universe(u, h)
    if not math.isfinite(p):
        p = g
    return min(p, g) if math.isfinite(g) else p


def node_table(pa, pb):
    """Per node: abundance and coverage-corrected content comparison of two profiles."""
    ra, rb = report(pa), report(pb)
    kpath = common.work("m6", "tmp", "kmers.tsv")
    common.bonsai("profile", "compare", "-K", kpath, pa, pb)
    _, _, rows = common.read_tsv(kpath)
    shared = {int(r["taxid"]): float(r["shared"]) for r in rows}
    out = {}
    for t in set(ra) | set(rb):
        a, b = ra.get(t), rb.get(t)
        row = {"frac_a": a["frac"] if a else 0.0, "frac_b": b["frac"] if b else 0.0,
               "reads_a": a["r"] if a else 0, "reads_b": b["r"] if b else 0}
        if a and b and a["u"] > 0 and b["u"] > 0:
            g = a["g"] if math.isfinite(a["g"]) else b["g"]
            P_a, P_b = present_content(a["u"], a["h"], g), present_content(b["u"], b["h"], g)
            fa, fb = a["u"] / P_a, b["u"] / P_b
            sh = shared.get(t, 0.0)
            C = min(sh / (fa * fb), P_a, P_b) if fa > 0 and fb > 0 else 0.0
            row.update(P_a=P_a, P_b=P_b, C=C, a_only=P_a - C, b_only=P_b - C, J_content=C / (P_a + P_b - C),
                       J_raw=sh / (a["u"] + b["u"] - sh) if a["u"] + b["u"] > sh else 0.0)
        out[t] = row
    return out


def lineage(t):
    out = []
    while t and t in mock.PARENT and t != 1:
        out.append(t)
        t = mock.PARENT[t]
    return out


def well_covered(row):
    return row["reads_a"] >= MIN_READS and row["reads_b"] >= MIN_READS and "J_content" in row


def run():
    rng = random.Random(SEED)
    db, nodes, _ = mock.build_db("mock", mock.genomes("mock"))
    res = {"seed": SEED, "min_reads": MIN_READS}
    lines = ["# M6: tree-wide abundance + content comparison (every node)"]
    # E6a: same community, 1x vs 4x depth
    p = m4.base_props(rng)
    e6a = []
    for i in range(6):
        pi = simreads.dirichlet(rng, [m4.ALPHA * x for x in p])
        profs = []
        for depth in (30000, 120000):
            counts = simreads.multinomial(rng, depth, pi)
            comp = [(mock.mock_path(f), t, c) for (f, t), c in zip(m4.SPECIES, counts[:-1])]
            comp.append((mock.mock_path(m4.ECOLI_FILE[m4.STRAIN_A]), m4.STRAIN_A, counts[-1]))
            prof = common.work("m6", "depth", f"d{i}_{depth}.bprof")
            seed = rng.randrange(1 << 30)
            if not common.reuse(prof):
                stream.profile_stream(simreads.chunks(simreads.reads(comp, m4.READ_LEN, m4.ERR, seed)), prof, db, nodes)
            profs.append(prof)
        tab = node_table(*profs)
        e6a += [(row["J_content"], row["J_raw"]) for row in tab.values() if well_covered(row)]
    med = lambda xs: sorted(xs)[len(xs) // 2] if xs else float("nan")
    res["E6a"] = {"n_nodes": len(e6a), "median_J_content": med([a for a, _ in e6a]), "median_J_raw": med([b for _, b in e6a]),
                  "frac_J_content_ge_0.8": sum(1 for a, _ in e6a if a >= 0.8) / max(1, len(e6a))}
    lines.append(f"- E6a same community, 4x depth ratio ({len(e6a)} well-covered node comparisons): median J_content "
                 f"{common.fmt(res['E6a']['median_J_content'])} vs raw Jaccard {common.fmt(res['E6a']['median_J_raw'])}; "
                 f"J_content >= 0.8 at {common.fmt(res['E6a']['frac_J_content_ge_0.8'], 2)} of nodes (ideal: 1)")
    # E6b: M4 strain swaps, A0 vs B0 of each experiment
    ecoli_lineage = {562, 910004, 910005}
    pos, neg, top_hits = [], [], 0
    for r in range(m4.R_ALT):
        a, b = common.work("m4", f"swap{r}", "A0.bprof"), common.work("m4", f"swap{r}", "B0.bprof")
        if not (os.path.exists(a) and os.path.exists(b)):
            continue
        tab = node_table(a, b)
        cand = {t: row for t, row in tab.items() if well_covered(row) or (t in ecoli_lineage and "J_content" in row)}
        for t, row in cand.items():
            (pos if t in ecoli_lineage else neg).append(1 - row["J_content"])
        lowest = min(cand, key=lambda t: cand[t]["J_content"]) if cand else None
        top_hits += lowest in ecoli_lineage
    res["E6b"] = {"auroc": common.auroc(pos, neg), "lowest_node_in_ecoli_lineage": top_hits, "experiments": m4.R_ALT}
    lines.append(f"- E6b strain swap: content divergence (1 - J_content) AUROC E. coli lineage vs other nodes "
                 f"{common.fmt(res['E6b']['auroc'])}; lowest-J node in the E. coli lineage in {top_hits}/{m4.R_ALT} pairs")
    # E6c: M4 B. fragilis x4: abundance moves on 817's lineage, content does not
    bf_lineage = set(lineage(817))
    top_ab, jc_bf, jc_other = 0, [], []
    for r in range(m4.R_ALT):
        a, b = common.work("m4", f"shift4.0_{r}", "A0.bprof"), common.work("m4", f"shift4.0_{r}", "B0.bprof")
        if not (os.path.exists(a) and os.path.exists(b)):
            continue
        tab = node_table(a, b)
        cand = {t: row for t, row in tab.items() if well_covered(row)}
        lfc = {t: abs(math.log2((row["frac_b"] + 1e-9) / (row["frac_a"] + 1e-9))) for t, row in cand.items() if mock.rank_of(t) == "species"}
        if lfc:
            top_ab += max(lfc, key=lfc.get) == 817
        jc_bf += [cand[t]["J_content"] for t in cand if t in bf_lineage]
        jc_other += [cand[t]["J_content"] for t in cand if t not in bf_lineage]
    res["E6c"] = {"top_abundance_change_is_817": top_ab, "median_J_content_817_lineage": med(jc_bf),
                  "median_J_content_other": med(jc_other)}
    lines.append(f"- E6c abundance shift (B. fragilis x4): largest species |log2 fraction ratio| is 817 in {top_ab}/{m4.R_ALT}; "
                 f"median J_content on its lineage {common.fmt(res['E6c']['median_J_content_817_lineage'])} vs elsewhere "
                 f"{common.fmt(res['E6c']['median_J_content_other'])} (content unchanged: both ~ same)")
    common.write_results("m6", res, lines)
