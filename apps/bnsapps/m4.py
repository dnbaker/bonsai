"""M4: differential abundance, calibration and power (see apps/README.md).

Replicate mock samples: base species proportions p (log-normal, fixed by seed),
per sample p_i ~ Dirichlet(ALPHA p), reads ~ Multinomial(N, p_i), 150 bp at 1%
error; E. coli is one strain per group. Each experiment: 6 vs 6 samples,
`profile diff -K -r species` (exact enumeration: C(12,6) = 924 relabelings,
two-sided p >= 2/924), plus a Wilcoxon rank-sum on species proportions.
E4a null (R_NULL experiments): FPR at 0.05 per test, KS against uniform.
E4b abundance shift: B. fragilis x f in group B; power = fraction of
    experiments with q < 0.1 for it (read test, Wilcoxon BH).
E4c strain swap: group A JM109, group B B766, equal E. coli proportion; at
    species 562: read-test rejection rate (expect ~ alpha), k-mer power (q < 0.1)."""
import math
import random

from . import common, mock, simreads, stream

SEED = 20261023
N_READS, ALPHA, READ_LEN, ERR = 60000, 50.0, 150, 0.01
N_PER_GROUP = 6
R_NULL, R_ALT = 20, 10
FOLDS = (1.5, 2.0, 4.0)
SHIFT_TAXON = 817            # Bacteroides fragilis
STRAIN_A, STRAIN_B = 910005, 910004   # JM109, B766
SPECIES = [(f, t) for f, t in mock.MOCK if t not in mock.STRAINS]   # 14 species + E. coli (a strain)
ECOLI_FILE = {t: f for f, t in mock.MOCK if t in mock.STRAINS}


def base_props(rng):
    w = [math.exp(rng.gauss(0, 1)) for _ in range(len(SPECIES) + 1)]   # + E. coli
    s = sum(w)
    return [x / s for x in w]


def sample(rng, p, strain, fold=1.0):
    q = list(p)
    if fold != 1.0:
        i = [t for _, t in SPECIES].index(SHIFT_TAXON)
        q[i] *= fold
        s = sum(q)
        q = [x / s for x in q]
    pi = simreads.dirichlet(rng, [ALPHA * x for x in q])
    counts = simreads.multinomial(rng, N_READS, pi)
    comp = [(mock.mock_path(f), t, c) for (f, t), c in zip(SPECIES, counts[:-1])]
    comp.append((mock.mock_path(ECOLI_FILE[strain]), strain, counts[-1]))
    return comp


def experiment(tag, rng, p, db, nodes, strain_b=STRAIN_A, fold=1.0):
    groups = common.work("m4", tag, "groups.tsv")
    paths = []
    with open(groups, "w") as fh:
        for g in ("A", "B"):
            for i in range(N_PER_GROUP):
                comp = sample(rng, p, STRAIN_A if g == "A" else strain_b, fold if g == "B" else 1.0)
                prof = common.work("m4", tag, f"{g}{i}.bprof")
                stream.profile_stream(simreads.chunks(simreads.reads(comp, READ_LEN, ERR, rng.randrange(1 << 30))),
                                      prof, db, nodes)
                fh.write(f"{prof}\t{g}\n")
                paths.append((prof, g))
    _, _, rows = common.bonsai_tsv("profile", "diff", "-K", "-r", "species", "-N", 9999, "-s", rng.randrange(1 << 30), groups)
    diff = {int(r["taxid"]): r for r in rows}
    # Wilcoxon on species proportions
    _, _, trows = common.bonsai_tsv("profile", "table", "-r", "species", *[p for p, _ in paths])
    names = [p.rsplit("/", 1)[1][:-6] for p, _ in paths]
    tot = {n: sum(float(r[n]) for r in trows) for n in names}
    wil = {}
    for r in trows:
        a = [float(r[n]) / tot[n] for n, (_, g) in zip(names, paths) if g == "A"]
        b = [float(r[n]) / tot[n] for n, (_, g) in zip(names, paths) if g == "B"]
        wil[int(r["taxid"])] = common.ranksum_p(a, b)
    wq = dict(zip(wil, common.bh(list(wil.values()))))
    out = {}
    for t, r in diff.items():
        out[t] = {"p": float(r["p"]), "q": float(r["q"]), "kp": float(r["kmer_p"]), "kq": float(r["kmer_q"]),
                  "wp": wil.get(t, float("nan")), "wq": wq.get(t, float("nan")), "lfc": float(r["log2fc"])}
    return out


def rate(xs, thr):
    xs = [x for x in xs if not math.isnan(x)]
    return (sum(1 for x in xs if x < thr) / len(xs)) if xs else float("nan"), len(xs)


def run():
    rng = random.Random(SEED)
    db, nodes, _ = mock.build_db("mock", mock.genomes("mock"))
    p = base_props(rng)
    res = {"seed": SEED, "params": {"n_reads": N_READS, "alpha": ALPHA, "n_per_group": N_PER_GROUP, "r_null": R_NULL,
                                    "r_alt": R_ALT, "folds": FOLDS}, "base_props": p}
    lines = [f"# M4: differential abundance ({N_PER_GROUP} vs {N_PER_GROUP}, {N_READS} reads/sample, Dirichlet alpha {ALPHA})"]
    # E4a null
    null = [experiment(f"null{r}", rng, p, db, nodes) for r in range(R_NULL)]
    e4a = {}
    for key, name in (("p", "reads"), ("kp", "kmers"), ("wp", "wilcoxon")):
        ps = [v[key] for ex in null for v in ex.values()]
        fpr, n = rate(ps, 0.05)
        se = math.sqrt(0.05 * 0.95 / n) if n else float("nan")
        d, ksp = common.ks_uniform(ps)
        e4a[name] = {"fpr": fpr, "n_tests": n, "se_at_0.05": se, "ks_d": d, "ks_p": ksp, "pass": fpr <= 0.05 + 2 * se}
        lines.append(f"- E4a null {name}: FPR {common.fmt(fpr)} over {n} tests (pass <= {common.fmt(0.05 + 2 * se)}): "
                     f"{'PASS' if e4a[name]['pass'] else 'FAIL'}; KS D {common.fmt(d)} (p {common.fmt(ksp)})")
    res["E4a"] = e4a
    # E4b abundance shift
    e4b = {}
    for f in FOLDS:
        exs = [experiment(f"shift{f}_{r}", rng, p, db, nodes, fold=f) for r in range(R_ALT)]
        e4b[f] = {"reads_power": sum(ex[SHIFT_TAXON]["q"] < 0.1 for ex in exs) / R_ALT,
                  "wilcoxon_power": sum(ex[SHIFT_TAXON]["wq"] < 0.1 for ex in exs) / R_ALT,
                  "mean_log2fc": sum(ex[SHIFT_TAXON]["lfc"] for ex in exs) / R_ALT}
        lines.append(f"- E4b B. fragilis x{f}: power reads {common.fmt(e4b[f]['reads_power'], 2)}, Wilcoxon "
                     f"{common.fmt(e4b[f]['wilcoxon_power'], 2)} (q < 0.1, {R_ALT} experiments); mean log2fc "
                     f"{common.fmt(e4b[f]['mean_log2fc'], 2)} (true {common.fmt(math.log2(f), 2)}, before renormalization)")
    res["E4b"] = e4b
    # E4c strain swap
    exs = [experiment(f"swap{r}", rng, p, db, nodes, strain_b=STRAIN_B) for r in range(R_ALT)]
    e4c = {"reads_reject_p05": sum(ex[562]["p"] < 0.05 for ex in exs) / R_ALT,
           "wilcoxon_reject_p05": sum(ex[562]["wp"] < 0.05 for ex in exs) / R_ALT,
           "kmer_power_q01": sum(ex[562]["kq"] < 0.1 for ex in exs) / R_ALT,
           "kmer_p": [ex[562]["kp"] for ex in exs]}
    res["E4c"] = e4c
    lines.append(f"- E4c strain swap (JM109 -> B766, equal abundance) at species 562: read test rejects "
                 f"{common.fmt(e4c['reads_reject_p05'], 2)} (expect ~0.05), Wilcoxon {common.fmt(e4c['wilcoxon_reject_p05'], 2)}; "
                 f"k-mer test power {common.fmt(e4c['kmer_power_q01'], 2)} (q < 0.1; expect >= 0.8)")
    ok = e4a["reads"]["pass"] and e4a["kmers"]["pass"] and e4c["kmer_power_q01"] >= 0.8 and e4c["reads_reject_p05"] <= 0.3
    lines.append(f"- verdict: {'headline holds: calibrated, and the k-mer test sees what abundance cannot' if ok else 'headline does not hold as stated; see the numbers above'}")
    common.write_results("m4", res, lines)
