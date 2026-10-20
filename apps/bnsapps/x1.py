"""X1: query vs panel of normals, prototyped on the mock (apps/README.md X1/X2).

Panel: per family (children of a branching node), a Dirichlet-multinomial fit
over N normal profiles: mean child shares mu_c = sum_i x_ic / sum_i n_i and a
family-specific overdispersion rho_v by method of moments (sum of Pearson
residuals^2 = sum (1 + (n_i - 1) rho)). Query: per child, a beta-binomial
predictive e-value (null BB(n, mu_c, rho_v) vs a uniform-share mixture with the
same rho), and the pruned explanation of taxdiff -E (greedy within family against
the remaining siblings, renormalized; tau = e-BH threshold at ALPHA).

Mock data (no download): normals = M4 calibration profiles cal0..cal29 (360
samples, JM109 strain, Dirichlet alpha 50); held-out nulls = cal30..cal59 (360);
positives = M4 B. fragilis x1.5/x2/x4 B-group samples and strain-swap (B766)
B-group samples (60 each).
E11a calibration: fraction of null queries with >= 1 selected node (<= ALPHA).
E11b detection / localization: fraction of shift queries selecting 817 (or its
     lineage), and the mean number of other selections.
E11c strain swap: fraction of swap queries selecting a 562 strain leaf."""
import math
import os

from . import common, mock

ALPHA = 0.05
SEED = 20261028


def node_reads(profile):
    rep = profile[:-6] + ".report.tsv"
    if not common.reuse(rep):
        common.bonsai("profile", "report", profile, rep)
    _, _, rows = common.read_tsv(rep)
    return {int(r["taxid"]): float(r["clade_reads"]) for r in rows if r["taxid"] != "0"}


def families(samples):
    """Branching nodes (>= 2 children seen in the samples) -> their seen children. A lone
    child's counts equal its parent's, so chains need no handling."""
    seen = set()
    for s in samples:
        seen |= {t for t, v in s.items() if v > 0}
    kids = {}
    for t in seen:
        p = mock.PARENT.get(t)
        if p is not None and p != t:
            kids.setdefault(p, []).append(t)
    return {a: sorted(m) for a, m in kids.items() if len(m) >= 2}


def fit_panel(samples, fams):
    panel = {}
    for a, mem in fams.items():
        X = [[s.get(c, 0.0) for c in mem] for s in samples]
        n = [sum(x) for x in X]
        tot = sum(n)
        if tot <= 0:
            continue
        mu = [sum(x[j] for x in X) / tot for j in range(len(mem))]
        num = den = 0.0
        for x, ni in zip(X, n):
            if ni <= 1:
                continue
            for j, m in enumerate(mu):
                if 0 < m < 1:
                    num += (x[j] - ni * m) ** 2 / (ni * m * (1 - m)) - 1.0
                    den += ni - 1.0
        rho = max(0.0, num / den) if den > 0 else 0.0
        panel[a] = {"members": mem, "mu": mu, "rho": min(rho, 0.999)}
    return panel


def log_bb(k, n, p, rho):
    lc = math.lgamma(n + 1) - math.lgamma(k + 1) - math.lgamma(n - k + 1)
    p = min(max(p, 1e-12), 1 - 1e-12)
    if rho <= 1e-12:
        return lc + k * math.log(p) + (n - k) * math.log1p(-p)
    phi = 1 / rho - 1
    a, b = p * phi, (1 - p) * phi
    lb = lambda x, y: math.lgamma(x) + math.lgamma(y) - math.lgamma(x + y)
    return lc + lb(k + a, n - k + b) - lb(a, b)


def log10_e(k, n, p0, rho):
    """Beta-binomial null at p0 vs a uniform-share mixture (same rho): log10 e."""
    if n <= 0:
        return 0.0
    if rho <= 1e-12:
        lalt = -math.log(n + 1)
    else:
        c = (k + .5) / (n + 1)
        sd = math.sqrt(c * (1 - c) * (rho + 1 / (n + 1)))
        lo, hi = max(1e-9, c - 12 * sd), min(1 - 1e-9, c + 12 * sd)
        m = 400
        h = (hi - lo) / m
        vals = [log_bb(k, n, lo + i * h, rho) for i in range(m + 1)]
        mx = max(vals)
        acc = sum((0.5 if i in (0, m) else 1.0) * math.exp(v - mx) for i, v in enumerate(vals))
        lalt = mx + math.log(acc * h)
    return (lalt - log_bb(k, n, p0, rho)) / math.log(10)


def ebh_threshold(es, alpha):
    es = sorted(es, reverse=True)
    m, t = len(es), float("inf")
    for k in range(1, m + 1):
        if es[k - 1] >= m / (alpha * k):
            t = m / (alpha * k)
    return t


def explain(query, panel, alpha=ALPHA):
    """Pruned explanation of a query against the panel: [(taxid, log2 share ratio, log10 e)]."""
    first = []
    for a, f in panel.items():
        n = sum(query.get(c, 0.0) for c in f["members"])
        for c, m in zip(f["members"], f["mu"]):
            if n > 0:
                first.append(10 ** min(300, log10_e(query.get(c, 0.0), n, m, f["rho"])))
    tau = ebh_threshold(first, alpha)
    lt = math.log10(tau) if math.isfinite(tau) else float("inf")
    sel = []
    for a, f in panel.items():
        rem = list(zip(f["members"], f["mu"]))
        while len(rem) >= 2:
            n = sum(query.get(c, 0.0) for c, _ in rem)
            msum = sum(m for _, m in rem)
            if n <= 0 or msum <= 0:
                break
            best = max(rem, key=lambda cm: log10_e(query.get(cm[0], 0.0), n, cm[1] / msum, f["rho"]))
            le = log10_e(query.get(best[0], 0.0), n, best[1] / msum, f["rho"])
            if not le >= lt:
                break
            obs = (query.get(best[0], 0.0) + .5) / (n + 1)
            sel.append((best[0], math.log2(obs / (best[1] / msum)), le))
            rem.remove(best)
    return sel, tau


# ---- global parsimony: sparse multiplicative node effects (tree-guided aggregation) ----

def own_reads(profile):
    """Taxon-level (own) reads per node: the categories of the composition."""
    rep = profile[:-6] + ".report.tsv"
    if not common.reuse(rep):
        common.bonsai("profile", "report", profile, rep)
    _, _, rows = common.read_tsv(rep)
    return {int(r["taxid"]): float(r["taxon_reads"]) for r in rows if r["taxid"] != "0" and float(r["taxon_reads"]) > 0}


def subtree(t, cats):
    """Categories under node t (including t)."""
    out = []
    for c in cats:
        a = c
        while a is not None:
            if a == t:
                out.append(c)
                break
            p = mock.PARENT.get(a)
            a = None if p == a or p is None else p
    return out


class GlobalPanel:
    """Dirichlet-tree panel: per branching node a, a family = its children (seen) plus a pseudo
    child for a's own reads; family-specific concentration alpha_a (method of moments, as
    fit_panel). Category weights q_c (own reads) from the panel; node weight = sum of its
    subtree's category weights times multiplicative effects inherited from selected ancestors."""
    def __init__(self, normals_own, normals_clade):
        self.cats = sorted({c for s in normals_own for c in s})
        tot = {c: sum(s.get(c, 0.0) for s in normals_own) for c in self.cats}
        T = sum(tot.values())
        self.q = {c: (tot[c] + 0.5) / (T + 0.5 * len(self.cats)) for c in self.cats}
        self.nodes = sorted({a for c in self.cats for a in [c] + lineage_of(c)})
        self.sub = {v: subtree(v, self.cats) for v in self.nodes}
        kids = {}
        for v in self.nodes:
            p = mock.PARENT.get(v)
            if p is not None and p != v:
                kids.setdefault(p, []).append(v)
        self.fam = {a: sorted(m) for a, m in kids.items() if len(m) >= 2}
        self.alpha = {}
        for a, mem in self.fam.items():
            X = [[s.get(c, 0.0) for c in mem] + [o.get(a, 0.0)] for s, o in zip(normals_clade, normals_own)]
            n = [sum(x) for x in X]
            totn = sum(n)
            mu = [sum(x[j] for x in X) / totn for j in range(len(mem) + 1)] if totn > 0 else []
            num = den = 0.0
            for x, ni in zip(X, n):
                if ni <= 1:
                    continue
                for j, m in enumerate(mu):
                    if 0 < m < 1:
                        num += (x[j] - ni * m) ** 2 / (ni * m * (1 - m)) - 1.0
                        den += ni - 1.0
            rho = min(0.99, max(1e-6, num / den)) if den > 0 else 1e-6
            self.alpha[a] = 1.0 / rho - 1.0
        self.rho = {a: 1 / (1 + al) for a, al in self.alpha.items()}

    def loglik(self, x_own, x_clade, effects):
        w = dict(self.q)
        for v, b in effects.items():
            for c in self.sub[v]:
                w[c] *= math.exp(b)
        W = {v: sum(w[c] for c in self.sub[v]) for v in self.nodes}
        ll = 0.0
        for a, mem in self.fam.items():
            shares = [W[c] for c in mem] + [w.get(a, 0.0) if a in w else 0.0]
            Z = sum(shares)
            xs = [x_clade.get(c, 0.0) for c in mem] + [x_own.get(a, 0.0)]
            n = sum(xs)
            if n <= 0 or Z <= 0:
                continue
            al = self.alpha[a]
            ll += math.lgamma(al) - math.lgamma(n + al)
            for sh, xx in zip(shares, xs):
                if sh <= 0:
                    continue
                ac = al * sh / Z
                ll += math.lgamma(xx + ac) - math.lgamma(ac)
        return ll


def lineage_of(c):
    out, a = [], mock.PARENT.get(c)
    while a is not None and a != c:
        out.append(a)
        if a == 1:
            break
        nxt = mock.PARENT.get(a)
        if nxt == a:
            break
        a = nxt
    return out


def fit_beta(panel, x, effects, v):
    """x = (own reads, clade reads)."""
    """1-D maximum of the log likelihood in the effect at v (golden section on [-12, 12]),
    and the observed information by finite differences."""
    def f(b):
        e = dict(effects); e[v] = e.get(v, 0.0) + b
        return panel.loglik(x[0], x[1], e)
    lo, hi = -12.0, 12.0
    g = (math.sqrt(5) - 1) / 2
    c, d = hi - g * (hi - lo), lo + g * (hi - lo)
    fc, fd = f(c), f(d)
    for _ in range(60):
        if fc > fd:
            hi, d, fd = d, c, fc
            c = hi - g * (hi - lo); fc = f(c)
        else:
            lo, c, fc = c, d, fd
            d = lo + g * (hi - lo); fd = f(d)
    b = (lo + hi) / 2
    h = 1e-3
    info = -(f(b + h) - 2 * f(b) + f(b - h)) / (h * h)
    return b, f(b) - f(0.0), max(info, 1e-9)


def log10_e_effect(dll, b, info, sigma=2.0):
    """Mixture likelihood ratio over b ~ N(0, sigma^2), Laplace approximation: a valid e-value vs b = 0."""
    return (dll - 0.5 * math.log1p(sigma ** 2 * info) - b * b / (2 * (sigma ** 2 + 1 / info))) / math.log(10)


def explain_global(x, panel, alpha=ALPHA, max_steps=10):
    """Greedy sparse explanation: add the node whose multiplicative effect (inherited by its
    subtree, composition renormalized) has the largest e-value, while e >= tau; tau = the e-BH
    threshold at alpha over the first-step e-values of all nodes."""
    effects, sel = {}, []
    first = {}
    for v in panel.nodes:
        if v == 1:
            continue
        b, dll, info = fit_beta(panel, x, effects, v)
        first[v] = log10_e_effect(dll, b, info)
    tau = ebh_threshold([10 ** min(300, e) for e in first.values()], alpha)
    lt = math.log10(tau) if math.isfinite(tau) else float("inf")
    cand = first
    for _ in range(max_steps):
        v = max(cand, key=cand.get) if cand else None
        if v is None or not cand[v] >= lt:
            break
        b, dll, info = fit_beta(panel, x, effects, v)
        effects[v] = effects.get(v, 0.0) + b
        sel.append((v, b / math.log(2), cand[v], 1.96 / math.sqrt(info) / math.log(2)))
        cand = {}
        for u in panel.nodes:
            if u == 1 or u in effects:
                continue
            bu, dllu, infou = fit_beta(panel, x, effects, u)
            cand[u] = log10_e_effect(dllu, bu, infou)
    return sel, tau


def spike_profile(name, genome_file, taxid, n, seed):
    """A small classified profile of n simulated reads from one genome (the spike)."""
    from . import simreads, stream
    db, nodes, _ = mock.build_db("mock", mock.genomes("mock"))
    prof = common.work("x1", "spikes", f"{name}.bprof")
    if not common.reuse(prof):
        recs = simreads.reads([(mock.mock_path(genome_file), taxid, n)], 150, 0.01, seed)
        stream.profile_stream(simreads.chunks(recs), prof, db, nodes)
    return prof


def merged(null_prof, spike_prof, tag):
    out = common.work("x1", "queries", tag, os.path.basename(null_prof))
    if not common.reuse(out):
        common.bonsai("profile", "merge", "-n", os.path.basename(null_prof)[:-6], out, null_prof, spike_prof)
    return out


def run():
    W = common.work("m4")[:-len("m4")] + "m4"
    prof = lambda tag, g, i: os.path.join(W, tag, f"{g}{i}.bprof")
    panel_p = [prof(f"cal{r}", g, i) for r in range(30) for g in "AB" for i in range(6)]
    null_p = [prof(f"cal{r}", g, i) for r in range(30, 40) for g in "AB" for i in range(6)]    # 120 held out
    gpan = GlobalPanel([own_reads(p) for p in panel_p], [node_reads(p) for p in panel_p])
    fams = families([node_reads(p) for p in panel_p])
    fpanel = fit_panel([node_reads(p) for p in panel_p], fams)
    res = {"seed": SEED, "alpha": ALPHA, "n_panel": len(panel_p), "n_null": len(null_p),
           "dtm_rho": {str(a): r for a, r in gpan.rho.items()}}
    lines = [f"# X1: query vs panel of normals (mock; panel {len(panel_p)} normals, held-out nulls {len(null_p)})",
             "- Dirichlet-tree family rho: " + ", ".join(f"{a}:{r:.3g}" for a, r in sorted(gpan.rho.items())),
             "- positives = held-out nulls + an exactly merged spike profile (classified simulated reads); truth = the spiked taxon"]
    xq = lambda p: (own_reads(p), node_reads(p))
    def evaluate(qs, target, label, true_l2=None):
        g_hit = g_other = f_hit = f_other = 0
        effs = []
        for q in qs:
            sel, _ = explain_global(xq(q), gpan)
            ids = {t for t, _, _, _ in sel}
            g_hit += bool(ids & target); g_other += len(ids - target)
            effs += [e for t, e, _, _ in sel if t in target]
            fsel, _ = explain(node_reads(q), fpanel)
            fids = {t for t, _, _ in fsel}
            f_hit += bool(fids & target); f_other += len(fids - target)
        n = len(qs)
        out = {"global_detect": g_hit / n, "global_other": g_other / n, "family_detect": f_hit / n, "family_other": f_other / n,
               "global_mean_log2_effect": sum(effs) / len(effs) if effs else float("nan"), "true_log2_median": true_l2}
        lines.append(f"- {label}: global detect {common.fmt(out['global_detect'], 2)} (other selections/query "
                     f"{common.fmt(out['global_other'], 2)}, mean log2 effect {common.fmt(out['global_mean_log2_effect'], 2)}"
                     + (f", true median {common.fmt(true_l2, 2)}" if true_l2 is not None else "") +
                     f"); family-greedy detect {common.fmt(out['family_detect'], 2)} (other {common.fmt(out['family_other'], 2)})")
        return out
    # E11a: nulls (any selection = false)
    res["E11a"] = evaluate(null_p, set(), "E11a held-out nulls (any selection is false; pass <= 0.05)")
    res["E11a"]["global_any"] = res["E11a"]["global_other"]
    # E11b: B. fragilis spikes
    bf = {817, 816, 815}
    res["E11b"] = {}
    for n in (500, 2000, 8000):
        sp = spike_profile(f"bfrag{n}", "Bacteroides_fragilis.fasta", 817, n, SEED + n)
        qs = [merged(p, sp, f"bfrag{n}") for p in null_p[:60]]
        l2 = sorted(math.log2((node_reads(p).get(817, 0) + n) / max(1.0, node_reads(p).get(817, 0))) for p in null_p[:60])
        res["E11b"][n] = evaluate(qs, bf, f"E11b B. fragilis +{n} reads", l2[len(l2) // 2])
    # E11c: a new strain (B766) appears in a JM109 sample
    sp = spike_profile("b766_2000", "Escherichia_coli_B766.fasta", 910004, 2000, SEED + 7)
    qs = [merged(p, sp, "b766") for p in null_p[:60]]
    res["E11c"] = evaluate(qs, {910004, 562}, "E11c new strain B766 (+2000 reads)")
    common.write_results("x1", res, lines)
