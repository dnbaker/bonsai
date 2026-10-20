"""M3: strain resolution with the five D6331 E. coli strains (see apps/README.md).

Samples: 1-2 E. coli strains at per-strain coverage c on a mock background
(1500 reads per other mock species), 150 bp reads at 1% error. Pair truth: the
samples share >= 1 strain. Scores: k-mer Jaccard of clade 562 (profile compare
-K) and 1 - Bray-Curtis of strain-leaf read counts (profile table -r strain).
E3a: strains in the database (mock database).
E3b: leave one strain out of the database (5 folds); pairs where >= 1 sample
     carries the left-out strain, pooled over folds.
Endpoint: AUROC (bootstrap 95% CI over pairs) per coverage; detection limit =
lowest coverage with AUROC >= 0.95."""
import itertools
import random

from . import common, mock, simreads, stream

SEED = 20261022
COVERAGES = (0.1, 0.3, 1.0, 3.0)
N_PER_COV = 16
READ_LEN, ERR = 150, 0.01
BACKGROUND_READS = 1500
ECOLI = [(f, t) for f, t in mock.MOCK if t in mock.STRAINS]
OTHERS = [(f, t) for f, t in mock.MOCK if t not in mock.STRAINS]


def sample_plan(rng, cov):
    k = rng.choice((1, 2))
    strains = rng.sample([t for _, t in ECOLI], k)
    comp = [(mock.mock_path(f), t, BACKGROUND_READS) for f, t in OTHERS]
    for f, t in ECOLI:
        if t in strains:
            n = max(1, round(cov * len(simreads.genome(mock.mock_path(f))) / READ_LEN))
            comp.append((mock.mock_path(f), t, n))
    return strains, comp


def pair_scores(profiles, names):
    """k-mer Jaccard at clade 562 and 1 - Bray-Curtis over strain leaves, per pair."""
    kpath = common.work("m3", "tmp", "kmers.tsv")
    common.bonsai("profile", "compare", "-r", "species", "-K", kpath, *profiles)
    _, _, krows = common.read_tsv(kpath)
    jac = {}
    for r in krows:
        if r["taxid"] == "562":
            jac[frozenset((r["sample_a"], r["sample_b"]))] = float(r["jaccard"])
    _, header, trows = common.bonsai_tsv("profile", "table", "-r", "strain", *profiles)
    leaves = {n: {} for n in names}
    for r in trows:
        for n in names:
            leaves[n][r["taxid"]] = float(r[n])
    out = {}
    for a, b in itertools.combinations(names, 2):
        va, vb = leaves[a], leaves[b]
        sa, sb = sum(va.values()), sum(vb.values())
        if sa > 0 and sb > 0:
            bc = 1 - sum(min(va.get(t, 0) / sa, vb.get(t, 0) / sb) for t in set(va) | set(vb))
        else:
            bc = 1.0
        out[(a, b)] = {"kmer_jaccard": jac.get(frozenset((a, b)), 0.0), "leaf_abundance_sim": 1 - bc}
    return out


def evaluate(pairs, truth):
    res = {}
    for m in ("kmer_jaccard", "leaf_abundance_sim"):
        pos = [v[m] for p, v in pairs.items() if truth[p]]
        neg = [v[m] for p, v in pairs.items() if not truth[p]]
        res[m] = {"auroc": common.auroc(pos, neg), "ci95": common.bootstrap_ci(common.auroc, [pos, neg], n=1000, seed=SEED),
                  "n_pos": len(pos), "n_neg": len(neg)}
    return res


def run():
    rng = random.Random(SEED)
    dbs = {"mock": mock.build_db("mock", mock.genomes("mock"))}
    for t, s in mock.STRAINS.items():
        dbs[f"no_{s}"] = mock.build_db(f"mock_no_{s}", mock.genomes("mock", drop=(t,)))
    res = {"seed": SEED, "coverages": COVERAGES, "n_per_cov": N_PER_COV, "E3a": {}, "E3b": {}}
    lines = ["# M3: strain resolution (5 E. coli strains)",
             f"- {N_PER_COV} samples per coverage, 1-2 strains each, {READ_LEN} bp at {ERR:.0%} error; "
             "AUROC for 'pair shares >= 1 strain'"]
    for cov in COVERAGES:
        plans, profs = {}, {k: [] for k in dbs}
        for i in range(N_PER_COV):
            name = f"c{cov}_s{i}"
            strains, comp = sample_plan(rng, cov)
            plans[name] = strains
            data = list(simreads.chunks(simreads.reads(comp, READ_LEN, ERR, rng.randrange(1 << 30))))
            for k, (db, nodes, _) in dbs.items():
                p = common.work("m3", k, f"c{cov}", f"{name}.bprof")
                stream.profile_stream(iter(data), p, db, nodes)
                profs[k].append(p)
        names = list(plans)
        truth = {(a, b): bool(set(plans[a]) & set(plans[b])) for a, b in itertools.combinations(names, 2)}
        res["E3a"][cov] = evaluate(pair_scores(profs["mock"], names), truth)
        pooled, ptruth = {}, {}
        for t, s in mock.STRAINS.items():
            sc = pair_scores(profs[f"no_{s}"], names)
            for (a, b), v in sc.items():
                if t in plans[a] or t in plans[b]:
                    pooled[(s, a, b)] = v
                    ptruth[(s, a, b)] = truth[(a, b)]
        res["E3b"][cov] = evaluate(pooled, ptruth)
        e, f = res["E3a"][cov], res["E3b"][cov]
        lines.append(f"- c={cov}x: E3a k-mer {common.fmt(e['kmer_jaccard']['auroc'])} "
                     f"[{common.fmt(e['kmer_jaccard']['ci95'][0])}, {common.fmt(e['kmer_jaccard']['ci95'][1])}] vs leaf "
                     f"{common.fmt(e['leaf_abundance_sim']['auroc'])}; E3b (strain not in database) k-mer "
                     f"{common.fmt(f['kmer_jaccard']['auroc'])} [{common.fmt(f['kmer_jaccard']['ci95'][0])}, "
                     f"{common.fmt(f['kmer_jaccard']['ci95'][1])}] vs leaf {common.fmt(f['leaf_abundance_sim']['auroc'])} "
                     f"({f['kmer_jaccard']['n_pos']}+/{f['kmer_jaccard']['n_neg']}- pairs)")
    for cond in ("E3a", "E3b"):
        for m in ("kmer_jaccard", "leaf_abundance_sim"):
            lim = next((c for c in COVERAGES if res[cond][c][m]["auroc"] >= 0.95), None)
            res[cond][f"detection_limit_{m}"] = lim
            lines.append(f"- {cond} detection limit ({m}, AUROC >= 0.95): {lim if lim else '> ' + str(COVERAGES[-1])}x")
    b1 = res["E3b"][1.0]["kmer_jaccard"]["auroc"]
    lines.append(f"- verdict: E3b at 1x = {common.fmt(b1)}: " + (
        "k-mer sharing tracks strains absent from the database" if b1 > 0.9 else
        "k-mer sharing needs the strain in the database (state as a limitation)" if b1 < 0.7 else
        "partial: strains absent from the database are tracked only weakly"))
    common.write_results("m3", res, lines)
