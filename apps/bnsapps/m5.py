"""M5: nearest-profile search (see apps/README.md). Synthetic, so the result
depends on the design; two regimes are reported rather than one tuned to win.

Subjects: species proportions p_s ~ Dirichlet(BETA p_pop) (BETA large = subjects
alike in abundance; small = individual abundance profiles) and 1-2 E. coli
strains. Samples: 3 per subject, q = 0.8 p_s + 0.2 Dirichlet(1) (composition
drift), p_i ~ Dirichlet(50 q), reads ~ Multinomial(N, p_i). All samples are
indexed; each queries the index (every other sample ranked).
E5: top-1 same-subject retrieval (self excluded) and AUROC of same- vs
different-subject pairs, for tree EMD, Bray-Curtis, k-mer Jaccard, containment."""
import random

from . import common, m4, mock, simreads, stream

SEED = 20261024
N_SUBJECTS, PER_SUBJECT, N_READS = 20, 3, 100000
BETAS = (5.0, 100.0)
METRICS = {"tree_emd": -1, "bray_curtis": -1, "kmer_jaccard": 1, "cont_query_in_sample": 1}   # sign: higher = closer


def run():
    rng = random.Random(SEED)
    db, nodes, _ = mock.build_db("mock", mock.genomes("mock"))
    strains = list(mock.STRAINS)
    res = {"seed": SEED, "params": {"subjects": N_SUBJECTS, "per_subject": PER_SUBJECT, "reads": N_READS, "betas": BETAS}}
    lines = [f"# M5: nearest-profile search ({N_SUBJECTS} subjects x {PER_SUBJECT} samples, {N_READS} reads each)"]
    for beta in BETAS:
        p_pop = m4.base_props(rng)
        subj, profs = {}, []
        for s in range(N_SUBJECTS):
            ps = simreads.dirichlet(rng, [beta * x for x in p_pop])
            st = rng.sample(strains, rng.choice((1, 2)))
            for j in range(PER_SUBJECT):
                noise = simreads.dirichlet(rng, [1.0] * len(ps))
                q = [0.8 * a + 0.2 * b for a, b in zip(ps, noise)]
                pi = simreads.dirichlet(rng, [50.0 * x for x in q])
                counts = simreads.multinomial(rng, N_READS, pi)
                comp = [(mock.mock_path(f), t, c) for (f, t), c in zip(m4.SPECIES, counts[:-1])]
                for k, t in enumerate(st):   # E. coli share split over the subject's strains
                    comp.append((mock.mock_path(m4.ECOLI_FILE[t]), t, counts[-1] // len(st) + (k < counts[-1] % len(st))))
                name = f"b{int(beta)}_s{s}_{j}"
                prof = common.work("m5", f"b{int(beta)}", f"{name}.bprof")
                stream.profile_stream(simreads.chunks(simreads.reads(comp, 150, 0.01, rng.randrange(1 << 30))), prof, db, nodes)
                subj[name] = s
                profs.append(prof)
        idx = common.work("m5", f"b{int(beta)}", "index.bidx")
        common.bonsai("profile", "index", idx, *profs)
        _, _, rows = common.bonsai_tsv("profile", "search", "-n", len(profs), idx, *profs)
        hits = {}
        for r in rows:
            if r["query"] != r["sample"]:
                hits.setdefault(r["query"], []).append(r)
        e = {}
        for m, sign in METRICS.items():
            top1 = sum(1 for q, hs in hits.items() if subj[max(hs, key=lambda h: sign * float(h[m]))["sample"]] == subj[q])
            seen, pos, neg = set(), [], []
            for q, hs in hits.items():
                for h in hs:
                    key = frozenset((q, h["sample"]))
                    if key in seen:
                        continue
                    seen.add(key)
                    (pos if subj[q] == subj[h["sample"]] else neg).append(sign * float(h[m]))
            e[m] = {"top1": top1 / len(hits), "auroc": common.auroc(pos, neg),
                    "ci95": common.bootstrap_ci(common.auroc, [pos, neg], n=500, seed=SEED)}
        res[f"beta{beta}"] = e
        lines.append(f"- beta {beta} ({'individual abundance' if beta < 20 else 'subjects alike in abundance'}): " + "; ".join(
            f"{m} top-1 {common.fmt(v['top1'], 2)}, AUROC {common.fmt(v['auroc'])}" for m, v in e.items()))
    lines.append("- note: synthetic, and the regimes are chosen, not estimated; R1 (real longitudinal data) decides "
                 "whether k-mer or abundance metrics identify subjects better.")
    common.write_results("m5", res, lines)
