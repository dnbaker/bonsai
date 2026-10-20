"""M4 calibration follow-up: M4's first 20 null experiments gave a k-mer-test
FPR of 0.080 +/- 0.011 (experiment-level SE), above alpha = 0.05, although an
exact permutation test on exchangeable samples cannot exceed alpha in
expectation. This runs R more null experiments (own seed and directories, so
M4 stays reproducible) to tell chance from a defect."""
import math
import random

from . import common, m4, mock

SEED = 20261025
R = 60


def run():
    rng = random.Random(SEED)
    db, nodes, _ = mock.build_db("mock", mock.genomes("mock"))
    p = m4.base_props(rng)
    null = [m4.experiment(f"cal{r}", rng, p, db, nodes) for r in range(R)]
    res = {"seed": SEED, "R": R}
    lines = [f"# M4cal: {R} more null experiments (6 vs 6)"]
    for key, name in (("p", "reads"), ("kp", "kmers"), ("wp", "wilcoxon")):
        rates = []
        for ex in null:
            ps = [v[key] for v in ex.values() if not math.isnan(v[key])]
            if ps:
                rates.append(sum(1 for x in ps if x < 0.05) / len(ps))
        m = sum(rates) / len(rates)
        se = math.sqrt(sum((x - m) ** 2 for x in rates) / (len(rates) - 1)) / math.sqrt(len(rates))
        allp = [v[key] for ex in null for v in ex.values()]
        d, ksp = common.ks_uniform(allp)
        res[name] = {"fpr": m, "se": se, "ks_d": d, "ks_p": ksp}
        lines.append(f"- {name}: FPR {common.fmt(m)} +/- {common.fmt(se)} (experiment-level SE); "
                     f"z vs 0.05 = {common.fmt((m - 0.05) / se if se else float('nan'), 2)}; KS D {common.fmt(d)} p {common.fmt(ksp)}")
    res["null_pvalues"] = [{str(t): [v["p"], v["kp"], v["wp"]] for t, v in ex.items()} for ex in null]
    common.write_results("m4cal", res, lines)
