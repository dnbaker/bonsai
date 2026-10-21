"""Spaced-seed design (PI 2026-10-20: "how to choose spacing").

Sensitivity of a seed pattern = how often a homologous read of given divergence
gets hits. Estimated by Monte Carlo over match/mismatch strings of a 150 bp
read: P(>= 1 hit), P(>= 3 hits) and expected hits. Models:
  uniform  i.i.d. substitutions at rate d;
  codon    coding fraction CODING with codon-position rates proportional to
           CODON_RATIO (third positions ~4x the first: synonymous changes),
           reading frame random; non-coding positions at rate d. These rates are
           an ASSUMPTION (no aligner on this host to fit them); real held-out
           divergence decides between patterns.
Optimization: hill climbing over gap placements at fixed weight and span <=
SPAN_MAX (one 64-bit DNA window), maximizing P(>= 3 hits) at D_TARGET."""
import random

from . import common

L = 150
CODING = 0.88
CODON_RATIO = (1.0, 0.6, 4.0)
SPAN_MAX = 32
D_TARGET = 0.03
N_MC = 4000
SEED = 20261030


def pattern_str(offsets):
    span = offsets[-1] + 1
    s = ["0"] * span
    for o in offsets:
        s[o] = "1"
    return "".join(s)


def offsets_of(pattern):
    return [i for i, c in enumerate(pattern) if c == "1"]


def mismatch_string(rng, d, model):
    """True = mismatch, for one read."""
    if model == "uniform":
        return [rng.random() < d for _ in range(L)]
    w = sum(CODON_RATIO) / 3.0
    rates = [d * r / w for r in CODON_RATIO]
    out, i = [], 0
    while i < L:
        seg = rng.randint(60, 900)   # coding run (gene piece) or non-coding run
        coding = rng.random() < CODING
        phase = rng.randrange(3)
        for j in range(min(seg, L - i)):
            out.append(rng.random() < (rates[(phase + j) % 3] if coding else d))
        i += seg
    return out[:L]


def hits(offsets, mm):
    span = offsets[-1] + 1
    return sum(1 for s in range(L - span + 1) if not any(mm[s + o] for o in offsets))


def sensitivity(offsets, d, model, n=N_MC, seed=SEED):
    rng = random.Random(seed)   # common random numbers across patterns
    h = [hits(offsets, mismatch_string(rng, d, model)) for _ in range(n)]
    return {"p_ge1": sum(1 for x in h if x >= 1) / n, "p_ge3": sum(1 for x in h if x >= 3) / n, "mean_hits": sum(h) / n}


def contiguous(weight):
    return list(range(weight))


def optimize(weight, model, d=D_TARGET, iters=300, seed=SEED):
    """Hill climb: move one interior position at a time, keep the span <= SPAN_MAX,
    objective P(>= 3 hits) at d under the model (smaller n during the search)."""
    rng = random.Random(seed)
    span = min(SPAN_MAX, weight + max(4, (SPAN_MAX - weight)))
    # start: weight positions spread over the span, ends fixed
    cur = sorted({0, span - 1} | set(rng.sample(range(1, span - 1), weight - 2)))
    score = lambda o: sensitivity(o, d, model, n=800)["p_ge3"]
    best = score(cur)
    for _ in range(iters):
        cand = cur[:]
        i = rng.randrange(1, weight - 1)
        free = [p for p in range(1, span - 1) if p not in cand]
        cand[i] = rng.choice(free)
        cand.sort()
        sc = score(cand)
        if sc >= best:
            cur, best = cand, sc
    return cur


def run():
    res, lines = {"model": {"coding": CODING, "codon_ratio": CODON_RATIO, "d_target": D_TARGET, "span_max": SPAN_MAX}, "patterns": {}}, [
        "# Spaced-seed design (Monte Carlo sensitivity; codon rates are an assumption)",
        "| weight | pattern | name | P(>=3 hits) uniform 2% / 5% | P(>=3 hits) codon 2% / 5% | mean hits codon 3% |",
        "|---|---|---|---|---|---|"]
    for weight in (21, 25):
        pats = {"contiguous": contiguous(weight), "uniform_opt": optimize(weight, "uniform"), "codon_opt": optimize(weight, "codon")}
        for name, o in pats.items():
            su = [sensitivity(o, d, "uniform")["p_ge3"] for d in (0.02, 0.05)]
            sc = [sensitivity(o, d, "codon")["p_ge3"] for d in (0.02, 0.05)]
            mh = sensitivity(o, 0.03, "codon")["mean_hits"]
            res["patterns"][f"w{weight}_{name}"] = {"pattern": pattern_str(o), "uniform": su, "codon": sc, "mean_hits_codon3": mh}
            lines.append(f"| {weight} | `{pattern_str(o)}` | {name} | {su[0]:.3f} / {su[1]:.3f} | {sc[0]:.3f} / {sc[1]:.3f} | {mh:.1f} |")
    lines.append("- Next: validate on real held-out divergence (podar strain pairs, mock E. coli strains), after the spaced-seed "
                 "strand check passes.")
    common.write_results("seeds", res, lines)
