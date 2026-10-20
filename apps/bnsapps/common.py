"""Shared helpers for the profile applications: paths, running bonsai, TSV
parsing, and the statistics the endpoints use (stdlib only)."""
import json
import math
import os
import random
import subprocess
import time

APPS = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
REPO = os.path.dirname(APPS)
BONSAI = os.environ.get("BONSAI", os.path.join(REPO, "bonsai"))
RESULTS = os.path.join(APPS, "results")
# Work files (databases, profiles) live outside the repository.
WORK = os.environ.get("BNS_APPS_WORK", os.path.join(os.environ.get("TMPDIR", "/tmp"), "bns_apps"))


def work(*parts):
    p = os.path.join(WORK, *parts)
    os.makedirs(os.path.dirname(p), exist_ok=True)
    return p


def bonsai(*args, stdout=subprocess.DEVNULL, check=True, **kw):
    """Run bonsai; stderr captured and shown on failure."""
    r = subprocess.run([BONSAI, *map(str, args)], stdout=stdout, stderr=subprocess.PIPE, **kw)
    if check and r.returncode != 0:
        raise RuntimeError(f"bonsai {' '.join(map(str, args))} failed ({r.returncode}):\n" + r.stderr.decode()[-2000:])
    return r


def bonsai_tsv(*args):
    """Run a bonsai command that writes a TSV to stdout; returns (comment lines, header, rows as dicts)."""
    out = bonsai(*args, stdout=subprocess.PIPE).stdout.decode()
    return parse_tsv(out)


def parse_tsv(text):
    comments, header, rows = [], None, []
    for line in text.splitlines():
        if not line:
            continue
        if line.startswith("#"):
            if header is None and "\t" in line and not line.startswith(("#groups", "#permanova")):
                header = line[1:].split("\t")
            else:
                comments.append(line)
            continue
        f = line.split("\t")
        rows.append(dict(zip(header, f)) if header else f)
    return comments, header, rows


def read_tsv(path):
    with open(path) as fh:
        return parse_tsv(fh.read())


# ---- statistics ----

def auroc(pos, neg):
    """P(score_pos > score_neg) + 0.5 P(tie): Mann-Whitney AUROC."""
    if not pos or not neg:
        return float("nan")
    allv = sorted([(v, 1) for v in pos] + [(v, 0) for v in neg])
    ranks, i = {}, 0
    rank_sum = 0.0
    while i < len(allv):
        j = i
        while j < len(allv) and allv[j][0] == allv[i][0]:
            j += 1
        r = (i + j + 1) / 2.0  # average rank, 1-based
        rank_sum += r * sum(1 for k in range(i, j) if allv[k][1] == 1)
        i = j
    n1, n0 = len(pos), len(neg)
    return (rank_sum - n1 * (n1 + 1) / 2.0) / (n1 * n0)


def bootstrap_ci(stat, groups, n=2000, seed=1, level=0.95):
    """Percentile CI of stat(*resampled groups); groups resampled independently."""
    rng = random.Random(seed)
    vals = []
    for _ in range(n):
        res = [[g[rng.randrange(len(g))] for _ in g] for g in groups]
        v = stat(*res)
        if not math.isnan(v):
            vals.append(v)
    vals.sort()
    if not vals:
        return (float("nan"), float("nan"))
    lo = vals[int((1 - level) / 2 * len(vals))]
    hi = vals[min(len(vals) - 1, int((1 + level) / 2 * len(vals)))]
    return (lo, hi)


def ks_uniform(ps):
    """One-sample KS statistic against U(0,1) and its asymptotic p-value."""
    ps = sorted(p for p in ps if not math.isnan(p))
    n = len(ps)
    if n == 0:
        return float("nan"), float("nan")
    d = max(max((i + 1) / n - p, p - i / n) for i, p in enumerate(ps))
    # Kolmogorov distribution tail (Marsaglia-style series)
    x = (math.sqrt(n) + 0.12 + 0.11 / math.sqrt(n)) * d
    p = 2 * sum((-1) ** (k - 1) * math.exp(-2 * k * k * x * x) for k in range(1, 101))
    return d, max(0.0, min(1.0, p))


def ranksum_p(a, b):
    """Two-sided Wilcoxon rank-sum p (normal approximation, tie-corrected)."""
    n1, n2 = len(a), len(b)
    if n1 == 0 or n2 == 0:
        return float("nan")
    allv = sorted([(v, 0) for v in a] + [(v, 1) for v in b])
    rank_a, i, ties = 0.0, 0, 0.0
    while i < len(allv):
        j = i
        while j < len(allv) and allv[j][0] == allv[i][0]:
            j += 1
        t = j - i
        ties += t ** 3 - t
        r = (i + j + 1) / 2.0
        rank_a += r * sum(1 for k in range(i, j) if allv[k][1] == 0)
        i = j
    n = n1 + n2
    mu = n1 * (n + 1) / 2.0
    var = n1 * n2 / 12.0 * ((n + 1) - ties / (n * (n - 1)))
    if var <= 0:
        return 1.0
    z = (abs(rank_a - mu) - 0.5) / math.sqrt(var)
    return max(0.0, min(1.0, math.erfc(max(z, 0.0) / math.sqrt(2))))


def bh(ps):
    """Benjamini-Hochberg q-values; NaN stays NaN."""
    idx = sorted((i for i, p in enumerate(ps) if not math.isnan(p)), key=lambda i: ps[i])
    q = [float("nan")] * len(ps)
    run = 1.0
    m = len(idx)
    for r in range(m - 1, -1, -1):
        i = idx[r]
        run = min(run, ps[i] * m / (r + 1))
        q[i] = run
    return q


def fmt(x, d=3):
    return "NA" if x is None or (isinstance(x, float) and math.isnan(x)) else f"{x:.{d}f}"


# ---- results ----

def write_results(aid, data, verdict_lines):
    os.makedirs(RESULTS, exist_ok=True)
    data = dict(data, analysis=aid, written=time.strftime("%Y-%m-%d %H:%M:%S"))
    with open(os.path.join(RESULTS, aid + ".json"), "w") as fh:
        json.dump(data, fh, indent=1, sort_keys=True, default=str)
    with open(os.path.join(RESULTS, aid + ".md"), "w") as fh:
        fh.write("\n".join(verdict_lines) + "\n")
    print("\n".join(verdict_lines))
