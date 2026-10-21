"""Read simulation for mock analyses, streamed (FASTQ bytes into a pipe; no
files). Reads are drawn uniformly from either strand of a source genome, with
Binomial(L, error) substitutions; read names carry the source taxid
(@t<taxid>_<n>) so per-read truth stays available."""
import math
import random

from . import mock

_COMP = str.maketrans("ACGTN", "TGCAN")
_GENOMES = {}


def genome(path):
    if path not in _GENOMES:
        _GENOMES[path] = mock.read_fasta(path)
    return _GENOMES[path]


def revcomp(s):
    return s.translate(_COMP)[::-1]


def reads(composition, length=150, error=0.01, seed=1, strand=None):
    """composition: [(genome path, taxid, n reads)]. Yields FASTQ records (str).
    strand: None (random), '+' (forward only) or '-' (reverse complement only)."""
    rng = random.Random(seed)
    qual = "I" * length
    log1m = math.log1p(-error) if 0 < error < 1 else 0.0
    alt = {"A": "CGT", "C": "AGT", "G": "ACT", "T": "ACG"}
    for path, taxid, n in composition:
        g = genome(path)
        span = len(g) - length
        for i in range(n):
            o = rng.randrange(span)
            s = g[o:o + length]
            flip = rng.random() < 0.5
            if strand == "-" or (strand is None and flip):
                s = revcomp(s)
            # independent substitutions at rate `error` (Binomial(L, error) of them):
            # geometric gaps between error positions, ~1 draw per error, not per base
            if error > 0:
                pos = int(math.log(1.0 - rng.random()) / log1m) if log1m else length
                if pos < length:
                    s = list(s)
                    while pos < length:
                        c = s[pos]
                        s[pos] = alt[c][rng.randrange(3)] if c in alt else "A"
                        pos += 1 + int(math.log(1.0 - rng.random()) / log1m)
                    s = "".join(s)
            yield f"@t{taxid}_{i}\n{s}\n+\n{qual}\n"


def chunks(records, size=1 << 20):
    """Group FASTQ text into ~size-byte bytes chunks for a pipe."""
    buf, n = [], 0
    for r in records:
        buf.append(r)
        n += len(r)
        if n >= size:
            yield "".join(buf).encode()
            buf, n = [], 0
    if buf:
        yield "".join(buf).encode()


def dirichlet(rng, alphas):
    xs = [rng.gammavariate(a, 1.0) if a > 0 else 0.0 for a in alphas]
    s = sum(xs)
    return [x / s for x in xs]


def multinomial(rng, n, ps):
    """Counts of n draws over probabilities ps (sequential binomials via normal or exact)."""
    counts, left, rest = [], n, 1.0
    for i, p in enumerate(ps):
        if i == len(ps) - 1:
            counts.append(left)
            break
        q = min(1.0, p / rest) if rest > 0 else 0.0
        c = binomial(rng, left, q)
        counts.append(c)
        left -= c
        rest -= p
    return counts


def binomial(rng, n, p):
    if n <= 0 or p <= 0:
        return 0
    if p >= 1:
        return n
    if n < 1000:
        return sum(1 for _ in range(n) if rng.random() < p)
    mu, sd = n * p, (n * p * (1 - p)) ** 0.5
    return max(0, min(n, int(round(rng.gauss(mu, sd)))))
