"""M1: profile sufficiency, streaming and cost (see apps/README.md).

E1a merge(profiles of K random read shards) == profile(all reads), byte for byte
E1b `profile report` == `classify -U`, byte for byte (HyperLogLog and MinHash)
E1c profile of the file read directly == profile of the same file through a pipe
E1d cost: reads/s, Mbp/s, peak RSS, profile bytes / input bytes (pass: <= 1%)"""
import filecmp
import gzip
import os
import random

from . import common, mock, simreads, stream

HIFI = os.path.join(common.REPO, "test", "zymo", "SRR13128014.head150M.fastq.gz")
SEED = 20261020


def fastq_records(path):
    """Complete FASTQ records of a (possibly truncated) gzip file."""
    out = []
    try:
        with gzip.open(path, "rt") as fh:
            while True:
                rec = [fh.readline() for _ in range(4)]
                if not rec[3].endswith("\n"):
                    break
                out.append("".join(rec))
    except EOFError:
        pass   # truncated gzip member: keep the complete records read so far
    return out


def records_source(recs):
    return simreads.chunks(iter(recs))


def run():
    db, nodes, dbinfo = mock.build_db("full", mock.genomes("full"))
    recs = fastq_records(HIFI)
    nbp = sum(len(r.split("\n")[1]) for r in recs)
    res = {"seed": SEED, "db": dbinfo, "hifi_reads": len(recs), "hifi_bp": nbp, "E1a": {}, "E1b": {}}
    lines = ["# M1: sufficiency, streaming, cost", f"- data: {len(recs)} HiFi reads ({nbp / 1e6:.1f} Mbp), "
             f"head of SRR13128014; database full ({dbinfo['genomes']} genomes, {dbinfo['bytes'] / 1e6:.0f} MB)"]
    ok_all = True
    # E1a: shards (MinHash profiles; the merge must reproduce the whole sample's name)
    whole = common.work("m1", "whole", "hifi.bprof")
    cost = stream.profile_stream(records_source(recs), whole, db, nodes)
    for K in (2, 4, 8):
        rng = random.Random(SEED + K)
        assign = [rng.randrange(K) for _ in recs]
        shards = []
        for s in range(K):
            sp = common.work("m1", f"K{K}", f"shard{s}.bprof")
            stream.profile_stream(records_source([r for r, a in zip(recs, assign) if a == s]), sp, db, nodes)
            shards.append(sp)
        merged = common.work("m1", f"K{K}", "merged", "hifi.bprof")
        common.bonsai("profile", "merge", merged, *shards)
        same = filecmp.cmp(merged, whole, shallow=False)
        res["E1a"][K] = same
        ok_all &= same
    lines.append(f"- E1a merge(K shards) == whole: {res['E1a']} -> {'PASS' if all(res['E1a'].values()) else 'FAIL'}")
    # E1b: report == -U, both counter kinds
    for kind, mh in (("hll", False), ("minhash", True)):
        prof = common.work("m1", f"rep_{kind}", "hifi.bprof")
        u = common.work("m1", f"rep_{kind}", "classify_U.tsv")
        stream.profile_stream(records_source(recs), prof, db, nodes, minhash=mh, uniq=u)
        rep = common.work("m1", f"rep_{kind}", "profile_report.tsv")
        common.bonsai("profile", "report", prof, rep)
        same = filecmp.cmp(u, rep, shallow=False)
        res["E1b"][kind] = same
        ok_all &= same
    lines.append(f"- E1b profile report == classify -U: {res['E1b']} -> {'PASS' if all(res['E1b'].values()) else 'FAIL'}")
    # E1c: direct file vs pipe of the same file
    direct = common.work("m1", "direct", "hifi_file.bprof")
    piped = common.work("m1", "piped", "hifi_file.bprof")
    common.bonsai(*stream.classify_args(direct, db, nodes)[:-1], HIFI)
    stream.profile_stream(stream.file_source(HIFI), piped, db, nodes)
    res["E1c"] = filecmp.cmp(direct, piped, shallow=False)
    ok_all &= res["E1c"]
    lines.append(f"- E1c file == pipe: {res['E1c']} -> {'PASS' if res['E1c'] else 'FAIL'}")
    # E1d: cost, real HiFi and simulated 150 bp reads (1% error, mock genomes uniform)
    gz_bytes = os.path.getsize(HIFI)
    hifi_cost = dict(cost, input_gz_bytes=gz_bytes, reads_per_s=len(recs) / cost["seconds"],
                     mbp_per_s=nbp / 1e6 / cost["seconds"], profile_frac_of_gz=cost["profile_bytes"] / gz_bytes)
    comp = [(mock.mock_path(f), t, 20000) for f, t in mock.MOCK]
    sim = common.work("m1", "sim", "sim150.bprof")
    sc = stream.profile_stream(simreads.chunks(simreads.reads(comp, 150, 0.01, SEED)), sim, db, nodes)
    nsim = sum(n for _, _, n in comp)
    sim_cost = dict(sc, reads=nsim, reads_per_s=nsim / sc["seconds"], mbp_per_s=nsim * 150 / 1e6 / sc["seconds"],
                    profile_frac_of_plain=sc["profile_bytes"] / sc["bytes"])
    res["E1d"] = {"hifi": hifi_cost, "sim150": sim_cost}
    frac_ok = hifi_cost["profile_frac_of_gz"] <= 0.01
    ok_all &= frac_ok
    lines.append(f"- E1d HiFi: {hifi_cost['reads_per_s']:.0f} reads/s, {hifi_cost['mbp_per_s']:.1f} Mbp/s, "
                 f"peak RSS {hifi_cost['peak_rss_bytes'] / 1e6:.0f} MB, profile {hifi_cost['profile_bytes'] / 1e3:.0f} kB = "
                 f"{100 * hifi_cost['profile_frac_of_gz']:.3f}% of the gzipped input -> {'PASS' if frac_ok else 'FAIL'}")
    lines.append(f"- E1d simulated 150 bp: {sim_cost['reads_per_s']:.0f} reads/s ({sim_cost['mbp_per_s']:.1f} Mbp/s, "
                 f"includes the Python simulator), peak RSS {sim_cost['peak_rss_bytes'] / 1e6:.0f} MB, profile "
                 f"{sim_cost['profile_bytes'] / 1e3:.0f} kB = {100 * sim_cost['profile_frac_of_plain']:.3f}% of the plain FASTQ")
    lines.append(f"- M1 overall: {'PASS' if ok_all else 'FAIL'}")
    common.write_results("m1", res, lines)
