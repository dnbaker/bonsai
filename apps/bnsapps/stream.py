"""Profile a stream of reads: a byte source (simulator chunks, a local file
through `cat`, or `curl` of an ENA URL) piped into `bonsai classify -P`, with
per-read output discarded. Only the profile reaches disk. Returns cost stats;
bytes are counted as they pass, for the download tally."""
import json
import os
import subprocess
import tempfile
import time
import urllib.parse
import urllib.request

from . import common

CURL = "/usr/bin/curl"   # the miniconda curl on PATH lacks libssl.1.1


def classify_args(profile, db, nodes, minhash=True, uniq=None, threads=2, extra=()):
    a = ["classify", "-k", "-p", threads]
    if minhash:
        a.append("-m")
    if uniq:
        a += ["-U", uniq]
    a += ["-P", profile, *extra, db, nodes, "/dev/stdin"]
    return a


def profile_stream(source, profile, db, nodes, **kw):
    """source: iterable of bytes chunks, or a command list whose stdout is the reads.
    Returns {seconds, bytes, peak_rss_bytes, profile_bytes}."""
    # /usr/bin/time -l (macOS) reports this run's peak RSS on stderr
    cmd = ["/usr/bin/time", "-l", common.BONSAI, *map(str, classify_args(profile, db, nodes, **kw))]
    t0 = time.time()
    # stderr to a file: classify logs per chunk, and a full stderr pipe would
    # block it while we block writing its stdin
    errf = tempfile.TemporaryFile()
    p = subprocess.Popen(cmd, stdin=subprocess.PIPE, stdout=subprocess.DEVNULL, stderr=errf)
    nbytes = 0
    producer = None
    try:
        if isinstance(source, (list, tuple)) and source and isinstance(source[0], str):
            producer = subprocess.Popen(source, stdout=subprocess.PIPE)
            it = iter(lambda: producer.stdout.read(1 << 20), b"")
        else:
            it = source
        for chunk in it:
            nbytes += len(chunk)
            p.stdin.write(chunk)
        p.stdin.close()
    except BrokenPipeError:
        pass
    rc = p.wait()
    errf.seek(0)
    err = errf.read()
    errf.close()
    if producer is not None:
        producer.wait()
    if rc != 0:
        raise RuntimeError(f"classify failed ({rc}): {err.decode()[-2000:]}")
    rss = None
    for line in err.decode(errors="replace").splitlines():
        if "maximum resident set size" in line:
            rss = int(line.split()[0])
    return {"seconds": round(time.time() - t0, 3), "bytes": nbytes, "profile_bytes": os.path.getsize(profile),
            "peak_rss_bytes": rss}


def file_source(path):
    return ["cat", path]


def url_source(url, max_bytes=None):
    """curl a URL; max_bytes stops the stream after that many bytes (the head of a
    run; bonsai stops cleanly at a truncated gzip member)."""
    cmd = [CURL, "-sSfL", url]
    if max_bytes:
        cmd = ["sh", "-c", f"{CURL} -sSfL '{url}' | head -c {int(max_bytes)}"]
    return cmd


def ena_runs(accession):
    """ENA read-run metadata (fastq URLs and sizes) for a run/sample/study accession.
    Metadata only: no reads are fetched."""
    q = urllib.parse.urlencode({"accession": accession, "result": "read_run", "format": "json",
                                "fields": "run_accession,fastq_ftp,fastq_bytes,read_count,instrument_platform,library_layout"})
    with urllib.request.urlopen("https://www.ebi.ac.uk/ena/portal/api/filereport?" + q, timeout=60) as r:
        return json.load(r)
