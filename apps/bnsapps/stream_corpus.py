"""Stream a selected corpus (corpus.py manifests) into profiles against the gut
database (gutdb.py): per sample, the R1 head (<= CORPUS_MAX_BYTES) from ENA into
classify -m -P; resume (existing profiles are skipped); memory gate before each
run (the night gate: BNS_APPS_SWAP_MAX_GB=2, BNS_APPS_FREE_MIN_GB=4); failures
are logged and skipped. Light: one classify at a time, network-bound.

  CORPUS_SET=normals|crc  CORPUS_LIMIT=<n>  run.py stream_corpus"""
import csv
import json
import os
import time

from . import common, gutdb, stream

MAX_BYTES = int(float(os.environ.get("CORPUS_MAX_BYTES", 60e6)))


def run():
    which = os.environ.get("CORPUS_SET", "normals")
    limit = int(os.environ.get("CORPUS_LIMIT", 0)) or None
    manifest = os.path.join(common.RESULTS, f"corpus_{which}.tsv")
    with open(manifest, newline="") as fh:
        rows = list(csv.DictReader(fh, delimiter="\t"))[:limit]
    sp = gutdb.species_table(gutdb.metadata())[:gutdb.N_SPECIES]
    db, nodes, _ = gutdb.build("full", sp)
    log_path = common.work("corpus", which, "stream_log.jsonl")
    done = streamed = failed = 0
    for r in rows:
        name = f"{r['study']}__{r['sample_id']}".replace("/", "_")
        prof = common.work("corpus", which, "profiles", f"{name}.bprof")
        if common.reuse(prof):
            done += 1
            continue
        try:
            cost = stream.profile_stream(stream.url_source(r["r1_url"], MAX_BYTES), prof, db, nodes, threads=2)
        except Exception as e:   # one bad run must not stop the night
            failed += 1
            if os.path.exists(prof):
                os.remove(prof)
            cost = {"error": str(e)[-300:]}
        else:
            streamed += cost["bytes"]
            done += 1
        with open(log_path, "a") as fh:
            fh.write(json.dumps(dict(sample=name, study=r["study"], condition=r["condition"], time=time.strftime("%H:%M:%S"), **cost)) + "\n")
    lines = [f"# Corpus stream ({which})", f"- {done}/{len(rows)} profiles; {failed} failed; streamed this run "
             f"{streamed / 1e9:.2f} GB (heads <= {MAX_BYTES / 1e6:.0f} MB); log {log_path}"]
    common.write_results(f"stream_{which}", {"done": done, "failed": failed, "bytes": streamed, "n": len(rows)}, lines)
