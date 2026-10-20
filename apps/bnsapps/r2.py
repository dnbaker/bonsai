"""R2: the sourmash gather paper's samples (podar SRR606249, Zymo D6300
SRR12324253, Alaska oil reservoir SRR1976948), profiled with bonsai against
each sample's gather-matched genomes (gatherdata.py), R1 streamed from ENA
(head of R2_MAX_BYTES). Reference abundances: the paper's read mapping
(coverage x genome length, which is proportional to reads) and gather's
f_unique_weighted. Circular for presence (database = gather matches): this
evaluates abundance estimation and per-genome scores, not database search.

E10a abundance: Spearman and log10 Pearson of per-genome abundance against
     mapping, for bonsai (genome-leaf reads; and species-level sums) and gather.
E10b genome-level resolution: fraction of classified reads that land on a
     genome leaf, vs stopping at species or above (near-identical genomes, e.g.
     podar's strain sets, push reads up the tree).
E10c gather false positives: genomes gather reports whose mapping covers < 50%
     of the genome (likely spurious matches). AUROC of bonsai per-genome coverage
     (cov) and read count at separating them from well-covered genomes."""
import math
import os

from . import common, gatherdata, r0, stream

MAX_BYTES = int(float(os.environ.get("R2_MAX_BYTES", 300e6)))


def run():
    res, lines = {"max_bytes": MAX_BYTES}, [f"# R2: sourmash gather paper samples (R1 head <= {MAX_BYTES / 1e6:.0f} MB, gather-matched genomes as database)"]
    only = os.environ.get("R2_SAMPLES")
    for sample, (run_acc, label) in gatherdata.SAMPLES.items():
        if only and sample not in only.split(","):
            continue
        d = gatherdata.build(sample)
        meta = r0.ena_search(f"run_accession={run_acc}")[0]
        prof = common.work("gather", sample, f"{sample}.bprof")
        rep = common.work("gather", sample, f"{sample}.report.tsv")
        cost = None
        if not (common.reuse(prof) and common.reuse(rep)):
            cost = stream.profile_stream(stream.url_source(meta["r1"], MAX_BYTES), prof, d["db"], d["nodes"], uniq=rep, threads=4)
        _, _, rows = common.read_tsv(rep)
        by = {int(r["taxid"]): r for r in rows if r["taxid"] != "0"}
        classified = sum(int(r["taxon_reads"]) for r in rows if r["taxid"] != "0")
        mapping = gatherdata.mapping_rows(sample)
        recs = []
        for g in d["rows"]:
            acc = gatherdata.accession(g)
            m = mapping.get(acc.split(".")[0])
            if m is None:
                continue
            leaf = by.get(d["genome_id"][acc])
            sp = by.get(d["species_id"][acc])
            cov_bonsai = float(leaf["cov"]) if leaf and leaf["cov"] != "NA" else 0.0
            recs.append({"acc": acc, "map_abund": float(m["coverage"]) * float(m["genome bp"]),
                         "map_covered_frac": 1 - float(m["percent missed"]) / 100.0,
                         "gather": float(g["f_unique_weighted"]),
                         "bonsai_leaf_reads": float(leaf["clade_reads"]) if leaf else 0.0,
                         "bonsai_species_reads": float(sp["clade_reads"]) if sp else 0.0, "bonsai_cov": cov_bonsai,
                         "species": d["species_id"][acc]})
        lg = lambda v: [math.log10(x + 1e-9) for x in v]
        ma = [r["map_abund"] for r in recs]
        e = {}
        for key in ("bonsai_leaf_reads", "gather"):
            v = [r[key] for r in recs]
            e[key] = {"spearman": common.spearman(v, ma), "pearson_log10": common.pearson(lg(v), lg(ma))}
        # species level: sum mapping and gather over genomes of a species; bonsai species clade reads
        sp_map, sp_g, sp_b = {}, {}, {}
        for r in recs:
            sp_map[r["species"]] = sp_map.get(r["species"], 0) + r["map_abund"]
            sp_g[r["species"]] = sp_g.get(r["species"], 0) + r["gather"]
            sp_b[r["species"]] = r["bonsai_species_reads"]
        keys = list(sp_map)
        e["species_bonsai"] = {"spearman": common.spearman([sp_b[k] for k in keys], [sp_map[k] for k in keys]),
                               "pearson_log10": common.pearson(lg([sp_b[k] for k in keys]), lg([sp_map[k] for k in keys]))}
        e["species_gather"] = {"spearman": common.spearman([sp_g[k] for k in keys], [sp_map[k] for k in keys]),
                               "pearson_log10": common.pearson(lg([sp_g[k] for k in keys]), lg([sp_map[k] for k in keys]))}
        leaf_reads = sum(r["bonsai_leaf_reads"] for r in recs)
        fp = [r for r in recs if r["map_covered_frac"] < 0.5]
        ok = [r for r in recs if r["map_covered_frac"] >= 0.5]
        e10c = {s: common.auroc([r[s] for r in ok], [r[s] for r in fp]) for s in ("bonsai_cov", "bonsai_leaf_reads", "gather")} if fp and ok else {}
        res[sample] = {"label": label, "genomes": len(recs), "classified": classified, "leaf_frac": leaf_reads / classified if classified else float("nan"),
                       "E10a": e, "E10c": e10c, "n_gather_fp": len(fp), "cost": cost, "downloaded_bytes": d["downloaded_bytes"],
                       "missing_genomes": d["missing"], "records": recs}
        lines.append(f"- {label} ({sample}, {len(recs)} genomes{', %d unavailable' % len(d['missing']) if d['missing'] else ''}, {classified} classified reads"
                     f"{'' if not cost else ', streamed %.0f MB' % (cost['bytes'] / 1e6)}):")
        lines.append(f"  - E10a genome level vs mapping: bonsai Spearman {common.fmt(e['bonsai_leaf_reads']['spearman'])} "
                     f"(log Pearson {common.fmt(e['bonsai_leaf_reads']['pearson_log10'])}); gather Spearman "
                     f"{common.fmt(e['gather']['spearman'])} (log Pearson {common.fmt(e['gather']['pearson_log10'])})")
        lines.append(f"  - E10a species level vs mapping: bonsai Spearman {common.fmt(e['species_bonsai']['spearman'])}; gather "
                     f"{common.fmt(e['species_gather']['spearman'])} ({len(keys)} species)")
        lines.append(f"  - E10b reads resolved to a genome leaf: {common.fmt(res[sample]['leaf_frac'], 2)} of classified")
        if e10c:
            lines.append(f"  - E10c {len(fp)} gather genomes with mapping covering < 50%: AUROC (well-covered vs those) "
                         + "; ".join(f"{k} {common.fmt(v)}" for k, v in e10c.items()))
    common.write_results("r2", res, lines)
