#!/usr/bin/env python3
"""How much of each gene's selected transcript is duplicated, base for base, by the
selected transcript of ANOTHER gene in the transcriptome reference.

    python code/clustering/validate_cluster_reference_duplication.py \\
        --clusters results/clustering/HeLa.post_dedup.clusters_k4.tsv \\
        --output results/clustering [--gtf G --appris A]

For every entry of the transcriptome reference (the APPRIS lengths table: one selected
transcript per gene id), its exons are taken from the GENCODE v34 GTF. Per chromosome and
strand, every exonic genomic base is counted over all reference entries. Then per entry:

    exonic_nt              merged exon length of the selected transcript
    duplicated_nt          its exonic bases also inside an exon of >= 1 OTHER reference
                           entry on the same strand
    duplicated_fraction    duplicated_nt / exonic_nt
    n_duplicating_entries  how many other entries share >= 1 base with it
    duplicating_entries    their transcript ids (gene names in the next column)
    duplicated_nt_any_strand, duplicated_fraction_any_strand
                           the same counting both strands (antisense overlaps; the
                           pipeline aligns with --norc, so these do NOT tie in bowtie2)

Same strand is the headline because bowtie2 `--norc` only ever aligns a read to the
forward strand of an entry, so only same-strand duplicates can produce the MAPQ 0-1 tie
that a read-level audit of the genome-only reads found.

Outputs, under <output> with the cluster table's stem:

    <stem>.reference_duplication_entries.tsv     every reference entry, with the cluster
                                                 where the gene is clustered (else NA)
    <stem>.reference_duplication_by_cluster.tsv  per cluster (and all clustered, and every
                                                 reference entry): n, % with any duplicated
                                                 base, % with >= 50 % duplicated, quartiles,
                                                 mean; odds ratio of "any duplication" for
                                                 the cluster vs the other clusters (Fisher)
    <stem>.reference_duplication_tests.tsv       Kruskal-Wallis over clusters and Holm-
                                                 adjusted pairwise Mann-Whitney U
    <stem>.reference_duplication_sources.json

The entries table is Figure 6G's input (`panels/plot_read_fate_clusters.py`).
"""
from __future__ import annotations

import json
import os
from collections import defaultdict
from datetime import datetime, timezone
from pathlib import Path

import numpy as np
import pandas as pd

import validate_common as vc
from validate_common import inputs, merge

log = inputs.make_log("clustering/duplication")


def reference_entries(appris_path):
    """transcript id -> gene name, from the APPRIS lengths table's headers."""
    entries = {}
    with open(appris_path) as handle:
        for line in handle:
            fields = line.split("\t")[0].split("|")
            if len(fields) > 6 and fields[0].startswith("ENST"):
                entries[fields[0]] = fields[5]
    if not entries:
        raise SystemExit("%s holds no reference entries" % appris_path)
    return entries


def parse_exons(gtf, tids):
    exons = defaultdict(list)
    where = {}
    for f in vc.gtf_records(gtf, ("exon",)):
        tid = f[8].split('transcript_id "', 1)[1].split('"', 1)[0]
        if tid not in tids:
            continue
        exons[tid].append((int(f[3]) - 1, int(f[4])))
        where[tid] = (f[0], f[6])
    return {t: merge(v) for t, v in exons.items()}, where


def duplication(exons, where, strand_aware=True):
    """Per entry: duplicated_nt and the set of other entries sharing a base."""
    groups = defaultdict(list)
    for tid, (chrom, strand) in where.items():
        groups[(chrom, strand if strand_aware else ".")].append(tid)
    dup_nt, partners = {}, defaultdict(set)
    for key, tids in groups.items():
        # sweep: coverage count over sorted boundaries
        events = []
        for tid in tids:
            for s, e in exons[tid]:
                events.append((s, 1, tid))
                events.append((e, -1, tid))
        events.sort(key=lambda x: (x[0], x[1]))
        active = set()
        last = None
        # per-entry duplicated length: accumulate segments where >= 2 active
        for pos, delta, tid in events:
            if last is not None and pos > last and len(active) >= 2:
                length = pos - last
                for t in active:
                    dup_nt[t] = dup_nt.get(t, 0) + length
                    partners[t].update(active - {t})
            if delta == 1:
                active.add(tid)
            else:
                active.discard(tid)
            last = pos
    return dup_nt, partners


def summarize(table):
    clustered = table[table["cluster"].notna()].copy()
    clustered["cluster"] = clustered["cluster"].astype(int)
    clusters = sorted(clustered["cluster"].unique())
    rows = []
    groups = [("all_reference_entries", table), ("all_clustered", clustered)] + \
             [(str(c), clustered[clustered["cluster"] == c]) for c in clusters]
    for label, g in groups:
        v = g["duplicated_fraction"].to_numpy(float)
        row = {"cluster": label, "n": len(g),
               "pct_any_duplication": 100.0 * np.mean(v > 0) if len(v) else np.nan,
               "pct_half_or_more": 100.0 * np.mean(v >= 0.5) if len(v) else np.nan,
               **vc.quartiles(v),
               "duplicated_nt_median": g["duplicated_nt"].median(),
               "pct_any_duplication_any_strand": 100.0 * (g["duplicated_fraction_any_strand"] > 0).mean()
               if len(g) else np.nan}
        if label in [str(c) for c in clusters]:
            inside = g["duplicated_fraction"] > 0
            rest = clustered.loc[clustered["cluster"] != int(label), "duplicated_fraction"] > 0
            orr, lo, hi, p, _cells = vc.odds_ratio(inside, rest)
            row.update({"any_duplication_or": orr, "any_duplication_or_ci_low": lo,
                        "any_duplication_or_ci_high": hi, "any_duplication_fisher_p": p})
        rows.append(row)
    summary = pd.DataFrame(rows)

    per = {c: clustered.loc[clustered["cluster"] == c, "duplicated_fraction"].to_numpy(float) for c in clusters}
    return summary, pd.DataFrame(vc.cluster_tests(per, clusters))


def main(argv=None):
    p = vc.parser(__doc__)
    p.add_argument("--appris", type=Path)
    a = p.parse_args(argv)
    gtf = vc.resolve_gtf(a)
    appris = a.appris or os.environ.get("RIBOFLOW_PAPER_APPRIS") or inputs._local_config().get("appris")
    if not appris or not Path(appris).exists():
        inputs.die("appris not configured or missing")

    a.output.mkdir(parents=True, exist_ok=True)
    stem, clusters = vc.load_clusters(a.clusters)
    entries = reference_entries(appris)
    log("%d reference entries; %d clustered genes" % (len(entries), len(clusters)))
    exons, where = parse_exons(gtf, set(entries))
    missing = [t for t in entries if t not in exons]
    if missing:
        log("note: %d reference entries have no exon in the GTF and are skipped: %s"
            % (len(missing), ", ".join(missing[:5])))
    dup_same, partners_same = duplication(exons, where, strand_aware=True)
    dup_any, _partners_any = duplication(exons, where, strand_aware=False)

    cluster_of = dict(zip(clusters["transcript_id"], clusters["cluster"]))
    rows = []
    for tid in exons:
        total = sum(e - s for s, e in exons[tid])
        d = dup_same.get(tid, 0)
        d_any = dup_any.get(tid, 0)
        rows.append({"gene": entries[tid], "transcript_id": tid, "chrom": where[tid][0],
                     "strand": where[tid][1], "cluster": cluster_of.get(tid, np.nan),
                     "exonic_nt": total, "duplicated_nt": d, "duplicated_fraction": d / total,
                     "n_duplicating_entries": len(partners_same.get(tid, ())),
                     "duplicating_entries": ";".join(sorted(partners_same.get(tid, ()))),
                     "duplicating_genes": ";".join(entries[t] for t in sorted(partners_same.get(tid, ()))),
                     "duplicated_nt_any_strand": d_any, "duplicated_fraction_any_strand": d_any / total})
    table = pd.DataFrame(rows).sort_values(["duplicated_fraction", "gene"], ascending=[False, True])
    not_found = [t for t in clusters["transcript_id"] if t not in exons]
    if not_found:
        log("note: %d clustered transcripts are not reference entries with exons" % len(not_found))
    summary, tests = summarize(table)

    entries_tsv = a.output / ("%s.reference_duplication_entries.tsv" % stem)
    by_cluster_tsv = a.output / ("%s.reference_duplication_by_cluster.tsv" % stem)
    tests_tsv = a.output / ("%s.reference_duplication_tests.tsv" % stem)
    sources_json = a.output / ("%s.reference_duplication_sources.json" % stem)
    table.to_csv(entries_tsv, sep="\t", index=False, float_format="%.6f", na_rep="NA")
    summary.to_csv(by_cluster_tsv, sep="\t", index=False, float_format="%.4g", na_rep="NA")
    tests.to_csv(tests_tsv, sep="\t", index=False, float_format="%.4g")
    sources_json.write_text(json.dumps({
        "generated_utc": datetime.now(timezone.utc).isoformat(timespec="seconds"),
        "clusters": vc.describe(a.clusters), "gtf": vc.describe(gtf), "appris": vc.describe(appris),
        "n_reference_entries": len(entries), "n_entries_with_exons": len(exons)}, indent=2) + "\n")
    for path in (entries_tsv, by_cluster_tsv, tests_tsv, sources_json):
        log("wrote %s" % path)
    with pd.option_context("display.width", 200, "display.max_columns", 30):
        print(summary.to_string(index=False))
        print(tests.to_string(index=False))


if __name__ == "__main__":
    main()
