#!/usr/bin/env python3
"""Per-read gene partition dump -> the seven-segment table (`gene_partition_route7.{tsv,json}`).

Folds the dump through `panels/plot_gene_read_partition.prepare_route_explicit`, so the
fold uses the same semantics as the panel. Run with `python` (3.9).

Segments (unit: read IDs; denominator: the union of read IDs at the gene on either route).
"Shared" and "genome-only" are BAM-presence terms over the whole library;
"transcriptome-only" is gene-local.

    r7_shared_unique        shared, genome-unique
    r7_shared_multi_pp      shared, genome-multimapping, protein-coding/pseudogene tie
    r7_shared_multi_other   shared, other genome-multimapping
    r7_gonly_unique_omit    genome-only, genome-unique, on an exon the selected isoform omits
    r7_gonly_unique_other   genome-only, other genome-unique
    r7_gonly_multi          genome-only, genome-multimapping
    r7_txonly               transcriptome-only at the gene

Validated counts (GSM2100602), in the order above: COMT 1084/0/34/105/37/26/9 (union
1,295); GAPDH 1057/2207/805/63/49/64/115 (4,360); LRRFIP1 281/40/18/755/27/388/16 (1,525).
Counts other than these are rejected.
"""
from __future__ import annotations

import argparse
import csv
import json
import os
import sys

HERE = os.path.dirname(os.path.abspath(__file__))
CODE = os.path.dirname(HERE)
sys.path.insert(0, os.path.join(CODE, "common"))
import inputs as paths  # noqa: E402

DEFAULT_OUTPUT = os.path.join(paths.REPO, "results", "alignment_fate", "gene_partition_route7")

#: The validated counts (Fig6A segment semantics audit). Refusing to write anything else
#: is the point: a short bar is the one failure a stacked figure cannot show you.
EXPECTED_COUNTS = {
    "COMT":    {"r7_shared_unique": 1084, "r7_shared_multi_pp": 0,
                "r7_shared_multi_other": 34, "r7_gonly_unique_omit": 105,
                "r7_gonly_unique_other": 37, "r7_gonly_multi": 26, "r7_txonly": 9},
    "GAPDH":   {"r7_shared_unique": 1057, "r7_shared_multi_pp": 2207,
                "r7_shared_multi_other": 805, "r7_gonly_unique_omit": 63,
                "r7_gonly_unique_other": 49, "r7_gonly_multi": 64, "r7_txonly": 115},
    "LRRFIP1": {"r7_shared_unique": 281, "r7_shared_multi_pp": 40,
                "r7_shared_multi_other": 18, "r7_gonly_unique_omit": 755,
                "r7_gonly_unique_other": 27, "r7_gonly_multi": 388, "r7_txonly": 16},
}
EXPECTED_UNION = {"COMT": 1295, "GAPDH": 4360, "LRRFIP1": 1525}
GENE_ORDER = ("COMT", "GAPDH", "LRRFIP1")

COLUMNS = ("sample", "gsm", "gene_order", "gene_name", "transcript_id", "n_union",
           "segment_order", "segment_key", "segment_label", "n_reads", "pct_of_union")


def fold(reads_path, sample, genes):
    sys.path.insert(0, os.path.join(CODE, "panels"))
    import pandas as pd
    import plot_gene_read_partition as root

    prepared = root.prepare_route_explicit(reads_path, sample=sample, genes=genes)
    # Re-check the partition invariants from the raw frame, independently of the root module.
    frame = pd.read_csv(reads_path, sep="\t")
    frame = frame[frame["sample"].astype(str) == str(sample)]
    for entry in prepared["entries"]:
        gene = entry["gene_name"]
        rows = frame[frame["gene_name"] == gene]
        if rows["read_id"].duplicated().any():
            paths.die("%s: duplicated read id in %s" % (gene, reads_path))
        if sum(entry["counts"].values()) != entry["n_union"]:
            paths.die("%s: segments sum to %d, union is %d"
                      % (gene, sum(entry["counts"].values()), entry["n_union"]))
        if entry["n_union"] != len(rows):
            paths.die("%s: union %d but %d reads in the dump"
                      % (gene, entry["n_union"], len(rows)))
        if abs(sum(entry["pct"].values()) - 100.0) > 1e-9:
            paths.die("%s: percentages sum to %.9f" % (gene, sum(entry["pct"].values())))
    segments = [{"key": k, "label": l, "colour": c, "text_colour": t, "hatch": h}
                for k, l, c, t, h in root.ROUTE7_SEGMENTS]
    mapping = {"shared_unique": list(root._R7_UNIQUE_SHARED),
               "multi": list(root._R7_MULTI), "multi_pseudogene_tie": list(root._R7_MULTI_PP),
               "absent_omitted_exon": list(root._R7_ABSENT_OMIT),
               "absent_other": list(root._R7_ABSENT_OTHER), "txonly": list(root._R7_TXONLY)}
    return prepared, segments, mapping


def check_expected(prepared):
    order = [e["gene_name"] for e in prepared["entries"]]
    if tuple(order) != GENE_ORDER:
        paths.die("gene order is %r, expected %r" % (order, GENE_ORDER))
    for entry in prepared["entries"]:
        gene = entry["gene_name"]
        if entry["n_union"] != EXPECTED_UNION[gene]:
            paths.die("%s union %d != expected %d"
                      % (gene, entry["n_union"], EXPECTED_UNION[gene]))
        for key, expected in EXPECTED_COUNTS[gene].items():
            if entry["counts"][key] != expected:
                paths.die("%s %s = %d, expected %d"
                          % (gene, key, entry["counts"][key], expected))


def write(prepared, segments, mapping, reads_path, sample, gsm, output_stem,
          record_input_paths):
    label_of = {s["key"]: s["label"] for s in segments}
    with open(output_stem + ".tsv", "w", newline="") as handle:
        writer = csv.writer(handle, delimiter="\t", lineterminator="\n")
        writer.writerow(COLUMNS)
        for g_index, entry in enumerate(prepared["entries"]):
            for s_index, seg in enumerate(segments):
                key = seg["key"]
                writer.writerow([sample, gsm, g_index, entry["gene_name"],
                                 entry["transcript_id"], entry["n_union"], s_index, key,
                                 label_of[key], entry["counts"][key],
                                 "%.10f" % entry["pct"][key]])
    meta = {
        "sample": sample, "gsm": gsm,
        "gene_order": [e["gene_name"] for e in prepared["entries"]],
        "transcripts": {e["gene_name"]: e["transcript_id"] for e in prepared["entries"]},
        "n_union": {e["gene_name"]: e["n_union"] for e in prepared["entries"]},
        "segments": segments,
        "category_to_segment_group": mapping,
        "semantics": {
            "shared": "read id present in BOTH dedup BAMs (global, not gene-local)",
            "genome_only": "read id absent from the transcriptome dedup BAM entirely",
            "txonly": "no genome alignment AT THIS GENE; the read may align elsewhere",
            "unit": "read ids; denominator = union of read ids at the gene on either route"},
        "source": {"file": os.path.basename(reads_path),
                   "sha256": paths.sha256_of(reads_path),
                   "n_rows": sum(1 for _ in open(reads_path)) - 1,
                   "generator": "code/alignment_fate/build_gene_read_partition.py "
                                "--dump-reads"},
        "builder": "code/alignment_fate/build_gene_partition_data.py",
    }
    if record_input_paths:
        meta["source"]["path"] = os.path.abspath(reads_path)
    with open(output_stem + ".json", "w") as handle:
        handle.write(json.dumps(meta, indent=2, sort_keys=True) + "\n")


def main(argv=None):
    parser = argparse.ArgumentParser(
        description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument("--reads", required=True, help="*_gene_read_partition_reads.tsv")
    parser.add_argument("--sample", default="HeLa")
    parser.add_argument("--gsm", default="GSM2100602")
    parser.add_argument("--genes", default=",".join(GENE_ORDER))
    parser.add_argument("--output", default=DEFAULT_OUTPUT,
                        help="output stem (.tsv and .json are appended)")
    parser.add_argument("--record-input-paths", action="store_true")
    parser.add_argument("--force", action="store_true")
    args = parser.parse_args(argv)

    reads = str(paths.repo_path(args.reads))
    if not os.path.exists(reads):
        paths.die("%s is missing; run code/alignment_fate/build_gene_read_partition.py "
                  "--dump-reads first" % reads)
    if os.path.exists(args.output + ".tsv") and not args.force:
        paths.die("%s.tsv exists; pass --force" % args.output)
    genes = [g.strip() for g in args.genes.split(",") if g.strip()]
    prepared, segments, mapping = fold(reads, args.sample, genes)
    check_expected(prepared)
    os.makedirs(os.path.dirname(args.output), exist_ok=True)
    write(prepared, segments, mapping, reads, args.sample, args.gsm, args.output,
          args.record_input_paths)
    for entry in prepared["entries"]:
        print("[tables] %-8s union %5d  %s"
              % (entry["gene_name"], entry["n_union"],
                 "  ".join("%s=%d" % (k.replace("r7_", ""), entry["counts"][k])
                           for k in entry["counts"])))
    print("[tables] wrote %s.tsv / .json" % args.output)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
