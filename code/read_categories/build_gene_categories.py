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

Validated counts (GSM2100602), in the order above: COMT 1084/0/34/88/54/14/9 (union
1,283); GAPDH 1057/2326/670/53/59/40/117 (4,322); LRRFIP1 281/44/13/747/35/380/16 (1,516).
Genome membership is a top-score placement (the primary, or a secondary tied with its AS) on
the gene's span; every read was re-derived by an independent pysam/GTF recomputation.
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

DEFAULT_OUTPUT = os.path.join(paths.REPO, "results", "read_categories", "gene_partition_route7")

#: validated counts from the gene partition segment semantics audit  refusing to write
#: anything else is the point  a short bar is the one failure a stacked figure cannot show
EXPECTED_COUNTS = {
    "COMT":    {"r7_shared_unique": 1084, "r7_shared_multi_pp": 0,
                "r7_shared_multi_other": 34, "r7_gonly_unique_omit": 88,
                "r7_gonly_unique_other": 54, "r7_gonly_multi": 14, "r7_txonly": 9},
    "GAPDH":   {"r7_shared_unique": 1057, "r7_shared_multi_pp": 2326,
                "r7_shared_multi_other": 670, "r7_gonly_unique_omit": 53,
                "r7_gonly_unique_other": 59, "r7_gonly_multi": 40, "r7_txonly": 117},
    "LRRFIP1": {"r7_shared_unique": 281, "r7_shared_multi_pp": 44,
                "r7_shared_multi_other": 13, "r7_gonly_unique_omit": 747,
                "r7_gonly_unique_other": 35, "r7_gonly_multi": 380, "r7_txonly": 16},
}
EXPECTED_UNION = {"COMT": 1283, "GAPDH": 4322, "LRRFIP1": 1516}
GENE_ORDER = ("COMT", "GAPDH", "LRRFIP1")

COLUMNS = ("sample", "gsm", "gene_order", "gene_name", "transcript_id", "n_union",
           "segment_order", "segment_key", "segment_label", "n_reads", "pct_of_union")


def fold(reads_path, sample, genes):
    import pandas as pd
    import gene_read_partition_lib as root

    prepared = root.prepare_route_explicit(reads_path, sample=sample, genes=genes)
    # recheck partition invariants from the raw frame  independent of the root module
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


def write(prepared, segments, mapping, reads_path, sample, gsm, output_stem):
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
            "txonly": "no top-score genome placement AT THIS GENE; the read may align "
                      "elsewhere",
            "genome_membership": "a top-score genome placement (the primary, or a secondary "
                                 "with the primary's AS) overlapping the gene's full span; a "
                                 "read counts once per gene, and may be in several genes",
            "pseudogene_tie": "a top-score placement with its 5' base in an exon of THIS "
                              "gene and another on a processed pseudogene",
            "omitted_exon": "genome-only unique read with an aligned block on this gene's "
                            "exonic sequence absent from its selected transcript",
            "unit": "read ids; denominator = union of read ids at the gene on either route"},
        "source": {"file": os.path.basename(reads_path),
                   "sha256": paths.sha256_of(reads_path),
                   "n_rows": sum(1 for _ in open(reads_path)) - 1,
                   "generator": "code/read_categories/build_gene_categories.py"},
        "builder": "code/read_categories/build_gene_categories.py",
    }
    with open(output_stem + ".json", "w") as handle:
        handle.write(json.dumps(meta, indent=2, sort_keys=True) + "\n")


def dump_reads(sample, genome_bam, txome_bam, gene_ids, coverage_path, output_dir):
    """Run the chain over the three genes and write the per-read dump it is folded from."""
    sys.path.insert(0, HERE)
    import gene_read_partition_lib as lib

    coverage = None
    if coverage_path:
        sys.path.insert(0, os.path.join(CODE, "coverage"))
        import coverage_schema
        coverage = coverage_schema.open_coverage(coverage_path)
    try:
        print("[tables] %s: resolving %d gene(s)" % (sample, len(gene_ids)), flush=True)
        _wide, _tidy, dump = lib.compute_partition(
            sample, genome_bam, txome_bam, gene_ids=gene_ids, coverage=coverage,
            log=lambda m: print("[tables] %s" % m, flush=True))
    finally:
        if coverage is not None:
            coverage.close()
    os.makedirs(output_dir, exist_ok=True)
    dump_path = os.path.join(output_dir, "%s.gene_read_partition_reads.tsv" % sample)
    dump.to_csv(dump_path, sep="\t", index=False, lineterminator="\n")
    print("[tables] wrote %s (%d rows)" % (dump_path, len(dump)), flush=True)
    return dump_path


def main(argv=None):
    parser = argparse.ArgumentParser(
        description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument("--sample", default="HeLa")
    parser.add_argument("--gsm", default="GSM2100602")
    parser.add_argument("--genome-bam", required=True,
                        help="coordinate-sorted and INDEXED: the gene side is a region fetch")
    parser.add_argument("--transcriptome-bam", required=True)
    parser.add_argument("--gene-id", required=True, help="comma-separated gene IDs")
    parser.add_argument("--coverage", default=None,
                        help="a shared_coverage.h5, used only to resolve gene IDs and names")
    parser.add_argument("--output", default=DEFAULT_OUTPUT,
                        help="output stem (.tsv and .json are appended)")
    parser.add_argument("--force", action="store_true")
    args = parser.parse_args(argv)

    if os.path.exists(args.output + ".tsv") and not args.force:
        paths.die("%s.tsv exists; pass --force" % args.output)
    gene_ids = [g.strip() for g in args.gene_id.split(",") if g.strip()]
    reads = dump_reads(args.sample, args.genome_bam, args.transcriptome_bam, gene_ids,
                       args.coverage, os.path.dirname(args.output) or ".")
    prepared, segments, mapping = fold(reads, args.sample, list(GENE_ORDER))
    check_expected(prepared)
    os.makedirs(os.path.dirname(args.output), exist_ok=True)
    write(prepared, segments, mapping, reads, args.sample, args.gsm, args.output)
    for entry in prepared["entries"]:
        print("[tables] %-8s union %5d  %s"
              % (entry["gene_name"], entry["n_union"],
                 "  ".join("%s=%d" % (k.replace("r7_", ""), entry["counts"][k])
                           for k in entry["counts"])))
    print("[tables] wrote %s.tsv / .json" % args.output)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
