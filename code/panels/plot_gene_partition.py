#!/usr/bin/env python3
"""The route-explicit seven-segment gene partition, from the compact table.

Reads `gene_partition_route7.tsv` (+ .json) and reuses `plot_gene_read_partition.draw()`
with the route-explicit fold (see ROUTE7_SEGMENTS there). Run with `python` (3.9).
"""
from __future__ import annotations

import argparse
import csv
import json
import os
import sys
from pathlib import Path

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, HERE)


def die(message):
    raise SystemExit("error: %s" % message)

GENE_ORDER = ("COMT", "GAPDH", "LRRFIP1")


def load_compact(table, meta_path, genes=None):
    """The compact table -> {"entries": [...]} in the `draw()` shape, plus the JSON."""
    with open(meta_path) as handle:
        meta = json.load(handle)
    keys = [s["key"] for s in meta["segments"]]
    by_gene = {}
    with open(table) as handle:
        for row in csv.DictReader(handle, delimiter="\t"):
            entry = by_gene.setdefault(row["gene_name"], {
                "gene_name": row["gene_name"], "transcript_id": row["transcript_id"],
                "sample": row["sample"], "n_union": int(row["n_union"]),
                "pct": {}, "counts": {}, "_order": int(row["gene_order"])})
            entry["counts"][row["segment_key"]] = int(row["n_reads"])
            entry["pct"][row["segment_key"]] = float(row["pct_of_union"])
    order = sorted(by_gene, key=lambda g: by_gene[g]["_order"])
    if genes:
        unknown = [g for g in genes if g not in by_gene]
        if unknown:
            die("%r not in %s" % (unknown, table))
        order = list(genes)
    entries = []
    for gene in order:
        entry = by_gene[gene]
        if set(entry["counts"]) != set(keys):
            die("%s: segments %r != %r" % (gene, sorted(entry["counts"]), sorted(keys)))
        if sum(entry["counts"].values()) != entry["n_union"]:
            die("%s: counts sum to %d, union %d"
                      % (gene, sum(entry["counts"].values()), entry["n_union"]))
        # Percentages recomputed from counts, not trusted from the text column.
        entry["pct"] = {k: 100.0 * entry["counts"][k] / entry["n_union"] for k in keys}
        if abs(sum(entry["pct"].values()) - 100.0) > 1e-9:
            die("%s: percentages sum to %.9f" % (gene, sum(entry["pct"].values())))
        del entry["_order"]
        entries.append(entry)
    return {"entries": entries}, meta


def main(argv=None):
    parser = argparse.ArgumentParser(description=__doc__,
                                     formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument("--derived-table", required=True, help="gene_partition_route7.tsv")
    parser.add_argument("--derived-meta", required=True, help="gene_partition_route7.json")
    parser.add_argument("--genes", default=None,
                        help="comma-separated gene names, in drawing order")
    parser.add_argument("--title", help="default: the GSM recorded in the JSON")
    parser.add_argument("--xlabel")
    parser.add_argument("--figsize", nargs=2, type=float)
    parser.add_argument("--font-size", type=float)
    parser.add_argument("--title-size", type=float)
    parser.add_argument("--bar-height", type=float)
    parser.add_argument("--compact", action="store_true")
    parser.add_argument("--output", required=True, type=Path, help="stem, no extension")
    parser.add_argument("--format", dest="formats", default="pdf")
    parser.add_argument("--force", action="store_true")
    args = parser.parse_args(argv)

    genes = [g.strip() for g in args.genes.split(",") if g.strip()] if args.genes else None
    prepared, meta = load_compact(args.derived_table, args.derived_meta, genes)
    if tuple(e["gene_name"] for e in prepared["entries"]) != GENE_ORDER:
        die("gene order %r != %r" % ([e["gene_name"] for e in prepared["entries"]],
                                            GENE_ORDER))

    import panel_style as ps
    import plot_gene_read_partition as root

    if args.font_size:
        ps.FONT_TITLE = ps.FONT_LABEL = args.font_size
        ps.FONT_TICK = ps.FONT_ANNOTATION = ps.FONT_INSET = args.font_size

    for entry in prepared["entries"]:
        print("[panel] %-8s union %5d  %s"
              % (entry["gene_name"], entry["n_union"],
                 "  ".join("%s=%d" % (k.replace("r7_", ""), entry["counts"][k])
                           for k, _l, _c, _t, _h in root.ROUTE7_SEGMENTS)))
    figure, axis, extra = root.draw(
        prepared, title=args.title if args.title is not None else meta["gsm"],
        figsize=tuple(args.figsize) if args.figsize else None,
        xlabel=args.xlabel, compact=args.compact, title_size=args.title_size,
        bar_height=args.bar_height, grouped_key=True,
        segments=[(k, l, c, x, h) for k, l, c, x, h in root.ROUTE7_SEGMENTS])
    # A full-width panel has room for a tick every 10 %.
    axis.set_xticks(range(0, 101, 10))
    ps.save(figure, args.output, ps.resolve_formats(args.formats), force=args.force,
            extra_artists=extra)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
