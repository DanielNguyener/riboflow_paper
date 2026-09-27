#!/usr/bin/env python3
"""Figure 5A: the per-gene read categories, one route-explicit seven-segment bar per gene.

Reads `gene_partition_route7.tsv` (+ .json); the ten-to-seven fold lives with the chain
in `read_categories/gene_read_partition_lib.py` (`ROUTE7_SEGMENTS`, `prepare_route_explicit`).
"""
from __future__ import annotations

import argparse
import csv
import json
import os
import sys
from pathlib import Path

import numpy as np

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, HERE)
sys.path.insert(0, os.path.join(os.path.dirname(HERE), "read_categories"))
import categories  # noqa: E402
from gene_read_partition_lib import ROUTE7_SEGMENTS  # noqa: E402
from panel_style import die  # noqa: E402

GENE_ORDER = ("COMT", "GAPDH", "LRRFIP1")

#: two section key for the hatched design  colour = route and uniqueness  hatch = mechanism
ROUTE7_KEY = (
    tuple((abbr, colour, None) for abbr, colour in categories.KEY),
    # mechanism entries name the biology alone  wording matches figure 4C panel title
    (("Protein-coding–pseudogene ties", "#ffffff", "//"),
     ("Alternative exon", "#ffffff", "..")),
)


BAR_HEIGHT = 0.62
ROW_HEIGHT = 1.05
ROW_MARGIN = 2.3
PAGE_WIDTH = 9.5


def draw(prepared, title=None, figsize=None, label_threshold=6.0, xlabel=None,
         compact=False, title_size=None, bar_height=None):
    import matplotlib.pyplot as plt
    from matplotlib.patches import Patch
    import panel_style as ps

    ps.apply_rcparams()
    entries = prepared["entries"]
    figsize = figsize or (PAGE_WIDTH, ROW_HEIGHT * len(entries) + ROW_MARGIN)
    figure, axis = plt.subplots(figsize=figsize)
    y = np.arange(len(entries))[::-1]

    for yi, entry in zip(y, entries):
        left = 0.0
        for key, _label, colour, text_colour, hatch in ROUTE7_SEGMENTS:
            width = entry["pct"][key]
            axis.barh(yi, width, left=left, color=colour, edgecolor="white",
                      linewidth=0.6, height=bar_height or BAR_HEIGHT)
            if hatch and width > 0:
                # hatch colour rides the artist EDGE colour  so draw a second fill less
                # bar  linewidth=0 keeps the overlay from doubling the boundary
                axis.barh(yi, width, left=left, fill=False, hatch=hatch,
                          edgecolor="white", linewidth=0.0,
                          height=bar_height or BAR_HEIGHT)
            if width >= label_threshold:
                axis.text(left + width / 2, yi, "%.0f%%" % width, va="center",
                          ha="center", fontsize=ps.FONT_ANNOTATION, color=text_colour,
                          zorder=6, bbox=dict(boxstyle="round,pad=0.12", fc=colour,
                                              ec="none") if hatch else None)
            left += width

    axis.set_yticks(list(y))
    # compact uses one line labels and drops the read count  belongs in the caption
    axis.set_yticklabels(
        [e["gene_name"] if compact
         else "%s\n%s reads" % (e["gene_name"], format(e["n_union"], ","))
         for e in entries], fontsize=ps.FONT_TICK)
    axis.set_ylim(-0.6, len(entries) - 0.4)
    axis.set_xlim(0, 100)
    axis.set_xlabel(
        xlabel or ("% of read IDs at this gene" if compact
                   else "% of the read IDs aligning to this gene on either route"),
        fontsize=ps.FONT_LABEL)
    axis.grid(axis="x", alpha=0.15)
    if title is None:
        title = entries[0]["sample"] if entries else ""
    axis.set_title(title, fontsize=title_size or ps.FONT_TITLE, loc="left",
                   fontweight="normal", pad=2.0)

    figure.tight_layout()
    # one band  the five fates then the two mechanism hatches
    handles = [Patch(facecolor=colour, edgecolor="#666666" if hatch else "white",
                     hatch=hatch, linewidth=0.6 if hatch else 0.4, label=label)
               for members in ROUTE7_KEY for label, colour, hatch in members]
    legend = ps.legend_below(
        axis, handles=handles, ncol=len(handles),
        fontsize=ps.FONT_ANNOTATION, handlelength=1.3, columnspacing=1.2,
        handletextpad=0.5, labelspacing=0.3, borderpad=0.0, pad_pt=2.0)
    return figure, axis, [legend]


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
        # percentages recomputed from counts  not trusted from the text column
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

    if args.font_size:
        ps.FONT_TITLE = ps.FONT_LABEL = args.font_size
        ps.FONT_TICK = ps.FONT_ANNOTATION = ps.FONT_INSET = args.font_size

    for entry in prepared["entries"]:
        print("[panel] %-8s union %5d  %s"
              % (entry["gene_name"], entry["n_union"],
                 "  ".join("%s=%d" % (k.replace("r7_", ""), entry["counts"][k])
                           for k, _l, _c, _t, _h in ROUTE7_SEGMENTS)))
    figure, axis, extra = draw(
        prepared, title=args.title if args.title is not None else meta["gsm"],
        figsize=tuple(args.figsize) if args.figsize else None,
        xlabel=args.xlabel, compact=args.compact, title_size=args.title_size,
        bar_height=args.bar_height)
    # full width panel has room for a tick every 10 %
    axis.set_xticks(range(0, 101, 10))
    ps.save(figure, args.output, ps.resolve_formats(args.formats), force=args.force,
            extra_artists=extra)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
