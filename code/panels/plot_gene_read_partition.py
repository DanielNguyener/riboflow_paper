#!/usr/bin/env python3
"""Test figure -- every read at a gene, on either route, in one bar.

Denominator is the UNION of read IDs aligning to the gene on either route.
Exploratory, not a manifest panel; table from code/alignment_fate/build_gene_read_partition.py.
"""
from __future__ import annotations

import argparse
import sys
from pathlib import Path

import numpy as np
import pandas as pd

HERE = Path(__file__).resolve().parent
REQUIRED = ("sample", "gene_name", "transcript_id", "category", "n_reads", "pct_of_union")

#: The drawn segments: (key, label, colour, text colour, the categories folded into it).
#: Folds the sixteen table categories to six; the table/read dump keep all sixteen, so a
#: finer fold needs only a different SEGMENTS. Colours follow Figure 5's vocabulary.
SEGMENTS = (
    ("genome_unique", "genome-unique",
     "#a6d96a", "black", ("genome_unique_shared_concordant",
                          "genome_unique_shared_discordant",
                          "genome_unique_txome_other_transcript",
                          "genome_unique_txome_multimapped")),
    ("multi_pseudogene", "genome-multi, pseudogene tie",
     "#1a7d1a", "white", ("genome_multi_primary_in_gene_pseudogene_tie",
                          "genome_multi_primary_pseudogene")),
    ("multi_other", "genome-multi",
     "#8fbf8f", "black", ("genome_multi_primary_in_gene_no_pseudogene_tie",
                          "genome_multi_primary_elsewhere_other",
                          "genome_multi_primary_lost_in_dedup")),
    ("reach_isoform", "genome-only, exon of a nonselected isoform",
     "#7fb9da", "black", ("genome_unique_absent_nonselected_isoform_exon",)),
    ("reach_other", "genome-only, other",
     "#cfe3ef", "black", ("genome_unique_absent_splice_junction",
                          "genome_unique_absent_representable",
                          "genome_unique_absent_pseudogene",
                          "genome_unique_absent_other")),
    ("txome_only", "transcriptome-only",
     "#cc3d3d", "white", ("txome_only_genome_absent",
                          "txome_only_genome_elsewhere")),
)

#: Published key wording. Caveats (the "Both routes; genome-multi" folds are not literally
#: all both-routes; "Transcriptome only" means only AT THIS GENE) belong in the caption.
ROUTE_LABELS = {
    "genome_unique": "Both routes: genome-unique",
    "multi_pseudogene": "Both routes: genome-multi, pseudogene tie",
    "multi_other": "Both routes: genome-multi, other",
    "reach_isoform": "Genome only: omitted exon",
    "reach_other": "Genome only: other",
    "txome_only": "Transcriptome only",
}

SHORT_LABELS = {
    "genome_unique": "genome-unique",
    "multi_pseudogene": "multi: pseudogene",
    "multi_other": "multi: other",
    "reach_isoform": "genome-only: alt. isoform",
    "reach_other": "genome-only: other",
    "txome_only": "transcriptome-only",
}

#: The route-explicit seven-segment fold, folded from the PER-READ dump -- the tidy table
#: never records a genome multimapper's transcriptome status. "Shared" is read-level
#: presence in both BAMs, not "assigned to this gene by both routes".
ROUTE7_SEGMENTS = (
    ("r7_shared_unique", "Genome-unique", "#a6d96a", "black", None),
    ("r7_shared_multi_pp", "Genome-multi, pseudogene tie", "#1a7d1a", "white", "//"),
    ("r7_shared_multi_other", "Genome-multi, other", "#1a7d1a", "white", None),
    ("r7_gonly_unique_omit", "Genome-unique, omitted exon", "#7fb9da", "black", ".."),
    ("r7_gonly_unique_other", "Genome-unique, other", "#7fb9da", "black", None),
    ("r7_gonly_multi", "Genome-multi", "#0d57a1", "white", None),
    ("r7_txonly", "Transcriptome only", "#cc3d3d", "white", None),
)

#: The two-section key for the hatched design: colour = route/uniqueness, hatch = mechanism.
ROUTE7_KEY = (
    (("Shared, unique", "#a6d96a", None),
     ("Shared, multimapping", "#1a7d1a", None),
     ("Genome-only, unique", "#7fb9da", None),
     ("Genome-only, multimapping", "#0d57a1", None),
     ("Transcriptome-only at gene", "#cc3d3d", None)),
    # Mechanism entries name the biology alone; wording matches Figure 5C's panel title.
    (("Protein-coding–pseudogene ties", "#ffffff", "//"),
     ("Alternative exon", "#ffffff", "..")),
)

#: Raw chain category -> genome-status half of the seven-way fold; the transcriptome half
#: comes from the dump's per-read `txome_primary_transcript`.
_R7_MULTI = ("genome_multi_primary_in_gene_pseudogene_tie", "genome_multi_primary_pseudogene",
             "genome_multi_primary_in_gene_no_pseudogene_tie",
             "genome_multi_primary_elsewhere_other", "genome_multi_primary_lost_in_dedup")
_R7_MULTI_PP = _R7_MULTI[:2]
_R7_UNIQUE_SHARED = ("genome_unique_shared_concordant", "genome_unique_shared_discordant",
                     "genome_unique_txome_other_transcript",
                     "genome_unique_txome_multimapped")
_R7_ABSENT_OMIT = ("genome_unique_absent_nonselected_isoform_exon",)
_R7_ABSENT_OTHER = ("genome_unique_absent_splice_junction",
                    "genome_unique_absent_representable",
                    "genome_unique_absent_pseudogene", "genome_unique_absent_other")
_R7_TXONLY = ("txome_only_genome_absent", "txome_only_genome_elsewhere")


def _route7_segment(category, txome_present):
    """One read -> one of the seven keys. Raises on an unknown category."""
    if category in _R7_UNIQUE_SHARED:
        return "r7_shared_unique"          # every such category conditions on presence
    if category in _R7_MULTI:
        if not txome_present:
            return "r7_gonly_multi"
        return ("r7_shared_multi_pp" if category in _R7_MULTI_PP
                else "r7_shared_multi_other")
    if category in _R7_ABSENT_OMIT:
        return "r7_gonly_unique_omit"
    if category in _R7_ABSENT_OTHER:
        return "r7_gonly_unique_other"
    if category in _R7_TXONLY:
        return "r7_txonly"
    raise SystemExit("unknown chain category %r" % category)


def prepare_route_explicit(reads_path, sample=None, genes=None):
    """The per-read dump -> seven-way entries, same shape `draw` consumes; re-asserts the
    partition invariants."""
    frame = pd.read_csv(reads_path, sep="\t")
    needed = ("sample", "gene_name", "transcript_id", "read_id", "category",
              "txome_primary_transcript")
    missing = [c for c in needed if c not in frame.columns]
    if missing:
        raise SystemExit("%s lacks column(s) %s -- pass the *_reads.tsv dump, not the "
                         "tidy table" % (reads_path, ", ".join(missing)))
    if sample:
        frame = frame[frame["sample"].astype(str) == str(sample)]
    txp = frame["txome_primary_transcript"].fillna("").astype(str) != ""
    frame = frame.assign(_seg=[_route7_segment(c, p)
                               for c, p in zip(frame["category"], txp)])

    order = list(dict.fromkeys(frame["gene_name"]))
    if genes:
        unknown = [g for g in genes if g not in set(order)]
        if unknown:
            raise SystemExit("%r not in %s" % (unknown, reads_path))
        order = list(dict.fromkeys(genes))

    entries = []
    for gene in order:
        rows = frame[frame["gene_name"] == gene]
        if rows["read_id"].duplicated().any():
            raise SystemExit("%s: duplicated read id in the dump" % gene)
        n_union = len(rows)
        counts = rows["_seg"].value_counts()
        pct = {key: 100.0 * int(counts.get(key, 0)) / n_union
               for key, _l, _c, _t, _h in ROUTE7_SEGMENTS}
        if abs(sum(pct.values()) - 100.0) > 1e-9:
            raise SystemExit("%s: the seven segments sum to %.6f %%" % (gene, sum(pct.values())))
        entries.append({"transcript_id": rows["transcript_id"].iloc[0],
                        "gene_name": gene, "sample": str(rows["sample"].iloc[0]),
                        "n_union": n_union, "pct": pct,
                        "counts": {key: int(counts.get(key, 0))
                                   for key, _l, _c, _t, _h in ROUTE7_SEGMENTS}})
    return {"entries": entries}


BAR_HEIGHT = 0.62
ROW_HEIGHT = 1.05
ROW_MARGIN = 2.3
PAGE_WIDTH = 9.5


def prepare(partition_path, sample=None, genes=None):
    """The table -> one entry per gene, in drawing order, with the folded percentages.

    A missing declared category is an error, never a silent gap.
    """
    sys.path.insert(0, str(HERE))
    import panel_style as ps

    frame = pd.read_csv(partition_path, sep="\t")
    ps.require_columns(frame, REQUIRED, partition_path)
    if sample:
        frame = frame[frame["sample"].astype(str) == str(sample)]
        if frame.empty:
            raise SystemExit("no rows for sample %r in %s" % (sample, partition_path))

    order = list(dict.fromkeys(frame["transcript_id"]))
    if genes:
        wanted, lookup = [], {}
        for tid in order:
            row = frame[frame["transcript_id"] == tid].iloc[0]
            lookup[tid] = tid
            lookup[str(row["gene_name"])] = tid
        for name in genes:
            if name not in lookup:
                raise SystemExit("%r is not in %s" % (name, partition_path))
            wanted.append(lookup[name])
        order = list(dict.fromkeys(wanted))

    declared = set()
    for _key, _label, _colour, _text, categories in SEGMENTS:
        declared |= set(categories)

    entries = []
    for tid in order:
        rows = frame[frame["transcript_id"] == tid]
        present = set(rows["category"])
        missing = declared - present
        if missing:
            raise SystemExit("%s: the table has no row for %s"
                             % (tid, ", ".join(sorted(missing))))
        extra = present - declared
        if extra:
            raise SystemExit("%s: the table has category/ies this panel cannot draw: %s. "
                             "Add them to SEGMENTS rather than dropping them silently."
                             % (tid, ", ".join(sorted(extra))))
        pct = dict(zip(rows["category"], rows["pct_of_union"].astype(float)))
        folded = {key: float(sum(pct[c] for c in categories))
                  for key, _l, _c, _t, categories in SEGMENTS}
        total = sum(folded.values())
        if abs(total - 100.0) > 1e-3:
            raise SystemExit("%s: the drawn segments sum to %.4f %%, not 100" % (tid, total))
        entries.append({
            "transcript_id": tid,
            "gene_name": str(rows["gene_name"].iloc[0]) or tid,
            "sample": str(rows["sample"].iloc[0]),
            "n_union": int(rows["n_reads"].sum()),
            "pct": folded})
    return {"entries": entries}


def draw(prepared, title=None, figsize=None, label_threshold=6.0, xlabel=None,
         compact=False, title_size=None, bar_height=None, route_labels=False,
         segments=None, grouped_key=False):
    import matplotlib.pyplot as plt
    from matplotlib.patches import Patch
    import panel_style as ps

    ps.apply_rcparams()
    entries = prepared["entries"]
    # Default six-way fold, or a caller-supplied segment list (the seven-way fold does this).
    segs = (segments if segments is not None
            else [(k, l, c, x) for k, l, c, x, _cats in SEGMENTS])
    figsize = figsize or (PAGE_WIDTH, ROW_HEIGHT * len(entries) + ROW_MARGIN)
    figure, axis = plt.subplots(figsize=figsize)
    y = np.arange(len(entries))[::-1]

    for yi, entry in zip(y, entries):
        left = 0.0
        for seg in segs:
            key, _label, colour, text_colour = seg[:4]
            hatch = seg[4] if len(seg) > 4 else None
            width = entry["pct"][key]
            axis.barh(yi, width, left=left, color=colour, edgecolor="white",
                      linewidth=0.6, height=bar_height or BAR_HEIGHT)
            if hatch and width > 0:
                # Hatch colour rides the artist's EDGE colour, so draw a second fill-less
                # bar; linewidth=0 keeps the overlay from doubling the boundary.
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
    # `compact` uses one-line labels and drops the read count (belongs in the caption).
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
    if grouped_key:
        # Two columns, filled column-major; the short one is padded so it does not borrow
        # entries from its neighbour.
        rows = max(len(members) for members in ROUTE7_KEY)
        handles = []
        for members in ROUTE7_KEY:
            for label, colour, hatch in members:
                handles.append(Patch(facecolor=colour,
                                     edgecolor="#666666" if hatch else "white",
                                     hatch=hatch, linewidth=0.6 if hatch else 0.4,
                                     label=label))
            handles += [Patch(facecolor="none", edgecolor="none", label=" ")
                        ] * (rows - len(members))
        legend = ps.legend_below(
            axis, handles=handles, ncol=len(ROUTE7_KEY),
            fontsize=ps.FONT_ANNOTATION, handlelength=1.3, columnspacing=2.0,
            handletextpad=0.5, labelspacing=0.3, borderpad=0.0, pad_pt=2.0)
        return figure, axis, [legend]
    # Compact mode uses short labels in a tight three-column block.
    if route_labels:
        labels = [ROUTE_LABELS[key] for key, _l, _c, _t, _cats in SEGMENTS]
    else:
        labels = [SHORT_LABELS[key] if compact else label
                  for key, label, _c, _t, _cats in SEGMENTS]
    handles = [Patch(color=colour, label=text)
               for (_k, _l, colour, _t, _c), text in zip(SEGMENTS, labels)]
    if route_labels:
        # Matplotlib fills a multi-column legend COLUMN-major; interleave so the rows read
        # in bar order left to right.
        handles = [handles[i] for i in (0, 3, 1, 4, 2, 5)]
    legend = ps.legend_below(
        axis,
        handles=handles,
        fontsize=ps.FONT_ANNOTATION, handlelength=1.0 if compact else 1.1, ncol=3,
        columnspacing=0.8 if compact else 2.0,
        handletextpad=0.4 if compact else 0.8,
        labelspacing=0.15 if compact else 0.5,
        borderpad=0.0 if compact else 0.4,
        pad_pt=2.0 if compact else 8.0)
    return figure, axis, [legend]


def main(argv=None):
    parser = argparse.ArgumentParser(
        description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument("--partition", type=Path,
                        help="a *.gene_read_partition.tsv from "
                             "code/alignment_fate/build_gene_read_partition.py")
    parser.add_argument("--sample", help="restrict a multi-sample table to one sample")
    parser.add_argument("--gene", dest="genes", action="append",
                        help="gene name or transcript id; repeat to set the drawing order")
    parser.add_argument("--title")
    parser.add_argument("--xlabel")
    parser.add_argument("--figsize", nargs=2, type=float)
    parser.add_argument("--output", required=True, type=Path)
    parser.add_argument("--format", dest="formats", default="pdf")
    parser.add_argument("--reads", type=Path,
                        help="the per-read *_gene_read_partition_reads.tsv dump; required "
                             "for --route-explicit")
    parser.add_argument("--route-explicit", action="store_true",
                        help="the seven-segment route-explicit fold (Fig6A semantics "
                             "audit, Set 2), folded from --reads with a grouped key")
    parser.add_argument("--route-labels", action="store_true",
                        help="key names that state each segment's verified route status")
    parser.add_argument("--font-size", type=float,
                        help="one point size for every string in this panel")
    parser.add_argument("--bar-height", type=float,
                        help="bar thickness as a fraction of the row step (default 0.62); "
                             "raise it so the in-bar percentages have room")
    parser.add_argument("--title-size", type=float,
                        help="title point size; defaults to the shared title size")
    parser.add_argument("--compact", action="store_true",
                        help="one-line y labels, for a short panel in an assembled figure")
    parser.add_argument("--force", action="store_true")
    args = parser.parse_args(argv)

    sys.path.insert(0, str(HERE))
    import panel_style as ps

    if args.font_size:
        # Set before `draw`, which reads these at call time.
        ps.FONT_TITLE = ps.FONT_LABEL = args.font_size
        ps.FONT_TICK = ps.FONT_ANNOTATION = ps.FONT_INSET = args.font_size

    if not args.partition and not args.route_explicit:
        raise SystemExit("give --partition (tidy fold) or --route-explicit with --reads")
    if args.route_explicit:
        if not args.reads:
            raise SystemExit("--route-explicit needs --reads (the per-read dump); the "
                             "tidy table does not record a multimapper's transcriptome "
                             "status")
        prepared = prepare_route_explicit(args.reads, sample=args.sample,
                                          genes=args.genes)
        for entry in prepared["entries"]:
            print("[panel] %-8s union %5d  %s"
                  % (entry["gene_name"], entry["n_union"],
                     "  ".join("%s=%d" % (k.replace("r7_", ""), entry["counts"][k])
                               for k, _l, _c, _t, _h in ROUTE7_SEGMENTS)))
        figure, _axis, extra = draw(
            prepared, title=args.title,
            figsize=tuple(args.figsize) if args.figsize else None,
            xlabel=args.xlabel, compact=args.compact, title_size=args.title_size,
            bar_height=args.bar_height, grouped_key=True,
            segments=[(k, l, c, x, h) for k, l, c, x, h in ROUTE7_SEGMENTS])
        # Finer grid helps read the 1-3 % segments this fold exists to show.
        _axis.set_xticks(range(0, 101, 10))
        ps.save(figure, args.output, ps.resolve_formats(args.formats), force=args.force,
                extra_artists=extra)
        return 0

    prepared = prepare(args.partition, sample=args.sample, genes=args.genes)
    figure, _axis, extra = draw(prepared, title=args.title,
                                figsize=tuple(args.figsize) if args.figsize else None,
                                xlabel=args.xlabel, compact=args.compact,
                                title_size=args.title_size,
                                bar_height=args.bar_height,
                                route_labels=args.route_labels)
    ps.save(figure, args.output, ps.resolve_formats(args.formats), force=args.force,
            extra_artists=extra)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
