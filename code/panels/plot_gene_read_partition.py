"""Library: every read at a gene, on either route, in one bar (Figure 5A).

Denominator is the UNION of read IDs aligning to the gene on either route. `draw` is the
bar drawer `plot_gene_partition.py` calls; `prepare_route_explicit` and `_route7_segment`
fold the per-read dump from code/alignment_fate/build_gene_read_partition.py.
"""
from __future__ import annotations

import numpy as np
import pandas as pd

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
    (("SH-U", "#a6d96a", None),
     ("SH-M", "#1a7d1a", None),
     ("GO-U", "#7fb9da", None),
     ("GO-M", "#0d57a1", None),
     ("TO", "#cc3d3d", None)),
    # Mechanism entries name the biology alone; wording matches Figure 4C's panel title.
    (("Protein-coding–pseudogene ties", "#ffffff", "//"),
     ("Alternative exon", "#ffffff", "..")),
)

#: Raw chain category -> genome-status half of the seven-way fold; the transcriptome half
#: comes from the dump's per-read `txome_primary_transcript`.
#: A genome-multimapped read at a gene has a top-score placement there (its primary or a
#: score-tied secondary); any other multi category is unknown and raises.
_R7_MULTI = ("genome_multi_top_at_gene_pseudogene_tie",
             "genome_multi_top_at_gene_no_pseudogene_tie")
_R7_MULTI_PP = _R7_MULTI[:1]
_R7_UNIQUE_SHARED = ("genome_unique_txome_present",)
_R7_ABSENT_OMIT = ("genome_unique_absent_omitted_exon",)
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
    # One band: the five fates, then the two mechanism hatches.
    handles = [Patch(facecolor=colour, edgecolor="#666666" if hatch else "white",
                     hatch=hatch, linewidth=0.6 if hatch else 0.4, label=label)
               for members in ROUTE7_KEY for label, colour, hatch in members]
    legend = ps.legend_below(
        axis, handles=handles, ncol=len(handles),
        fontsize=ps.FONT_ANNOTATION, handlelength=1.3, columnspacing=1.2,
        handletextpad=0.5, labelspacing=0.3, borderpad=0.0, pad_pt=2.0)
    return figure, axis, [legend]
