#!/usr/bin/env python3
"""Figure 4 B -- composition of the union of read IDs across the two routes."""
from __future__ import annotations

import sys
from pathlib import Path

import numpy as np

HERE = Path(__file__).resolve().parent

SEGMENTS = (
    ("both_genome_unique", "SH-U", "#a6d96a"),
    ("both_genome_multi", "SH-M", "#1a7d1a"),
    ("genome_only_unique", "GO-U", "#7fb9da"),
    ("genome_only_multi", "GO-M", "#0d57a1"),
    ("txome_only", "TO", "#cc3d3d"))
NOT_IN_UNION = "#dddddd"

def _format_count(value):
    value = int(value)
    if value >= 1_000_000:
        return "%.1fM" % (value / 1_000_000)
    if value >= 1_000:
        return "%.1fK" % (value / 1_000)
    return str(value)

def prepare(taxonomy_path, samples_csv=None):
    sys.path.insert(0, str(HERE))
    import cohort_common as common

    frame = common.load_taxonomy(taxonomy_path)
    order = common.sample_order(frame)
    provenance = "shared cohort order, derived from %s" % taxonomy_path
    frame = frame.set_index("sample").loc[order].reset_index()
    for column, _label, _colour in SEGMENTS:
        frame["pct_" + column] = 100.0 * frame[column] / frame["n_universe"]
    labels = common.load_labels(samples_csv)
    return {"frame": frame, "order": order, "order_provenance": provenance,
            "labels": [labels.get(s, s) for s in order]}

def draw(prepared, figsize=(6.4, 8.6), show_labels=True, show_key=True):
    import matplotlib.pyplot as plt
    from matplotlib.patches import Rectangle
    import cohort_common as common
    import panel_style as ps

    ps.apply_rcparams()
    frame = prepared["frame"]
    y = np.arange(len(frame))
    figure, axis = plt.subplots(figsize=figsize)

    left = np.zeros(len(frame))
    for column, label, colour in SEGMENTS:
        values = frame["pct_" + column].to_numpy()
        axis.barh(y, values, left=left, color=colour, label=label,
                  edgecolor="white", linewidth=0.3)
        left += values
    for yi, total in zip(y, frame["n_universe"].to_numpy()):
        axis.text(100.5, yi, _format_count(total), va="center", ha="left",
                  fontsize=ps.FONT_ANNOTATION)

    axis.set_yticks(list(y))
    axis.set_yticklabels(prepared["labels"] if show_labels else [""] * len(y))
    axis.set_ylim(-0.6, len(y) - 0.4)
    #: 100 % of the bar, plus just enough room for the widest count drawn at 100.5.
    axis.set_xlim(0, 109)
    axis.set_xticks([0, 20, 40, 60, 80, 100])
    axis.set_xlabel("Read IDs in union of read alignments (%)", fontsize=ps.FONT_LABEL)
    axis.grid(axis="x", alpha=0.15)

    cells = {(0, 0): SEGMENTS[0][2], (0, 1): SEGMENTS[2][2],
             (1, 0): SEGMENTS[1][2], (1, 1): SEGMENTS[3][2],
             (2, 0): SEGMENTS[4][2], (2, 1): NOT_IN_UNION}
    # The key hangs below the shared cohort_common box. Its anchor is a FRACTION of the axes
    # height, so at small sizes it lands on the x label -- draw small with the key off.
    if not show_key:
        figure.tight_layout()
        bottom, top = common.stack_axes_fractions(figsize[1])
        figure.subplots_adjust(bottom=bottom, top=top)
        return figure, axis, []
    inset = axis.inset_axes([0.36, -0.185, 0.28, 0.085])
    inset.set_xlim(0, 2.0)
    inset.set_ylim(0, 3.5)
    inset.axis("off")
    for row, name in enumerate(("unique", "multimapping", "absent")):
        y0 = 3.0 - (row + 1)
        for column in range(2):
            inset.add_patch(Rectangle((column, y0), 1.0, 1.0, facecolor=cells[(row, column)],
                                      edgecolor="white", linewidth=1.0))
        inset.text(-0.12, y0 + 0.5, name, ha="right", va="center", fontsize=ps.FONT_INSET)
    for column, name in enumerate(("present", "absent")):
        inset.text(column + 0.5, 3.05, name, ha="center", va="bottom", fontsize=ps.FONT_INSET)
    inset.text(1.0, 3.50, "TRANSCRIPTOME", ha="center", va="bottom", fontsize=ps.FONT_INSET)
    inset.text(-0.72, 0.5, "GENOME", transform=inset.transAxes, ha="center",
               va="center", rotation=90, fontsize=ps.FONT_INSET)
    figure.tight_layout()
    bottom, top = common.stack_axes_fractions(figsize[1])
    figure.subplots_adjust(bottom=bottom, top=top)     # the box panel A also draws in
    return figure, axis, [inset]
