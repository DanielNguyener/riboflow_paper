#!/usr/bin/env python3
"""Figure 4 A -- distinct mapped read IDs per alignment route."""
from __future__ import annotations

import sys
from pathlib import Path

import numpy as np

HERE = Path(__file__).resolve().parent

def prepare(taxonomy_path, samples_csv=None):
    sys.path.insert(0, str(HERE))
    import cohort_common as common

    frame = common.load_taxonomy(taxonomy_path).sort_values("delta_reads")
    labels = common.load_labels(samples_csv)
    return {"frame": frame,
            "order": frame["sample"].tolist(),
            "labels": [labels.get(s, s) for s in frame["sample"]],
            "genome_m": (frame["genome_present"] / 1e6).to_numpy(float),
            "txome_m": (frame["txome_present"] / 1e6).to_numpy(float)}

def draw(prepared, figsize=(4.6, 8.0)):
    import matplotlib.pyplot as plt
    from matplotlib.lines import Line2D
    import cohort_common as common
    import panel_style as ps

    ps.apply_rcparams()
    genome, txome = prepared["genome_m"], prepared["txome_m"]
    y = np.arange(len(genome))
    figure, axis = plt.subplots(figsize=figsize)

    axis.hlines(y, txome, genome, color="#bbbbbb", lw=1.8, zorder=1)
    axis.scatter(txome, y, s=26, color=ps.TXOME, zorder=3, linewidths=0.5,
                 edgecolors="white")
    axis.scatter(genome, y, s=26, color=ps.GENOME, zorder=3, linewidths=0.5,
                 edgecolors="white")
    xmax = genome.max() * 1.30
    for yi, gi, delta in zip(y, genome, genome - txome):
        axis.text(gi + xmax * 0.015, yi, "\u0394%.2fM" % delta, va="center", ha="left",
                  fontsize=ps.FONT_ANNOTATION)
    axis.set_xlim(0, xmax)
    axis.set_xticks(np.arange(0, xmax, 2.0))
    axis.set_yticks(list(y))
    axis.set_yticklabels(prepared["labels"])
    axis.set_ylim(-0.6, len(y) - 0.4)
    gap = genome - txome
    q1, median, q3 = np.percentile(gap, [25, 50, 75])
    axis.set_xlabel("Distinct mapped read IDs\n(millions)\n"
                    "gap median %.2fM  IQR [%.2f, %.2f]" % (median, q1, q3),
                    fontsize=ps.FONT_LABEL)
    axis.grid(axis="x", alpha=0.15)
    axis.legend(handles=[
        Line2D([], [], marker="o", ls="none", color=ps.GENOME, label="genome"),
        Line2D([], [], marker="o", ls="none", color=ps.TXOME, label="transcriptome")],
        fontsize=ps.FONT_TICK, frameon=False, loc="lower right")
    figure.tight_layout()
    bottom, top = common.stack_axes_fractions(figsize[1])
    figure.subplots_adjust(bottom=bottom, top=top)     # the box panel B also draws in
    return figure, axis
