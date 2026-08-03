#!/usr/bin/env python3
"""The horizontal-bar idiom shared by Figure 5 panels C and D."""
from __future__ import annotations

import sys
from pathlib import Path

import numpy as np

HERE = Path(__file__).resolve().parent

def draw_side_panel(values, labels, colour, title, xlabel, figsize=(3.0, 8.0),
                    show_labels=False):
    import matplotlib.pyplot as plt
    sys.path.insert(0, str(HERE))
    import fig05_common as common
    import panel_style as ps

    ps.apply_rcparams()
    values = np.asarray(values, dtype=float)
    y = np.arange(len(values))
    figure, axis = plt.subplots(figsize=figsize)
    axis.barh(y, values, color=colour, edgecolor="white", linewidth=0.3)

    finite = values[np.isfinite(values)]
    xmax = (finite.max() * 1.30) if finite.size else 1.0
    for yi, value in zip(y, values):
        if np.isfinite(value):
            axis.text(value + xmax * 0.03, yi, "%.1f" % value, va="center", ha="left",
                      fontsize=ps.FONT_ANNOTATION)
    axis.set_xlim(0, xmax)
    axis.set_yticks(list(y))
    axis.set_yticklabels(labels if show_labels else [""] * len(y))
    axis.set_ylim(-0.6, len(y) - 0.4)
    axis.grid(axis="x", alpha=0.15)
    axis.set_title(title, fontsize=ps.FONT_TITLE, loc="left", fontweight="normal")

    subtext = common.median_iqr(values)
    axis.set_xlabel(xlabel + ("\n" + subtext.replace("\n", "  ") if subtext else ""),
                    fontsize=ps.FONT_LABEL)
    figure.tight_layout()
    return figure, axis
