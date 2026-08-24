#!/usr/bin/env python3
"""Pooled per-cell-line genome-vs-transcriptome concordance."""
from __future__ import annotations

import argparse
import sys
from pathlib import Path

import numpy as np
import pandas as pd

HERE = Path(__file__).resolve().parent
METRICS = (("spearman_all", "Spearman $\\rho$"), ("pearson_all", "Pearson $r$"))
#: The same two metrics with one-glyph tick labels, for a narrow slot beside panel C.
SHORT_LABELS = {"spearman_all": "$\\rho$", "pearson_all": "$r$"}
REQUIRED = ("sample", "spearman_all", "pearson_all")
#: The box is drawn over the dots, so its face has to let them through.
BOX_FACE_ALPHA = 0.35
DOT_ALPHA = 0.55

def prepare(psite_path, footprint_path):
    """Both per-sample tables, validated, plus the medians the panel prints.

    No cell line is singled out; the panel's claim is about the distribution.
    """
    sys.path.insert(0, str(HERE))
    import panel_style as ps

    frames, summary = {}, {}
    for name, path in (("P-site", psite_path), ("footprint", footprint_path)):
        frame = pd.read_csv(path, sep="\t")
        ps.require_columns(frame, REQUIRED, str(path))
        frames[name] = frame
        summary[name] = {column: float(np.nanmedian(frame[column]))
                         for column, _label in METRICS}
    return {"frames": frames, "medians": summary,
            "samples": sorted(frames["P-site"]["sample"].tolist()),
            "sources": {"P-site": str(psite_path), "footprint": str(footprint_path)}}

def draw(prepared, ylim=None, figsize=(5.2, 4.4), seed=0, short_labels=False,
         points=True):
    import matplotlib.pyplot as plt
    from matplotlib.colors import to_rgba
    import panel_style as ps

    ps.apply_rcparams()
    # (fill, line) per metric: box face = tint, box strokes and dots = dark shade.
    colours = {"spearman_all": (ps.SPEARMAN_FILL, ps.SPEARMAN_LINE),
               "pearson_all": (ps.PEARSON_FILL, ps.PEARSON_LINE)}

    if ylim is None:
        lowest = min(float(frame[[c for c, _ in METRICS]].min().min())
                     for frame in prepared["frames"].values())
        ylim = (np.floor((lowest - 0.01) * 20) / 20, 1.005)

    rng = np.random.RandomState(seed)
    figure, axes = plt.subplots(1, 2, figsize=figsize, sharey=True)
    positions = np.arange(1, len(METRICS) + 1)

    for axis, (name, frame) in zip(axes, prepared["frames"].items()):
        series = [frame[column].dropna().values for column, _label in METRICS]
        # Box and dots only (a violin over n=24 overstates the density). Whiskers/caps come
        # back two per box, hence the pairwise slice. The box sits ON TOP of the dots, so
        # translucency goes on the FACE COLOUR, not artist `alpha` (which would fade strokes).
        box = axis.boxplot(series, positions=positions, widths=0.42, showfliers=False,
                           patch_artist=True, zorder=12)
        for k, (column, _label) in enumerate(METRICS):
            fill, line = colours[column]
            box["boxes"][k].set(facecolor=to_rgba(fill, BOX_FACE_ALPHA), edgecolor=line,
                                linewidth=1.4, zorder=12)
            box["medians"][k].set(color=line, linewidth=2.0, zorder=13)
            for stroke in box["whiskers"][2 * k:2 * k + 2] + box["caps"][2 * k:2 * k + 2]:
                stroke.set(color=line, linewidth=1.4, zorder=12)

        for position, (column, _label) in zip(positions, METRICS):
            values = frame[column].values
            # UNDER the box, faint. `facecolor=`, not `color=`, or the black edge is lost.
            if points:
                axis.scatter(rng.normal(position, 0.06, len(values)), values, s=18,
                             facecolor=colours[column][1], edgecolors="black",
                             linewidths=0.35, alpha=DOT_ALPHA, zorder=5)
            axis.text(position, np.nanmax(values) + 0.004,
                      "%.3f" % np.nanmedian(values), va="bottom", ha="center",
                      fontsize=ps.FONT_TICK, zorder=14,
                      bbox=dict(boxstyle="round,pad=0.15", fc="white", ec="none",
                                alpha=0.75))
        axis.set_xticks(positions)
        axis.set_xticklabels([SHORT_LABELS[c] if short_labels else label
                              for c, label in METRICS], fontsize=ps.FONT_TICK)
        axis.grid(axis="y", alpha=0.15)
        axis.set_title(name, fontsize=ps.FONT_TITLE)
        axis.set_ylim(*ylim)

    from matplotlib.ticker import FormatStrFormatter, MultipleLocator
    axes[0].set_ylabel("correlation", fontsize=ps.FONT_LABEL)
    axes[0].yaxis.set_major_locator(MultipleLocator(0.05))
    axes[0].yaxis.set_major_formatter(FormatStrFormatter("%.2f"))
    figure.tight_layout()
    return figure, axes

def main(argv=None):
    parser = argparse.ArgumentParser(
        description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument("--psite", required=True, type=Path)
    parser.add_argument("--footprint", required=True, type=Path)
    parser.add_argument("--ylim", nargs=2, type=float, default=None)
    parser.add_argument("--figsize", nargs=2, type=float, default=(5.2, 4.4))
    parser.add_argument("--no-points", action="store_true",
                        help="box only, without the per-cell-line dots (narrow slots)")
    parser.add_argument("--short-labels", action="store_true",
                        help="tick labels rho / r instead of Spearman rho / Pearson r")
    parser.add_argument("--output", required=True, type=Path)
    parser.add_argument("--format", dest="formats", default="pdf")
    parser.add_argument("--force", action="store_true")
    args = parser.parse_args(argv)

    sys.path.insert(0, str(HERE))
    import panel_style as ps

    prepared = prepare(args.psite, args.footprint)
    print("[panel] %d cell lines" % len(prepared["samples"]))
    for name, medians in prepared["medians"].items():
        print("[panel]   %-10s median Spearman %.3f   Pearson %.3f"
              % (name, medians["spearman_all"], medians["pearson_all"]))

    figure, _axes = draw(prepared, tuple(args.ylim) if args.ylim else None,
                         tuple(args.figsize), short_labels=args.short_labels, points=not args.no_points)
    written = ps.save(figure, args.output, ps.resolve_formats(args.formats), args.force)
    for path in written:
        print("[panel] wrote %s" % path)
    return 0

if __name__ == "__main__":
    sys.exit(main())
