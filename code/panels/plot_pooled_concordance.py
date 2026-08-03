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
REQUIRED = ("sample", "spearman_all", "pearson_all")

def prepare(psite_path, footprint_path, highlight_sample=None):
    """Both per-sample tables, validated, plus the medians the panel prints."""
    sys.path.insert(0, str(HERE))
    import panel_style as ps

    frames, summary, marks = {}, {}, {}
    for name, path in (("P-site", psite_path), ("footprint", footprint_path)):
        frame = pd.read_csv(path, sep="\t")
        ps.require_columns(frame, REQUIRED, str(path))
        frames[name] = frame
        summary[name] = {column: float(np.nanmedian(frame[column]))
                         for column, _label in METRICS}
        row = frame[frame["sample"] == highlight_sample] if highlight_sample else frame.iloc[0:0]
        if not row.empty:
            marks[name] = {column: float(row[column].iloc[0]) for column, _label in METRICS}
    return {"frames": frames, "medians": summary, "marks": marks,
            "highlight_sample": highlight_sample,
            "samples": sorted(frames["P-site"]["sample"].tolist()),
            "sources": {"P-site": str(psite_path), "footprint": str(footprint_path)}}

def draw(prepared, ylim=None, figsize=(5.2, 4.4), seed=0):
    import matplotlib.pyplot as plt
    import panel_style as ps

    ps.apply_rcparams()
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
        parts = axis.violinplot(series, positions=positions, showmedians=False,
                                showextrema=False, widths=0.8)
        for body, (column, _label) in zip(parts["bodies"], METRICS):
            body.set_facecolor(colours[column][0])
            body.set_alpha(0.30)
            body.set_edgecolor("none")
        axis.boxplot(series, positions=positions, widths=0.18, showfliers=False, zorder=10,
                     medianprops=dict(color="black", lw=1.6, zorder=11),
                     boxprops=dict(color="black", zorder=10),
                     whiskerprops=dict(color="black", zorder=10),
                     capprops=dict(color="black", zorder=10))
        for position, (column, _label) in zip(positions, METRICS):
            values = frame[column].values
            axis.scatter(rng.normal(position, 0.06, len(values)), values, s=18,
                         color=colours[column][1], alpha=0.5, zorder=5, linewidths=0)
            axis.text(position, np.nanmax(values) + 0.004,
                      "%.3f" % np.nanmedian(values), va="bottom", ha="center",
                      fontsize=ps.FONT_TICK, zorder=6,
                      bbox=dict(boxstyle="round,pad=0.15", fc="white", ec="none",
                                alpha=0.75))
        mark = prepared["marks"].get(name)
        if mark:
            for position, (column, _label) in zip(positions, METRICS):
                axis.scatter([position], [mark[column]], s=34,
                             color=colours[column][1], marker="o",
                             edgecolors="white", linewidths=0.9, zorder=12)
        axis.set_xticks(positions)
        axis.set_xticklabels([label for _c, label in METRICS], fontsize=ps.FONT_TICK)
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
    parser.add_argument("--highlight-sample", default=None,
                        help="mark this sample's dot in each metric -- no text label")
    parser.add_argument("--ylim", nargs=2, type=float, default=None)
    parser.add_argument("--figsize", nargs=2, type=float, default=(5.2, 4.4))
    parser.add_argument("--output", required=True, type=Path)
    parser.add_argument("--format", dest="formats", default="pdf")
    parser.add_argument("--force", action="store_true")
    args = parser.parse_args(argv)

    sys.path.insert(0, str(HERE))
    import panel_style as ps

    prepared = prepare(args.psite, args.footprint, args.highlight_sample)
    print("[panel] %d cell lines" % len(prepared["samples"]))
    for name, medians in prepared["medians"].items():
        print("[panel]   %-10s median Spearman %.3f   Pearson %.3f"
              % (name, medians["spearman_all"], medians["pearson_all"]))

    figure, _axes = draw(prepared, tuple(args.ylim) if args.ylim else None,
                         tuple(args.figsize))
    written = ps.save(figure, args.output, ps.resolve_formats(args.formats), args.force)
    for path in written:
        print("[panel] wrote %s" % path)
    return 0

if __name__ == "__main__":
    sys.exit(main())
