#!/usr/bin/env python3
"""Figure 5 A -- distinct mapped read IDs per alignment route."""
from __future__ import annotations

import argparse
import sys
from pathlib import Path

import numpy as np

HERE = Path(__file__).resolve().parent

def prepare(taxonomy_path, samples_csv=None):
    sys.path.insert(0, str(HERE))
    import fig05_common as common

    frame = common.load_taxonomy(taxonomy_path).sort_values("delta_reads")
    labels = common.load_labels(samples_csv)
    return {"frame": frame,
            "order": frame["sample"].tolist(),
            "labels": [labels.get(s, s) for s in frame["sample"]],
            "genome_m": (frame["genome_present"] / 1e6).to_numpy(float),
            "txome_m": (frame["txome_present"] / 1e6).to_numpy(float),
            "source": str(taxonomy_path)}

def draw(prepared, figsize=(4.6, 8.0)):
    import matplotlib.pyplot as plt
    from matplotlib.lines import Line2D
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
        axis.text(gi + xmax * 0.015, yi, "+%.2fM" % delta, va="center", ha="left",
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
    return figure, axis

def main(argv=None):
    parser = argparse.ArgumentParser(
        description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument("--taxonomy", required=True, type=Path)
    parser.add_argument("--samples-csv", type=Path)
    parser.add_argument("--figsize", nargs=2, type=float, default=(4.6, 8.0))
    parser.add_argument("--output", type=Path)
    parser.add_argument("--format", dest="formats", default="pdf")
    parser.add_argument("--force", action="store_true")
    args = parser.parse_args(argv)

    sys.path.insert(0, str(HERE))
    import panel_style as ps

    prepared = prepare(args.taxonomy, args.samples_csv)
    print("[panel] %d cell lines, ordered by ascending net genome excess"
          % len(prepared["order"]))
    gap = prepared["genome_m"] - prepared["txome_m"]
    print("[panel] gap median %.2fM, range [%.2f, %.2f]"
          % (np.median(gap), gap.min(), gap.max()))

    if not args.output:
        return 0
    figure, _axis = draw(prepared, tuple(args.figsize))
    written = ps.save(figure, args.output, ps.resolve_formats(args.formats), args.force)
    for path in written:
        print("[panel] wrote %s" % path)
    return 0

if __name__ == "__main__":
    sys.exit(main())
