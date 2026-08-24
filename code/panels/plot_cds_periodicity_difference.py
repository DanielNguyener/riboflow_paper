#!/usr/bin/env python3
"""Figure 2 B -- whole-CDS frame-0 periodicity, genome minus transcriptome."""
from __future__ import annotations

import argparse
import sys
from pathlib import Path

import numpy as np

HERE = Path(__file__).resolve().parent
#: Colour-scale limit, in percentage points; +/-2 fits the measured span without clipping.
#: Display only -- `prepare()` never sees it. `--vlim` overrides.
VLIM = 2.0
#: Fraction of VLIM past which in-cell text flips to white; relative so it tracks VLIM.
WHITE_TEXT_FRACTION = 0.75

def prepare(frame_genome, frame_txome, samples_csv):
    sys.path.insert(0, str(HERE))
    import _fig02_common as common

    samples, lengths, genome, txome = common.load_pair(
        frame_genome, frame_txome, "pct_frame0")
    difference = genome - txome
    finite = difference[~np.isnan(difference)]
    return {"samples": samples, "lengths": lengths, "difference": difference,
            "gsm": common.gsm_labels(samples_csv, samples),
            "n_shared_cells": int(finite.size),
            "median": float(np.median(finite)),
            "mean_abs": float(np.abs(finite).mean()),
            "max_abs": float(np.abs(finite).max()),
            "sources": {"frame_genome": str(frame_genome), "frame_txome": str(frame_txome),
                        "samples_csv": str(samples_csv)}}

def draw(prepared, axes_size=None, margins=None, vlim=VLIM, type_scale="large",
         show_ylabels=True):
    import matplotlib.pyplot as plt
    from matplotlib.colors import LinearSegmentedColormap
    sys.path.insert(0, str(HERE))
    import _fig02_common as common
    import panel_style as ps

    ps.apply_rcparams()
    sizes = common.grid_type(type_scale)
    samples, lengths = prepared["samples"], prepared["lengths"]
    difference = prepared["difference"]
    colormap = LinearSegmentedColormap.from_list(
        "genome_txome", ["#d62728", "white", "#2ca25f"])
    colormap.set_bad("#eeeeee")
    norm = plt.Normalize(vmin=-vlim, vmax=vlim)

    # Shared with fig02A (`_fig02_common.MARGINS`); right margin = colourbar gutter.
    width, height = axes_size or common.AXES_SIZE
    left, bottom, right, top = margins or common.MARGINS
    fig_w, fig_h = left + width + right, bottom + height + top
    figure, axis = plt.subplots(figsize=(fig_w, fig_h))
    axis.set_position([left / fig_w, bottom / fig_h, width / fig_w, height / fig_h])
    common.draw_grid(axis, difference, colormap, norm)
    white_above = WHITE_TEXT_FRACTION * vlim
    for i in range(difference.shape[0]):
        for j in range(difference.shape[1]):
            value = difference[i, j]
            if not np.isnan(value):
                axis.text(common.cell_centre(j), common.cell_centre(i), "%+.1f" % value,
                          ha="center", va="center", fontsize=sizes["annotation"],
                          color="white" if abs(value) > white_above else "black")
    common.style_grid(axis, samples, lengths, prepared["gsm"], sizes, show_ylabels)

    mappable = plt.cm.ScalarMappable(cmap=colormap, norm=norm)
    cax = figure.add_axes([(left + width + 0.15) / fig_w, (bottom + height * 0.15) / fig_h,
                           0.15 / fig_w, height * 0.70 / fig_h])
    bar = figure.colorbar(mappable, cax=cax, extend="both")
    # matplotlib rasterises the colorbar's solids by default -- force vector.
    if bar.solids is not None:
        bar.solids.set_rasterized(False)
    # Short on purpose; the sign convention is the caption's job.
    bar.set_label("Δ CDS periodicity (%)", fontsize=sizes["label"])
    bar.ax.tick_params(labelsize=sizes["tick"])
    return figure, axis

def main(argv=None):
    parser = argparse.ArgumentParser(
        description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument("--frame-genome", required=True, type=Path)
    parser.add_argument("--frame-txome", required=True, type=Path)
    parser.add_argument("--samples-csv", required=True, type=Path)
    parser.add_argument("--vlim", type=float, default=VLIM)
    parser.add_argument("--output", required=True, type=Path)
    parser.add_argument("--format", dest="formats", default="pdf")
    parser.add_argument("--force", action="store_true")
    parser.add_argument("--axes-size", nargs=2, type=float, metavar=("W", "H"),
                        help="grid size in inches (default: _fig02_common.AXES_SIZE)")
    parser.add_argument("--margins", nargs=4, type=float, metavar=("L", "B", "R", "T"),
                        help="margins in inches (default: _fig02_common.MARGINS)")
    parser.add_argument("--type-scale", choices=("large", "base"), default="large",
                        help="large: standalone panel type; base: journal-page type (8-12 pt)")
    parser.add_argument("--hide-ylabels", action="store_true",
                        help="drop the GSM row labels (when placed beside panel A, which has them)")
    args = parser.parse_args(argv)

    sys.path.insert(0, str(HERE))
    import panel_style as ps

    prepared = prepare(args.frame_genome, args.frame_txome, args.samples_csv)
    print("[panel] %d shared cells; median %+.2f, mean|Δ| %.2f, max|Δ| %.2f"
          % (prepared["n_shared_cells"], prepared["median"], prepared["mean_abs"],
             prepared["max_abs"]))

    figure, _axis = draw(prepared, tuple(args.axes_size) if args.axes_size else None,
                         tuple(args.margins) if args.margins else None, args.vlim,
                         args.type_scale, not args.hide_ylabels)
    written = ps.save(figure, args.output, ps.resolve_formats(args.formats), args.force,
                      tight=False)
    for path in written:
        print("[panel] wrote %s" % path)
    return 0

if __name__ == "__main__":
    sys.exit(main())
