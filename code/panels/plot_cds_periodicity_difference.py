#!/usr/bin/env python3
"""Figure 2 B -- whole-CDS frame-0 periodicity, genome minus transcriptome."""
from __future__ import annotations

import argparse
import sys
from pathlib import Path

import numpy as np

HERE = Path(__file__).resolve().parent
VLIM = 4.0

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

def draw(prepared, axes_size=(7.0, 7.8), vlim=VLIM):
    import matplotlib.pyplot as plt
    from matplotlib.colors import LinearSegmentedColormap
    sys.path.insert(0, str(HERE))
    import _fig02_common as common
    import panel_style as ps

    ps.apply_rcparams()
    samples, lengths = prepared["samples"], prepared["lengths"]
    difference = prepared["difference"]
    colormap = LinearSegmentedColormap.from_list(
        "genome_txome", ["#d62728", "white", "#2ca25f"])
    colormap.set_bad("#eeeeee")
    norm = plt.Normalize(vmin=-vlim, vmax=vlim)

    width, height = axes_size
    left, bottom, right, top = 1.25, 0.35, 1.25, 0.15
    fig_w, fig_h = left + width + right, bottom + height + top
    figure, axis = plt.subplots(figsize=(fig_w, fig_h))
    axis.set_position([left / fig_w, bottom / fig_h, width / fig_w, height / fig_h])
    common.draw_grid(axis, difference, colormap, norm)
    for i in range(difference.shape[0]):
        for j in range(difference.shape[1]):
            value = difference[i, j]
            if not np.isnan(value):
                axis.text(common.cell_centre(j), common.cell_centre(i), "%+.1f" % value,
                          ha="center", va="center", fontsize=ps.FONT_ANNOTATION,
                          color="white" if abs(value) > 9 else "black")
    common.style_grid(axis, samples, lengths, prepared["gsm"])

    mappable = plt.cm.ScalarMappable(cmap=colormap, norm=norm)
    cax = figure.add_axes([(left + width + 0.15) / fig_w, (bottom + height * 0.15) / fig_h,
                           0.15 / fig_w, height * 0.70 / fig_h])
    bar = figure.colorbar(mappable, cax=cax, extend="both")
    # matplotlib rasterises the colorbar's solids by default, which shows as a blemished
    # gradient in the exported PDF. Force it to vector.
    if bar.solids is not None:
        bar.solids.set_rasterized(False)
    bar.set_label("Δ whole-CDS frame-0 %  (genome − transcriptome)",
                  fontsize=ps.FONT_LABEL)
    bar.ax.tick_params(labelsize=ps.FONT_TICK)
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
    args = parser.parse_args(argv)

    sys.path.insert(0, str(HERE))
    import panel_style as ps

    prepared = prepare(args.frame_genome, args.frame_txome, args.samples_csv)
    print("[panel] %d shared cells; median %+.2f, mean|Δ| %.2f, max|Δ| %.2f"
          % (prepared["n_shared_cells"], prepared["median"], prepared["mean_abs"],
             prepared["max_abs"]))

    figure, _axis = draw(prepared, vlim=args.vlim)
    written = ps.save(figure, args.output, ps.resolve_formats(args.formats), args.force)
    for path in written:
        print("[panel] wrote %s" % path)
    return 0

if __name__ == "__main__":
    sys.exit(main())
