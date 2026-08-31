#!/usr/bin/env python3
"""Per-transcript genome-vs-transcriptome concordance, distribution per cell line."""
from __future__ import annotations

import argparse
import sys
from pathlib import Path

import numpy as np
import pandas as pd

HERE = Path(__file__).resolve().parent
REQUIRED = ("sample", "transcript_id", "spearman", "pearson")

def prepare(psite_path, footprint_path, samples_csv=None, highlight=None):
    """Load both tables, order the cell lines, and resolve the highlighted transcripts.

    Ordering (ascending median P-site Spearman) is applied once and shared by both sub-panels.
    """
    sys.path.insert(0, str(HERE))
    import panel_style as ps

    frames = {}
    for name, path in (("P-site", psite_path), ("footprint", footprint_path)):
        frame = pd.read_csv(path, sep="\t")
        ps.require_columns(frame, REQUIRED, str(path))
        frames[name] = frame

    order = (frames["P-site"].groupby("sample")["spearman"].median()
             .sort_values().index.tolist())

    labels = {}
    if samples_csv:
        gsm = pd.read_csv(samples_csv)
        ps.require_columns(gsm, ("cell_line", "ribo_GSM"), str(samples_csv))
        labels = {str(c).replace(" ", "_"): str(g)
                  for c, g in zip(gsm["cell_line"], gsm["ribo_GSM"])}

    marks = {}
    for name, frame in frames.items():
        marks[name] = {}
        for key, transcript_id in (highlight or {}).items():
            if key == "sample":
                continue
            row = frame[(frame["sample"] == highlight.get("sample"))
                        & (frame["transcript_id"] == transcript_id)]
            if not row.empty:
                marks[name][key] = {"transcript_id": transcript_id,
                                    "spearman": float(row["spearman"].iloc[0]),
                                    "pearson": float(row["pearson"].iloc[0])}
    return {"frames": frames, "order": order, "labels": labels, "marks": marks,
            "highlight_sample": (highlight or {}).get("sample"),
            "n_rows": {k: int(len(v)) for k, v in frames.items()},
            "sources": {"P-site": str(psite_path), "footprint": str(footprint_path)}}

def _half_box(axis, data, positions, fill, line):
    boxes = axis.boxplot(data, positions=positions, widths=0.4, showfliers=False,
                         showcaps=False, patch_artist=True,
                         medianprops=dict(color=line, lw=1.8),
                         boxprops=dict(edgecolor=line, lw=1.1),
                         whiskerprops=dict(color=line, lw=1.1))
    for box in boxes["boxes"]:
        box.set_facecolor(fill)
        box.set_alpha(0.6)

#: The leader from a marker to its gene name: a thin plain line, no head, no text box.
LEADER = dict(arrowstyle="-", color="#555", lw=0.6, shrinkA=0, shrinkB=2)

def draw(prepared, ylim=(0.0, 1.0), figsize=(11.0, 4.4), layout="side"):
    """`layout="side"`: P-site | footprint side by side; `"stacked"`: shared x axis,
    GSM labels appear once (the form that fits beside panel D)."""
    import matplotlib.pyplot as plt
    import matplotlib.patches as mpatches
    import panel_style as ps

    ps.apply_rcparams()
    order = prepared["order"]
    positions = np.arange(1, len(order) + 1)
    if layout == "stacked":
        figure, axes = plt.subplots(2, 1, figsize=figsize, sharex=True, sharey=True)
    else:
        figure, axes = plt.subplots(1, 2, figsize=figsize, sharey=True)

    for axis, (name, frame) in zip(axes, prepared["frames"].items()):
        spearman = [frame.loc[frame["sample"] == s, "spearman"].dropna().values
                    for s in order]
        pearson = [frame.loc[frame["sample"] == s, "pearson"].dropna().values
                   for s in order]
        _half_box(axis, spearman, positions - 0.2, ps.SPEARMAN_FILL, ps.SPEARMAN_LINE)
        _half_box(axis, pearson, positions + 0.2, ps.PEARSON_FILL, ps.PEARSON_LINE)

        marks = prepared["marks"].get(name, {})
        sample = prepared["highlight_sample"]
        if sample in order and marks:
            index = order.index(sample) + 1
            low, high = ylim
            for key, entry in marks.items():
                dot = dict(s=34, marker="o", edgecolors="white", linewidths=0.9,
                           zorder=12)
                diamond = dict(s=22, marker="D", edgecolors="black", linewidths=0.8,
                               zorder=13)
                # GAPDH is a diamond in EVERY sub-panel, on or off the axis, so the
                # discordant example reads as the same transcript in both views.
                style = diamond if key == "gapdh" else dot
                if entry["spearman"] < low:
                    floor = low + 0.012 * (high - low)
                    axis.scatter([index - 0.2], [floor], color=ps.SPEARMAN_LINE,
                                 clip_on=False, **diamond)
                    axis.scatter([index + 0.2], [floor], color=ps.PEARSON_LINE,
                                 clip_on=False, **diamond)
                    # Plain text on a thin leader, up-right into the empty band; no box.
                    axis.annotate(key.upper(), xy=(index + 0.2, floor), xytext=(16, 18),
                                  textcoords="offset points", ha="left", va="bottom",
                                  fontsize=ps.FONT_TICK, zorder=14,
                                  arrowprops=LEADER)
                else:
                    axis.scatter([index - 0.2], [entry["spearman"]],
                                 color=ps.SPEARMAN_LINE, **style)
                    axis.scatter([index + 0.2], [entry["pearson"]],
                                 color=ps.PEARSON_LINE, **style)
                    near_top = entry["spearman"] > high - 0.08 * (high - low)
                    # Down-right when the dot sits near the top; otherwise up-right.
                    xytext = (16, -22) if near_top else (16, 14)
                    va = "top" if near_top else "bottom"
                    axis.annotate(key.upper(), xy=(index + 0.2, entry["pearson"]),
                                  xytext=xytext, textcoords="offset points",
                                  ha="left", va=va, fontsize=ps.FONT_TICK, zorder=14,
                                  arrowprops=LEADER)

        axis.set_xticks(positions)
        axis.set_xticklabels([prepared["labels"].get(s, s) for s in order],
                             rotation=90, fontsize=ps.FONT_TICK)
        axis.set_xlim(0.4, len(order) + 0.6)
        axis.grid(axis="y", alpha=0.15)
        if layout == "stacked":
            # Sub-panel name inside the axes, bottom left; a title would cost a line each.
            axis.text(0.01, 0.04, name, transform=axis.transAxes, ha="left", va="bottom",
                      fontsize=ps.FONT_TITLE, zorder=14,
                      bbox=dict(boxstyle="round,pad=0.2", fc="white", ec="#555", lw=0.6))
        else:
            axis.set_title(name, fontsize=ps.FONT_TITLE)

    if layout == "stacked":
        axes[0].tick_params(labelbottom=False)
        figure.supylabel("per-transcript correlation", fontsize=ps.FONT_LABEL, x=0.012)
    else:
        axes[0].set_ylabel("per-transcript correlation", fontsize=ps.FONT_LABEL)
    axes[0].set_ylim(*ylim)
    axes[-1].legend(handles=[
        mpatches.Patch(facecolor=ps.SPEARMAN_FILL, alpha=0.6, edgecolor=ps.SPEARMAN_LINE,
                       lw=1.1, label="Spearman $\\rho$"),
        mpatches.Patch(facecolor=ps.PEARSON_FILL, alpha=0.6, edgecolor=ps.PEARSON_LINE,
                       lw=1.1, label="Pearson $r$")],
        loc="lower right", fontsize=ps.FONT_TICK, frameon=True)
    figure.tight_layout()
    if layout == "stacked":
        figure.subplots_adjust(hspace=0.08, left=0.13)
        for loc, label in zip(axes[0].get_yticks(), axes[0].get_yticklabels()):
            if loc <= ylim[0] + 1e-9:
                label.set_visible(False)
    return figure, axes

def main(argv=None):
    parser = argparse.ArgumentParser(
        description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument("--psite", required=True, type=Path)
    parser.add_argument("--footprint", required=True, type=Path)
    parser.add_argument("--samples-csv", type=Path,
                        help="the sample table; supplies the GSM axis labels")
    parser.add_argument("--highlight-sample", default="HeLa")
    parser.add_argument("--highlight-gapdh", default="ENST00000396861.5")
    parser.add_argument("--highlight-comt", default="ENST00000361682.11")
    parser.add_argument("--ylim", nargs=2, type=float, default=(0.0, 1.0))
    parser.add_argument("--figsize", nargs=2, type=float, default=(11.0, 4.4))
    parser.add_argument("--layout", choices=("side", "stacked"), default="side",
                        help="side: P-site | footprint; stacked: P-site over footprint, "
                             "one shared x axis")
    parser.add_argument("--output", required=True, type=Path)
    parser.add_argument("--format", dest="formats", default="pdf")
    parser.add_argument("--force", action="store_true")
    args = parser.parse_args(argv)

    sys.path.insert(0, str(HERE))
    import panel_style as ps

    highlight = {"sample": args.highlight_sample,
                 "gapdh": args.highlight_gapdh, "comt": args.highlight_comt}
    prepared = prepare(args.psite, args.footprint, args.samples_csv, highlight)
    print("[panel] %d cell lines; P-site %d rows, footprint %d rows"
          % (len(prepared["order"]), prepared["n_rows"]["P-site"],
             prepared["n_rows"]["footprint"]))
    print("[panel] order (ascending median Spearman): %s"
          % ", ".join(prepared["order"][:4] + ["..."] + prepared["order"][-2:]))

    figure, _axes = draw(prepared, tuple(args.ylim), tuple(args.figsize), args.layout)
    written = ps.save(figure, args.output, ps.resolve_formats(args.formats), args.force)
    for path in written:
        print("[panel] wrote %s" % path)
    return 0

if __name__ == "__main__":
    sys.exit(main())
