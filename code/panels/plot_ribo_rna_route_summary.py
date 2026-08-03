#!/usr/bin/env python3
"""Figure 4 C -- pooled Ribo-seq versus RNA-seq correlation, by alignment route."""
from __future__ import annotations

import argparse
import sys
from pathlib import Path

import numpy as np
import pandas as pd

HERE = Path(__file__).resolve().parent
ROUTES = ("genome", "transcriptome")
ROUTE_COLOUR = {"genome": "#3a923a", "transcriptome": "#cc3d3d"}
METRIC_COLUMN = {"spearman": "spearman_rho", "pearson": "pearson_log2_raw"}
METRIC_LABEL = {"spearman": "Spearman $\\rho$, Ribo-seq vs RNA-seq",
                "pearson": "Pearson $r$, Ribo-seq vs RNA-seq"}
COUNT_COLUMNS = {"genome": ("genome_ribo_reads", "genome_rna_reads"),
                 "transcriptome": ("txome_ribo_reads", "txome_rna_reads")}

def prepare(route_master, metric="spearman", check_counts=None, check_sample=None,
            tolerance=5e-3):
    sys.path.insert(0, str(HERE))
    import panel_style as ps

    column = METRIC_COLUMN[metric]
    frame = pd.read_csv(route_master, sep="\t")
    ps.require_columns(frame, ("sample", "route", column), str(route_master))
    wide = frame.pivot(index="sample", columns="route", values=column)
    missing = [r for r in ROUTES if r not in wide.columns]
    if missing:
        raise SystemExit("%s has no rows for route(s): %s"
                         % (route_master, ", ".join(missing)))

    checks = []
    if check_counts and check_sample:
        from scipy import stats
        counts = pd.read_csv(check_counts, sep="\t")
        for route in ROUTES:
            ribo_column, rna_column = COUNT_COLUMNS[route]
            ribo = counts[ribo_column].to_numpy(float)
            rna = counts[rna_column].to_numpy(float)
            recomputed = (float(stats.spearmanr(ribo, rna).correlation) if metric == "spearman"
                          else float(stats.pearsonr(np.log2(ribo + 1), np.log2(rna + 1))[0]))
            published = float(wide.loc[check_sample, route])
            delta = abs(recomputed - published)
            checks.append({"route": route, "recomputed": recomputed,
                           "master": published, "delta": delta, "ok": delta < tolerance})
            if delta >= tolerance:
                raise SystemExit(
                    "%s %s: the per-transcript counts give %.4f but the route master says "
                    "%.4f (difference %.4f >= %.4f). The two panels of this figure would "
                    "disagree." % (check_sample, route, recomputed, published, delta,
                                   tolerance))

    return {"wide": wide, "metric": metric, "column": column,
            "medians": {r: float(np.nanmedian(wide[r])) for r in ROUTES},
            "samples": wide.index.tolist(), "checks": checks,
            "source": str(route_master)}

def draw(prepared, ylim=None, highlight=None, figsize=(4.6, 4.6), seed=0):
    import matplotlib.pyplot as plt
    import panel_style as ps

    ps.apply_rcparams()
    wide = prepared["wide"]
    if ylim is None:
        lowest = min(float(wide[r].min()) for r in ROUTES)
        ylim = (np.floor((lowest - 0.01) * 20) / 20, 1.005)

    figure, axis = plt.subplots(figsize=figsize)
    positions = np.arange(1, len(ROUTES) + 1)
    series = [wide[r].dropna().to_numpy() for r in ROUTES]
    parts = axis.violinplot(series, positions=positions, showmedians=False,
                            showextrema=False, widths=0.8)
    for body, route in zip(parts["bodies"], ROUTES):
        body.set_facecolor(ROUTE_COLOUR[route])
        body.set_alpha(0.30)
        body.set_edgecolor("none")
    axis.boxplot(series, positions=positions, widths=0.18, showfliers=False, zorder=10,
                 medianprops=dict(color="black", lw=1.6, zorder=11),
                 boxprops=dict(color="black", zorder=10),
                 whiskerprops=dict(color="black", zorder=10),
                 capprops=dict(color="black", zorder=10))

    jitter = np.random.RandomState(seed)
    for position, route in zip(positions, ROUTES):
        values = wide[route].to_numpy()
        x = jitter.normal(position, 0.06, len(values))
        axis.scatter(x, values, s=22, color=ROUTE_COLOUR[route], alpha=0.85, zorder=5,
                     linewidths=0)
        if highlight and highlight in wide.index:
            k = wide.index.get_loc(highlight)
            axis.scatter([x[k]], [values[k]], s=42, facecolors="white",
                         edgecolors=ROUTE_COLOUR[route], linewidths=1.6, zorder=11)
        axis.text(position, np.nanmax(values) + 0.004, "%.3f" % np.nanmedian(values),
                  va="bottom", ha="center", fontsize=ps.FONT_TICK, zorder=6,
                  bbox=dict(boxstyle="round,pad=0.15", fc="white", ec="none", alpha=0.75))

    axis.set_xticks(positions)
    axis.set_xticklabels(list(ROUTES), fontsize=ps.FONT_TICK)
    axis.set_ylim(*ylim)
    axis.set_ylabel(METRIC_LABEL[prepared["metric"]], fontsize=ps.FONT_LABEL)
    axis.grid(axis="y", alpha=0.15)
    axis.set_box_aspect(1.0)
    figure.tight_layout()
    return figure, axis

def main(argv=None):
    parser = argparse.ArgumentParser(
        description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument("--route-master", required=True, type=Path)
    parser.add_argument("--metric", choices=("spearman", "pearson"), default="spearman")
    parser.add_argument("--highlight", default="HeLa")
    parser.add_argument("--check-counts", type=Path,
                        help="per-transcript counts to cross-check the highlighted sample")
    parser.add_argument("--ylim", nargs=2, type=float, default=(0.8, 1.0))
    parser.add_argument("--figsize", nargs=2, type=float, default=(4.6, 4.6))
    parser.add_argument("--output", required=True, type=Path)
    parser.add_argument("--format", dest="formats", default="pdf")
    parser.add_argument("--force", action="store_true")
    args = parser.parse_args(argv)

    sys.path.insert(0, str(HERE))
    import panel_style as ps

    prepared = prepare(args.route_master, args.metric, args.check_counts, args.highlight)
    print("[panel] %d cell lines; median genome %.3f, transcriptome %.3f"
          % (len(prepared["samples"]), prepared["medians"]["genome"],
             prepared["medians"]["transcriptome"]))
    for check in prepared["checks"]:
        print("[panel]   cross-check %s: scatter %.4f vs master %.4f (Δ %.5f) %s"
              % (check["route"], check["recomputed"], check["master"], check["delta"],
                 "OK" if check["ok"] else "FAILED"))

    figure, _axis = draw(prepared, tuple(args.ylim), args.highlight, tuple(args.figsize))
    written = ps.save(figure, args.output, ps.resolve_formats(args.formats), args.force)
    for path in written:
        print("[panel] wrote %s" % path)
    return 0

if __name__ == "__main__":
    sys.exit(main())
