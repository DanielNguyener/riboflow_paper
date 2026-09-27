#!/usr/bin/env python3
"""Genome-versus-transcriptome per-base coverage for any transcript, from a coverage HDF5."""
from __future__ import annotations

import argparse
import json
import sys
from pathlib import Path

import numpy as np

HERE = Path(__file__).resolve().parent
REPO = HERE.parents[1]
COVERAGE_DIR = REPO / "code" / "coverage"

BUILD_HINT = """\
Build one with:

    python code/coverage/build_shared_coverage.py \\
        --sample SAMPLE \\
        --genome-bam        .../SAMPLE.post_dedup.bam \\
        --transcriptome-bam .../SAMPLE.transcriptome.post_dedup.bam \\
        --gtf     /path/to/gencode.gtf.gz \\
        --appris  /path/to/appris_transcript_lengths.tsv \\
        --qc-genome results/ribo_seq_qc/genome/tables/readlen_window_qc.csv \\
        --qc-txome  results/ribo_seq_qc/transcriptome/tables/readlen_window_qc.csv \\
        --output results/coverage

or for a whole cohort, code/coverage/build_cohort_coverage.py. This program does not
build anything itself: see docs/hdf5_schema.md."""

def _import_coverage_modules():
    for directory in (str(COVERAGE_DIR), str(HERE)):
        if directory not in sys.path:
            sys.path.insert(0, directory)
    import coverage_schema
    import compute_coverage_concordance as metrics
    return coverage_schema, metrics

def check_coverage_file(path, expect_sample=None, require_regions=True):
    """Validate the input before drawing anything. Returns the file's identity dict.

    Raises SystemExit naming everything that is wrong, plus `BUILD_HINT`.
    """
    coverage_schema, _metrics = _import_coverage_modules()

    path = Path(path)
    if not path.exists():
        raise SystemExit("coverage file does not exist: %s\n\n%s" % (path, BUILD_HINT))
    problems = coverage_schema.validate_file(path)
    if problems:
        raise SystemExit(
            "%s is not a usable coverage file:\n%s\n\n%s"
            % (path, "\n".join("  - %s" % p for p in problems), BUILD_HINT))

    with coverage_schema.open_coverage(path) as coverage:
        identity = coverage.identity()
        complaints = []
        if not identity["sample"]:
            complaints.append("it declares no sample")
        if expect_sample and identity["sample"] != expect_sample:
            complaints.append("it is sample %r, but %r was expected"
                              % (identity["sample"], expect_sample))
        if identity["coordinate_system"] != coverage_schema.COORDINATE_SYSTEM:
            complaints.append("coordinate_system is %r, expected %r"
                              % (identity["coordinate_system"],
                                 coverage_schema.COORDINATE_SYSTEM))
        missing_routes = set(coverage_schema.ROUTES) - set(identity["routes"])
        if missing_routes:
            complaints.append("it does not declare route(s): %s"
                              % ", ".join(sorted(missing_routes)))
        if not identity["psite_placement"]:
            complaints.append("it records no P-site placement rule")
        if require_regions and not (coverage.cds_start >= 0).any():
            complaints.append("it carries no CDS bounds to draw regions from")
    if complaints:
        raise SystemExit("%s cannot be plotted:\n%s\n\n%s"
                         % (path, "\n".join("  - %s" % c for c in complaints), BUILD_HINT))
    return identity

def load_tracks(coverage_path, gene_id):
    """Everything a plot needs for one transcript, and nothing about how to draw it.

    The published framing is fixed: the CDS window at the file's own trim, the raw
    counts, the canonical region overlay when the file carries regions.
    """
    coverage_schema, _metrics = _import_coverage_modules()

    with coverage_schema.open_coverage(coverage_path) as coverage:
        index = coverage.resolve_gene(gene_id)

        info = coverage.transcript_info(index)
        regions = coverage.regions_of(index)
        effective_trim = coverage.trim

        start, end = coverage.slice_region(index, "CDS", trim=effective_trim)
        if end <= start:
            raise SystemExit(
                "%s has a CDS of %d nt, which does not survive a %d nt trim at each "
                "end." % (info["transcript_id"], info["cds_len"], effective_trim))
        # CDS relative axis  a trimmed window runs [trim, cds_len - trim)
        axis_origin = regions["CDS"][0]

        tracks = {}
        for key, signal in (("g_ps", "genome_psite"), ("t_ps", "txome_psite"),
                            ("g_fp", "genome_footprint"), ("t_fp", "txome_footprint")):
            tracks[key] = coverage.get_track(index, signal)[start:end]

        counts = coverage.event_counts(index)
        states = {}
        for key, signal in (("g_ps", "genome_psite"), ("t_ps", "txome_psite"),
                            ("g_fp", "genome_footprint"), ("t_fp", "txome_footprint")):
            states[key] = coverage_schema.describe_coverage_state(
                counts[signal], int(tracks[key].sum()))
        sample = coverage.sample

    return {
        "sample": sample,
        "transcript_id": info["transcript_id"],
        "gene_id": info["gene_id"],
        "gene_name": info["gene_name"],
        "transcript_len": info["transcript_len"],
        "cds_len": info["cds_len"],
        "region": "cds",
        "trim": effective_trim,
        "x_start": start - axis_origin,
        "x_end": end - axis_origin,
        "axis_origin": axis_origin,
        "slice": [start, end],
        "x": np.arange(start - axis_origin, end - axis_origin),
        "raw": tracks,
        "values": tracks,
        "states": states,
        "regions": regions,
        "overlay": "canonical" if regions else "none",
        "requested_gene_id": gene_id,
    }

def overlay_intervals(tracks):
    """[(display label, start, end)] on the plotted axis, for the chosen overlay."""
    origin = tracks["axis_origin"]
    if tracks["overlay"] == "canonical":
        rows = [(label, start, end)
                for label, (start, end) in sorted(tracks["regions"].items(),
                                                  key=lambda kv: kv[1])]
    else:
        return []
    return [(label, start - origin, end - origin) for label, start, end in rows]

def annotate_correlations(tracks):
    """Spearman rho and the documented log-Pearson, from the RAW integer counts.

    Uses `_spear`/`_pe_log2`, never `scipy.stats.pearsonr` (differs at ~1e-15).
    """
    _coverage_schema, metrics = _import_coverage_modules()
    raw = tracks["raw"]
    return {
        "psite": {"spearman": metrics._spear(raw["g_ps"], raw["t_ps"]),
                  "pearson": metrics._pe_log2(raw["g_ps"], raw["t_ps"])},
        "footprint": {"spearman": metrics._spear(raw["g_fp"], raw["t_fp"]),
                      "pearson": metrics._pe_log2(raw["g_fp"], raw["t_fp"])},
    }

def _draw_track(axis, x, genome, txome, style, states, ylabel, mirrored):
    """One signal's axes. Mirrored bars for sparse P-sites, filled steps for depth."""
    import panel_style as ps

    if mirrored:
        axis.bar(x, genome, width=1.0, color=ps.GENOME, linewidth=0, zorder=1)
        axis.bar(x, -np.asarray(txome), width=1.0, color=ps.TXOME, linewidth=0, zorder=1)
        axis.axhline(0, color="black", lw=0.6, zorder=3)
        peak = max(float(np.max(genome)) if len(genome) else 0,
                   float(np.max(txome)) if len(txome) else 0, 1.0)
        axis.set_ylim(-peak * 1.08, peak * 1.08)
        from matplotlib.ticker import FuncFormatter, MaxNLocator
        axis.yaxis.set_major_locator(MaxNLocator(integer=True))
        axis.yaxis.set_major_formatter(FuncFormatter(lambda v, _p: "%g" % abs(v)))
    else:
        axis.fill_between(x, genome, step="mid", color=ps.GENOME, alpha=0.45,
                          linewidth=0, zorder=1)
        axis.fill_between(x, txome, step="mid", color=ps.TXOME, alpha=0.45,
                          linewidth=0, zorder=2)
        axis.step(x, genome, where="mid", color=ps.GENOME, lw=1.2, zorder=3)
        axis.step(x, txome, where="mid", color=ps.TXOME, lw=1.2, zorder=3)

    axis.set_ylabel(ylabel, fontsize=style["label"])
    axis.margins(x=0.005)
    axis.grid(axis="y", alpha=0.15)

    notes = []
    for key, name in ((states[0], "genome"), (states[1], "transcriptome")):
        if key == "no_reads_assigned":
            notes.append("no %s reads assigned to this transcript" % name)
        elif key == "reads_outside_requested_slice":
            notes.append("%s reads present, none in this window" % name)
    if notes:
        axis.axhspan(*axis.get_ylim(), facecolor="none", edgecolor=ps.MISSING_HATCH,
                     hatch="///", linewidth=0, zorder=0)
        axis.text(0.5, 0.5, "\n".join(notes), transform=axis.transAxes,
                  ha="center", va="center", fontsize=style["annotation"], color="#555",
                  bbox=dict(boxstyle="round,pad=0.35", fc="white", ec="#999", lw=0.8),
                  zorder=8)

def plot_coverage(tracks, correlations=None, figsize=None, title=None,
                  labels="full", title_correlations=False, route_legend=False):
    """Draw one transcript's coverage. Returns (figure, axes).

    `labels="minimal"` drops the in-axes text; boundary lines stay, numbers go to the record.
    `route_legend` adds an unframed genome/transcriptome colour key on the x-label line,
    flush right (returned as `figure._route_legend`; hand it to `save(..., extra_artists=[...])`).
    """
    import matplotlib.pyplot as plt
    import panel_style as ps

    ps.apply_rcparams()
    style = {"label": ps.FONT_LABEL, "title": ps.FONT_TITLE, "tick": ps.FONT_TICK,
             "annotation": ps.FONT_ANNOTATION}
    wanted = ["psite", "footprint"]
    figsize = figsize or (10.0, 2.0 * len(wanted) + 0.8)
    figure, axes = plt.subplots(len(wanted), 1, figsize=figsize, sharex=True, squeeze=False)
    axes = [a[0] for a in axes]

    x = tracks["x"]
    for axis, which in zip(axes, wanted):
        if which == "psite":
            _draw_track(axis, x, tracks["values"]["g_ps"], tracks["values"]["t_ps"], style,
                        (tracks["states"]["g_ps"], tracks["states"]["t_ps"]),
                        "P-site\ncoverage", mirrored=True)
        else:
            _draw_track(axis, x, tracks["values"]["g_fp"], tracks["values"]["t_fp"], style,
                        (tracks["states"]["g_fp"], tracks["states"]["t_fp"]),
                        "footprint\ncoverage", mirrored=False)
        if correlations and labels == "full":
            entry = correlations[which]
            axis.text(0.99, 0.94,
                      "Spearman $\\rho$ = %.3f\nPearson $r$ = %.3f"
                      % (entry["spearman"], entry["pearson"]),
                      transform=axis.transAxes, ha="right", va="top",
                      fontsize=style["annotation"],
                      bbox=dict(boxstyle="round,pad=0.35", fc="white", ec="#555", lw=0.9),
                      zorder=9)

    for boundary, text in ((tracks["x_start"], "start +%d nt" % tracks["trim"]),
                           (tracks["x_end"], "stop −%d nt" % tracks["trim"])):
        for axis in axes:
            axis.axvline(boundary, color="#555", ls="--", lw=1.0, zorder=5)
        if labels == "minimal":
            continue
        axes[0].annotate(text, xy=(boundary, 1.0), xycoords=("data", "axes fraction"),
                         xytext=(0, 2), textcoords="offset points", ha="center",
                         va="bottom", fontsize=style["tick"], color="#555")

    label_box = dict(boxstyle="round,pad=0.35", fc="white", ec="#555", lw=0.9)
    if labels == "full":
        axes[0].text(0.012, 0.90, "genome", transform=axes[0].transAxes, ha="left",
                     va="top", fontsize=style["title"], bbox=dict(label_box), zorder=9)
    if wanted[0] == "psite" and labels == "full":
        axes[0].text(0.012, 0.10, "transcriptome", transform=axes[0].transAxes,
                     ha="left", va="bottom", fontsize=style["title"],
                     bbox=dict(label_box), zorder=9)

    axes[-1].set_xlabel("CDS position", fontsize=style["label"])
    heading = title if title is not None else "%s (%s) - %s" % (
        tracks["gene_name"], tracks["transcript_id"], tracks["sample"])
    if title_correlations and correlations:
        # second smaller title line carrying what the in axes box would have said
        names = {"psite": "P-site", "footprint": "footprint"}
        stats = "; ".join("%s $\\rho$ = %.3f, $r$ = %.3f"
                          % (names[w], correlations[w]["spearman"],
                             correlations[w]["pearson"]) for w in wanted)
        axes[0].set_title(heading, fontsize=style["title"], pad=14)
        axes[0].annotate(stats, xy=(0.5, 1.0), xycoords="axes fraction", xytext=(0, 3),
                         textcoords="offset points", ha="center", va="bottom",
                         fontsize=style["annotation"], color="#333")
    else:
        axes[0].set_title(heading, fontsize=style["title"], pad=14)
    figure.tight_layout()
    figure._route_legend = None
    if route_legend:
        # colour key not a data annotation  drawn even with minimal labels
        # shares the x label line  vertically centred on the label  measured after
        # tight_layout which moves the axes  right aligned to the axes right edge
        from matplotlib.patches import Patch
        axis = axes[-1]
        figure.canvas.draw()
        label_box = axis.xaxis.label.get_window_extent(figure.canvas.get_renderer())
        y_mid = float(axis.transAxes.inverted().transform(
            (0.0, (label_box.y0 + label_box.y1) / 2.0))[1])
        figure._route_legend = axis.legend(
            [Patch(facecolor=ps.GENOME, edgecolor="none"),
             Patch(facecolor=ps.TXOME, edgecolor="none")],
            ["genome", "transcriptome"], loc="center right", bbox_to_anchor=(1.0, y_mid),
            bbox_transform=axis.transAxes, ncol=2, frameon=False,
            fontsize=style["annotation"], handlelength=1.2, handleheight=0.9,
            columnspacing=1.5, borderaxespad=0.0, borderpad=0.0)
    return figure, axes

def render(argv=None):
    """Draw the panel and RETURN the render record (the tests assert on it directly)."""
    parser = argparse.ArgumentParser(
        description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument("--coverage-h5", required=True, type=Path, dest="coverage")
    parser.add_argument("--gene-id", required=True)
    parser.add_argument("--title")
    parser.add_argument("--labels", choices=("full", "minimal"), default="full")
    parser.add_argument("--title-correlations", action="store_true")
    parser.add_argument("--record-json", type=Path)
    parser.add_argument("--figsize", nargs=2, type=float, metavar=("W", "H"))
    parser.add_argument("--output", required=True, type=Path)
    parser.add_argument("--format", dest="formats", default="pdf")
    parser.add_argument("--force", action="store_true")
    args = parser.parse_args(argv)

    sys.path.insert(0, str(HERE))
    import panel_style as ps

    identity = check_coverage_file(args.coverage)
    print("[panel] %s: sample %s, assay %s, routes %s, P-site %s, schema v%d"
          % (args.coverage.name, identity["sample"], identity["assay"],
             "+".join(identity["routes"]), identity["psite_placement"],
             identity["schema_version"]))

    tracks = load_tracks(args.coverage, args.gene_id)
    print("[panel] resolved %s -> %s (%s)"
          % (args.gene_id, tracks["transcript_id"], tracks["gene_name"]))
    if tracks["overlay"] != "none":
        print("[panel] region overlay (%s): %s"
              % (tracks["overlay"],
                 ", ".join("%s %d-%d" % row for row in overlay_intervals(tracks))))
    print("[panel] %s window [%d, %d) on the %s axis (transcript slice [%d, %d) of %d nt)"
          % (tracks["region"], tracks["x_start"], tracks["x_end"],
             "CDS-relative" if tracks["region"] == "cds" else "transcript",
             tracks["slice"][0], tracks["slice"][1], tracks["transcript_len"]))
    for key, state in sorted(tracks["states"].items()):
        if state != "covered":
            print("[panel]   %s: %s" % (key, state))

    correlations = annotate_correlations(tracks)
    figure, _axes = plot_coverage(tracks, correlations,
                                  tuple(args.figsize) if args.figsize else None,
                                  args.title, args.labels, args.title_correlations,
                                  route_legend=True)
    written = ps.save(figure, args.output, ps.resolve_formats(args.formats), args.force,
                      extra_artists=[figure._route_legend] if figure._route_legend else None)
    record = {
        "generator": "code/panels/plot_transcript_coverage.py",
        "coverage_file": args.coverage.name,
        "coverage_identity": identity,
        "sample": tracks["sample"],
        "requested": {"gene_id": args.gene_id, "transcript_id": None},
        "resolved": {"transcript_id": tracks["transcript_id"],
                     "gene_id": tracks["gene_id"], "gene_name": tracks["gene_name"]},
        "region": tracks["region"], "trim": tracks["trim"],
        "region_overlay": tracks["overlay"],
        "overlay_intervals": [list(row) for row in overlay_intervals(tracks)],
        "axis_window": [tracks["x_start"], tracks["x_end"]],
        "transcript_slice": tracks["slice"],
        "signal": "both", "normalize": "none",
        "coverage_states": tracks["states"],
        "correlations": correlations,
        "labels": args.labels,
        "route_legend": True,
        "outputs": [str(p) for p in written],
    }
    if args.record_json:
        import json
        args.record_json.parent.mkdir(parents=True, exist_ok=True)
        args.record_json.write_text(json.dumps(record, indent=2, default=str))
    for path in written:
        print("[panel] wrote %s" % path)
    return record

def main(argv=None):
    render(argv)
    return 0

if __name__ == "__main__":
    sys.exit(main())
