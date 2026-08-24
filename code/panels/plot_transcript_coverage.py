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

SIGNAL_CHOICES = ("psite", "footprint", "both")
NORMALIZE_CHOICES = ("none", "per-million", "max")
OVERLAY_CHOICES = ("auto", "canonical", "none")

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

def load_tracks(coverage_path, gene_id=None, transcript_id=None, region="whole",
                trim=None, normalize="none", overlay="auto"):
    """Everything a plot needs for one transcript, and nothing about how to draw it."""
    coverage_schema, _metrics = _import_coverage_modules()

    with coverage_schema.open_coverage(coverage_path) as coverage:
        if gene_id:
            index = coverage.resolve_gene(gene_id, transcript_id=transcript_id)
        elif transcript_id:
            index = coverage.index_of_transcript(transcript_id)
        else:
            raise SystemExit("give --gene-id and/or --transcript-id")

        info = coverage.transcript_info(index)
        regions = coverage.regions_of(index)
        file_trim = coverage.trim
        effective_trim = file_trim if trim is None else trim

        if region == "cds":
            start, end = coverage.slice_region(index, "CDS", trim=effective_trim)
            if end <= start:
                raise SystemExit(
                    "%s has a CDS of %d nt, which does not survive a %d nt trim at each "
                    "end. Use --region whole or a smaller --trim."
                    % (info["transcript_id"], info["cds_len"], effective_trim))
            # CDS-relative axis: a trimmed window runs [trim, cds_len - trim).
            axis_origin = regions["CDS"][0]
        elif region == "whole":
            start, end = 0, info["transcript_len"]
            axis_origin = 0
        else:
            raise SystemExit("unknown --region %r" % region)

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

    scale = 1.0
    if normalize == "per-million":
        total = sum(int(t.sum()) for t in tracks.values())
        scale = 1e6 / total if total else 1.0
    plotted = {}
    for key, values in tracks.items():
        if normalize == "max":
            peak = int(values.max()) if values.size else 0
            plotted[key] = values / peak if peak else values.astype(float)
        elif normalize == "per-million":
            plotted[key] = values * scale
        else:
            plotted[key] = values

    return {
        "sample": sample,
        "transcript_id": info["transcript_id"],
        "gene_id": info["gene_id"],
        "gene_name": info["gene_name"],
        "transcript_len": info["transcript_len"],
        "cds_len": info["cds_len"],
        "region": region,
        "trim": effective_trim,
        "file_trim": file_trim,
        "x_start": start - axis_origin,
        "x_end": end - axis_origin,
        "axis_origin": axis_origin,
        "slice": [start, end],
        "x": np.arange(start - axis_origin, end - axis_origin),
        "raw": tracks,
        "values": plotted,
        "states": states,
        "regions": regions,
        "overlay": resolve_overlay(overlay, regions),
        "normalize": normalize,
        "requested_gene_id": gene_id,
        "requested_transcript_id": transcript_id,
    }

def resolve_overlay(requested, regions):
    """Which region overlay to draw; explicit `canonical` without regions is an error."""
    if requested == "none":
        return "none"
    if requested == "auto":
        return "canonical" if regions else "none"
    if requested == "canonical":
        if not regions:
            raise SystemExit("--regions canonical was asked for, but this file carries "
                             "no CDS bounds for this transcript.")
        return "canonical"
    raise SystemExit("unknown --regions %r" % requested)

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

REGION_SHADING = {
    "UTR5": ("#dfe6ee", "UTR5"),
    "CDS": ("#ffffff", "CDS"),
    "UTR3": ("#dfe6ee", "UTR3"),
}

def _draw_region_overlay(axes, tracks, style):
    """Shade each region across every axes and label it once, along the top."""
    intervals = overlay_intervals(tracks)
    if not intervals:
        return
    span = max(end for _l, _s, end in intervals) - min(s for _l, s, _e in intervals)
    for label, start, end in intervals:
        if end <= start:
            continue
        colour, text = REGION_SHADING.get(label, ("#eeeeee", label))
        for axis in axes:
            if colour != "#ffffff":
                axis.axvspan(start, end, facecolor=colour, edgecolor="none",
                             alpha=0.55, zorder=0)
            axis.axvline(start, color="#7a7a7a", ls="--", lw=0.8, zorder=5)
        if span and (end - start) / span >= 0.04:
            axes[0].annotate(text, xy=((start + end) / 2.0, 1.0),
                             xycoords=("data", "axes fraction"), xytext=(0, 2),
                             textcoords="offset points", ha="center", va="bottom",
                             fontsize=style["tick"], color="#444")
    axes[0].annotate("regions: %s" % tracks["overlay"], xy=(1.0, 1.0),
                     xycoords="axes fraction", xytext=(0, 12),
                     textcoords="offset points", ha="right", va="bottom",
                     fontsize=style["annotation"], color="#888")

def plot_coverage(tracks, signal="both", correlations=None, figsize=None, title=None,
                  labels="full", title_correlations=False):
    """Draw one transcript's coverage. Returns (figure, axes).

    `labels="minimal"` drops the in-axes text; boundary lines stay, numbers go to the record.
    """
    import matplotlib.pyplot as plt
    import panel_style as ps

    ps.apply_rcparams()
    style = {"label": ps.FONT_LABEL, "title": ps.FONT_TITLE, "tick": ps.FONT_TICK,
             "annotation": ps.FONT_ANNOTATION}
    wanted = ["psite", "footprint"] if signal == "both" else [signal]
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

    if tracks["region"] == "whole":
        _draw_region_overlay(axes, tracks, style)
    else:
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

    unit = {"none": "", "per-million": ", per million",
            "max": ", scaled to max"}[tracks["normalize"]]
    axes[-1].set_xlabel(("CDS position" if tracks["region"] == "cds"
                         else "transcript position") + unit, fontsize=style["label"])
    heading = title if title is not None else "%s (%s) - %s" % (
        tracks["gene_name"], tracks["transcript_id"], tracks["sample"])
    if title_correlations and correlations:
        # A second, smaller title line carrying what the in-axes box would have said.
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
    return figure, axes

def render(argv=None):
    """Draw the panel and RETURN the render record (the tests assert on it directly)."""
    parser = argparse.ArgumentParser(
        description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument("--coverage-h5", required=True, type=Path, dest="coverage",
                        help="a shared-coordinate coverage HDF5 -- the only data input")
    parser.add_argument("--expect-sample",
                        help="fail unless the file declares this sample")
    parser.add_argument("--gene-id", help="versioned or unversioned")
    parser.add_argument("--transcript-id", help="versioned or unversioned")
    parser.add_argument("--signal", default="both", choices=SIGNAL_CHOICES)
    parser.add_argument("--region", default="whole", choices=("whole", "cds"),
                        help="which window to plot")
    parser.add_argument("--regions", default="auto", choices=OVERLAY_CHOICES,
                        dest="overlay",
                        help="which region overlay to draw on the whole-transcript view: "
                             "the canonical UTR5/CDS/UTR3, or none")
    parser.add_argument("--trim", type=int, default=None,
                        help="override the file's own paper_cds_trim")
    parser.add_argument("--normalize", default="none", choices=NORMALIZE_CHOICES)
    parser.add_argument("--annotate-correlation", action="store_true")
    parser.add_argument("--title")
    parser.add_argument("--labels", choices=("full", "minimal"), default="full",
                        help="minimal: no route names, correlation box or boundary "
                             "captions inside the axes (the numbers stay in the record)")
    parser.add_argument("--title-correlations", action="store_true",
                        help="add a second title line with each track's rho and r "
                             "(needs --annotate-correlation)")
    parser.add_argument("--record-json", type=Path,
                        help="also write the render record (resolved transcript, window, "
                             "correlations) to this JSON file")
    parser.add_argument("--figsize", nargs=2, type=float, metavar=("W", "H"))
    parser.add_argument("--output", required=True, type=Path)
    parser.add_argument("--format", dest="formats", default="pdf")
    parser.add_argument("--force", action="store_true")
    args = parser.parse_args(argv)

    if not args.gene_id and not args.transcript_id:
        raise SystemExit("give --gene-id and/or --transcript-id")

    sys.path.insert(0, str(HERE))
    import panel_style as ps

    identity = check_coverage_file(args.coverage, args.expect_sample)
    print("[panel] %s: sample %s, assay %s, routes %s, P-site %s, schema v%d"
          % (args.coverage.name, identity["sample"], identity["assay"],
             "+".join(identity["routes"]), identity["psite_placement"],
             identity["schema_version"]))

    tracks = load_tracks(args.coverage, args.gene_id, args.transcript_id,
                         args.region, args.trim, args.normalize, args.overlay)
    print("[panel] resolved %s -> %s (%s)"
          % (args.gene_id or args.transcript_id, tracks["transcript_id"],
             tracks["gene_name"]))
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

    correlations = annotate_correlations(tracks) if args.annotate_correlation else None
    figure, _axes = plot_coverage(tracks, args.signal, correlations,
                                  tuple(args.figsize) if args.figsize else None,
                                  args.title, args.labels, args.title_correlations)
    written = ps.save(figure, args.output, ps.resolve_formats(args.formats), args.force)
    record = {
        "generator": "code/panels/plot_transcript_coverage.py",
        "coverage_file": args.coverage.name,
        "coverage_identity": identity,
        "sample": tracks["sample"],
        "requested": {"gene_id": args.gene_id, "transcript_id": args.transcript_id},
        "resolved": {"transcript_id": tracks["transcript_id"],
                     "gene_id": tracks["gene_id"], "gene_name": tracks["gene_name"]},
        "region": tracks["region"], "trim": tracks["trim"],
        "region_overlay": tracks["overlay"],
        "overlay_intervals": [list(row) for row in overlay_intervals(tracks)],
        "axis_window": [tracks["x_start"], tracks["x_end"]],
        "transcript_slice": tracks["slice"],
        "signal": args.signal, "normalize": args.normalize,
        "coverage_states": tracks["states"],
        "correlations": correlations,
        "labels": args.labels,
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
