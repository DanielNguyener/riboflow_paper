#!/usr/bin/env python3
"""Figure 4 A-D: the per-library read categories across the 24-library cohort.

One panel per process (`--panel A`), rendered at PLOS page size with 8 pt type. All four
panels share one plot box (`AXES_TOP_OFFSET_PT` / `AXES_HEIGHT_PT` from the page top) and
one row order (libraries by ascending genome - transcriptome read-ID difference), so the
24 rows align when `assemble_figures.py` places the panels side by side; panels export
with `tight=False` so the offsets survive to the file.

A: distinct mapped read IDs, genome minus transcriptome.
B: composition of the read-ID union: SH-U, SH-M, GO-U, GO-M, TO (`categories.py`).
C: protein-coding-pseudogene ties, as % of shared genome-multimapping reads.
D: overlap with omitted alternative exons, as % of genome-only unique reads.
"""
from __future__ import annotations

import argparse
import os
import sys

import numpy as np
import pandas as pd

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, HERE)
sys.path.insert(0, os.path.join(os.path.dirname(HERE), "read_categories"))
import categories  # noqa: E402
import panel_style as ps  # noqa: E402
from panel_style import die  # noqa: E402

#: shared plot box  from page top  in points
AXES_HEIGHT_PT = 240.0
AXES_TOP_OFFSET_PT = 22.0

#: abbreviation and colour for the five categories  from the one key
UNION_KEY = categories.KEY

#: taxonomy column abbreviation colour  in manuscript category order
SEGMENTS = tuple(zip(("both_genome_unique", "both_genome_multi", "genome_only_unique",
                      "genome_only_multi", "txome_only"),
                     (categories.ABBR[k] for k in categories.KEYS),
                     (categories.COLOR[k] for k in categories.KEYS)))

TAXONOMY_REQUIRED = (
    "sample", "n_universe", "n_sh_u", "n_go_u", "n_sh_m", "n_go_m", "n_to",
    "n_genome_unique", "n_genome_multi", "n_genome_absent",
    "n_txome_present", "n_txome_absent")

#: panels C and D  the master flag  its required columns  the share
SHARE_PANELS = {
    "C": {"master": "tie_master",
          "required": ("sample", "pct_cross_pp_pc", "pct_cross_pc_pp"),
          "share": lambda f: f["pct_cross_pp_pc"] + f["pct_cross_pc_pp"]},
    "D": {"master": "reach_master",
          "required": ("sample", "n_omitted_exon_overlap", "n_go_u"),
          "share": lambda f: 100.0 * f["n_omitted_exon_overlap"] / f["n_go_u"]},
}


# ── shared loading ────────────────────────────────────────────────────────────

def load_taxonomy(path):
    """The per-library category master, with the union partition named and checked.

    The five Figure 4B segments are the five cells: genome status (unique by NH == 1,
    multimapping, absent) x transcriptome presence (a primary alignment in the post-dedup
    BAM, which RiboFlow_v2 filtered at MAPQ >= 10).
    """
    frame = pd.read_csv(path, sep="\t")
    ps.require_columns(frame, TAXONOMY_REQUIRED, str(path))
    frame = frame.copy()
    frame["both_genome_unique"] = frame["n_sh_u"]
    frame["both_genome_multi"] = frame["n_sh_m"]
    frame["genome_only_unique"] = frame["n_go_u"]
    frame["genome_only_multi"] = frame["n_go_m"]
    frame["txome_only"] = frame["n_to"]

    total = (frame["both_genome_unique"] + frame["both_genome_multi"]
             + frame["genome_only_unique"] + frame["genome_only_multi"]
             + frame["txome_only"])
    if not (total == frame["n_universe"]).all():
        bad = frame.loc[total != frame["n_universe"], "sample"].tolist()
        raise SystemExit("the union partition does not sum to n_universe for: %s"
                         % ", ".join(bad))

    frame["genome_present"] = frame["n_genome_unique"] + frame["n_genome_multi"]
    frame["txome_present"] = frame["n_txome_present"]
    for column, absent in (("genome_present", "n_genome_absent"),
                           ("txome_present", "n_txome_absent")):
        if not (frame[column] == frame["n_universe"] - frame[absent]).all():
            raise SystemExit("%s disagrees with n_universe - %s" % (column, absent))
    frame["delta_reads"] = frame["genome_present"] - frame["txome_present"]
    return frame


def load_labels(samples_csv=None):
    """{sample: GSM} from the sample table; keys normalised (spaces -> underscores)."""
    return ps.gsm_map(samples_csv) if samples_csv else {}


def sample_order(frame):
    """The shared cohort ordering: cell lines by ascending `delta_reads`."""
    return frame.sort_values("delta_reads")["sample"].tolist()


def axes_fractions(page_height_in):
    """`(bottom, top)` for `subplots_adjust` that reproduces the shared box."""
    page_pt = page_height_in * 72.0
    top = 1.0 - AXES_TOP_OFFSET_PT / page_pt
    bottom = 1.0 - (AXES_TOP_OFFSET_PT + AXES_HEIGHT_PT) / page_pt
    if bottom <= 0.0:
        raise SystemExit(
            "a %.2f in page cannot hold the shared plot box (%.0f pt from the top, %.0f pt "
            "tall) and anything below it" % (page_height_in, AXES_TOP_OFFSET_PT,
                                             AXES_HEIGHT_PT))
    return bottom, top


# ── preparation per panel ─────────────────────────────────────────────────────

def prepare_routes(taxonomy_path, samples_csv=None):
    """Panel A: per-library read IDs on each route, in the shared order."""
    frame = load_taxonomy(taxonomy_path).sort_values("delta_reads")
    labels = load_labels(samples_csv)
    return {"frame": frame,
            "order": frame["sample"].tolist(),
            "labels": [labels.get(s, s) for s in frame["sample"]],
            "genome_m": (frame["genome_present"] / 1e6).to_numpy(float),
            "txome_m": (frame["txome_present"] / 1e6).to_numpy(float)}


def prepare_union(taxonomy_path, samples_csv=None):
    """Panel B: the five category percentages per library, in the shared order."""
    frame = load_taxonomy(taxonomy_path)
    order = sample_order(frame)
    frame = frame.set_index("sample").loc[order].reset_index()
    for column, _label, _colour in SEGMENTS:
        frame["pct_" + column] = 100.0 * frame[column] / frame["n_universe"]
    labels = load_labels(samples_csv)
    return {"frame": frame, "order": order,
            "labels": [labels.get(s, s) for s in order]}


def prepare_share(letter, master, taxonomy, samples_csv=None):
    """Panels C and D: one share per library, in the shared order."""
    spec = SHARE_PANELS[letter]
    frame = pd.read_csv(master, sep="\t")
    ps.require_columns(frame, spec["required"], str(master))
    frame = frame.copy()
    frame["share"] = spec["share"](frame)

    order = sample_order(load_taxonomy(taxonomy))
    lookup = frame.set_index("sample")["share"]
    values = np.array([lookup.get(s, np.nan) for s in order], dtype=float)
    labels_map = load_labels(samples_csv)
    return {"order": order, "values": values,
            "labels": [labels_map.get(s, s) for s in order]}


# ── drawing ───────────────────────────────────────────────────────────────────

def _format_count(value):
    value = int(value)
    if value >= 1_000_000:
        return "%.1fM" % (value / 1_000_000)
    if value >= 1_000:
        return "%.1fK" % (value / 1_000)
    return str(value)


def draw_side_panel(values, labels, colour, title, xlabel, figsize):
    """The horizontal-bar idiom shared by panels A, C and D."""
    import matplotlib.pyplot as plt

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
    axis.set_yticklabels([""] * len(y))
    axis.set_ylim(-0.6, len(y) - 0.4)
    axis.grid(axis="x", alpha=0.15)
    axis.set_title(title, fontsize=ps.FONT_TITLE, loc="left", fontweight="normal")
    axis.set_xlabel(xlabel, fontsize=ps.FONT_LABEL)
    figure.tight_layout()
    bottom, top = axes_fractions(figsize[1])
    figure.subplots_adjust(bottom=bottom, top=top)   # the box A B C and D all draw in
    return figure, axis


def draw_union(prepared, figsize):
    """Panel B: one stacked bar of the five categories per library, plus its total."""
    import matplotlib.pyplot as plt

    ps.apply_rcparams()
    frame = prepared["frame"]
    y = np.arange(len(frame))
    figure, axis = plt.subplots(figsize=figsize)

    left = np.zeros(len(frame))
    for column, label, colour in SEGMENTS:
        values = frame["pct_" + column].to_numpy()
        axis.barh(y, values, left=left, color=colour, label=label,
                  edgecolor="white", linewidth=0.3)
        left += values
    for yi, total in zip(y, frame["n_universe"].to_numpy()):
        axis.text(100.5, yi, _format_count(total), va="center", ha="left",
                  fontsize=ps.FONT_ANNOTATION)

    axis.set_yticks(list(y))
    axis.set_yticklabels([""] * len(y))
    axis.set_ylim(-0.6, len(y) - 0.4)
    #: 100 % of the bar plus just enough room for widest count drawn at 100.5
    axis.set_xlim(0, 109)
    axis.set_xticks([0, 20, 40, 60, 80, 100])
    axis.set_xlabel("Read IDs in union of read alignments (%)", fontsize=ps.FONT_LABEL)
    axis.grid(axis="x", alpha=0.15)
    figure.tight_layout()
    bottom, top = axes_fractions(figsize[1])
    figure.subplots_adjust(bottom=bottom, top=top)
    return figure, axis, []


# ── the four page size renders ────────────────────────────────────────────────

def render_panel(letter, width_pt, height_in, inputs, font_pt, output, formats=("pdf",),
                 force=False):
    """One cohort panel to `output.<fmt>`, final labels applied in-process."""
    ps.FONT_TITLE = ps.FONT_LABEL = ps.FONT_TICK = font_pt
    ps.FONT_ANNOTATION = ps.FONT_INSET = font_pt
    figsize = (width_pt / 72.0, height_in)
    extras = []

    def rp(key):
        if not inputs.get(key):
            die("panel %s needs --%s" % (letter, key.replace("_", "-")))
        return inputs[key]

    if letter == "A":
        prepared = prepare_routes(rp("taxonomy"), rp("samples_csv"))
        values = prepared["genome_m"] - prepared["txome_m"]
        figure, axis = draw_side_panel(
            values, prepared["labels"], ps.GENOME,
            "Difference in aligned\nread IDs", "", figsize)
        axis.set_yticklabels(prepared["labels"])
        for artist in list(axis.texts):          # per-bar values: a table, not a panel
            artist.remove()
        # largest gap is 1.286 M  final tick at 1.5 keeps every bar short of it
        axis.set_xlim(0, 1.55)
        axis.set_xticks([0.0, 0.5, 1.0, 1.5])
        axis.set_xlabel("Genome −\ntranscriptome\n(millions)", fontsize=font_pt)
        figure.subplots_adjust(left=0.47, right=0.96)

    elif letter == "B":
        prepared = prepare_union(rp("taxonomy"), rp("samples_csv"))
        figure, axis, _ = draw_union(prepared, figsize)
        axis.set_title("Composition of the\nread-ID union", fontsize=font_pt, loc="left",
                       fontweight="normal")
        axis.set_xlabel("Read IDs in\nalignment-route union (%)", fontsize=font_pt)
        axis.set_xlim(0, 128)
        axis.set_xticks([0, 20, 40, 60, 80, 100])
        from matplotlib.patches import Patch
        handles = [Patch(facecolor=colour, edgecolor="white", linewidth=0.4, label=label)
                   for label, colour in UNION_KEY]
        extras.append(ps.legend_below(axis, handles=handles, ncol=3,
                                      fontsize=font_pt, handlelength=1.0,
                                      labelspacing=0.25, columnspacing=0.8,
                                      borderpad=0.0, pad_pt=4.0))

    elif letter in ("C", "D"):
        prepared = prepare_share(letter, rp(SHARE_PANELS[letter]["master"]),
                                 rp("taxonomy"))
        if letter == "C":
            colour, title = categories.COLOR["sh_m"], "Protein-coding–\npseudogene ties"
            xlabel = "Shared genome-\nmultimapping\nreads (%)"
        else:
            colour, title = categories.COLOR["go_u"], "Overlap with omitted\nalternative exons"
            xlabel = "Uniquely mapped\ngenome-only\nreads (%)"
        figure, axis = draw_side_panel(
            prepared["values"], prepared["labels"], colour, title, xlabel, figsize)
        if letter == "C":
            axis.set_xticks(range(0, 81, 20))
            axis.set_xticks(range(10, 80, 20), minor=True)
    else:
        die("unknown cohort panel %s" % letter)

    # tight=False  row alignment depends on shared box surviving to the file
    return ps.save(figure, output, formats, force=force,
                   extra_artists=extras or None, tight=False)


def main(argv=None):
    parser = argparse.ArgumentParser(description=__doc__,
                                     formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument("--panel", required=True, choices=tuple("ABCD"))
    parser.add_argument("--taxonomy", help="the per-library category master (A, B, C, D)")
    parser.add_argument("--samples-csv", help="S1 Table samples.csv (A, B)")
    parser.add_argument("--tie-master", help="pseudogene-tie master (C)")
    parser.add_argument("--reach-master", help="omitted alternative-exon master (D)")
    parser.add_argument("--width-pt", type=float, required=True, help="page width")
    parser.add_argument("--height-in", type=float, required=True, help="page height")
    parser.add_argument("--font-pt", type=float, default=8.0)
    parser.add_argument("--output", required=True, help="path stem, no extension")
    parser.add_argument("--format", dest="formats", default="pdf")
    parser.add_argument("--force", action="store_true")
    args = parser.parse_args(argv)

    inputs = {"taxonomy": args.taxonomy, "samples_csv": args.samples_csv,
              "tie_master": args.tie_master, "reach_master": args.reach_master}
    written = render_panel(args.panel, args.width_pt, args.height_in, inputs, args.font_pt,
                           args.output, ps.resolve_formats(args.formats), args.force)
    for path in written:
        print("wrote %s" % path)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
