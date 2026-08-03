#!/usr/bin/env python3
"""Figure 5 E -- where a transcript's transcriptome reads went in the genome route."""
from __future__ import annotations

import argparse
import sys
from pathlib import Path

import numpy as np
import pandas as pd

HERE = Path(__file__).resolve().parent
REQUIRED = ("sample", "transcript_id", "gene_name", "category", "n_reads",
            "pct_of_txome_assigned")

sys.path.insert(0, str(HERE.parent / "alignment_fate"))
from transcript_fate_lib import PANEL_TRANSCRIPTS
SEGMENTS = (
    ("genome_unique", "genome-unique (kept)", "#a6d96a", "black"),
    ("genome_multi_pseudogene_tie", "genome-multi: pseudogene tie (excluded)",
     "#1a7d1a", "white"),
    ("other", "other (genome-multi w/o tie, or unaligned)", "#dddddd", "black"))

def prepare(fates_path, sample=None):
    """The panel's two transcripts, in order, from the maintained fates master table.

    Every transcript in `PANEL_TRANSCRIPTS` must be present: a missing one is a broken
    input, not a smaller figure.
    """
    sys.path.insert(0, str(HERE))
    import panel_style as ps

    frame = pd.read_csv(fates_path, sep="\t")
    ps.require_columns(frame, REQUIRED, str(fates_path))
    if sample:
        frame = frame[frame["sample"] == sample]
        if frame.empty:
            raise SystemExit("no rows for sample %r in %s" % (sample, fates_path))

    available = list(frame["transcript_id"])
    ordered, absent = [], []
    for transcript_id, gene_name in PANEL_TRANSCRIPTS:
        base = transcript_id.split(".", 1)[0]
        hits = sorted(set(t for t in available
                          if t == transcript_id or t.split(".", 1)[0] == base))
        if not hits:
            absent.append("%s (%s)" % (transcript_id, gene_name))
            continue
        ordered.extend(h for h in hits if h not in ordered)
    if absent:
        raise SystemExit(
            "%s does not contain %d of Figure 5E's transcript(s): %s\nThis panel draws "
            "GAPDH and COMT; a table without them is the wrong table."
            % (fates_path, len(absent), ", ".join(absent)))

    entries = []
    for transcript_id in ordered:
        rows = frame[frame["transcript_id"] == transcript_id]
        counts = dict(zip(rows["category"], rows["n_reads"]))
        missing = [c for c, _l, _x, _t in SEGMENTS if c not in counts]
        if missing:
            raise SystemExit("%s is missing categor(y/ies): %s"
                             % (transcript_id, ", ".join(missing)))
        total = int(sum(counts[c] for c, _l, _x, _t in SEGMENTS))
        percentages = dict(zip(rows["category"], rows["pct_of_txome_assigned"]))
        if abs(sum(percentages.values()) - 100.0) > 1e-6:
            raise SystemExit("%s: the categories sum to %.6f %%, not 100"
                             % (transcript_id, sum(percentages.values())))
        entries.append({
            "transcript_id": transcript_id,
            "gene_name": rows["gene_name"].iloc[0],
            "sample": rows["sample"].iloc[0],
            "counts": {c: int(counts[c]) for c, _l, _x, _t in SEGMENTS},
            "pct": {c: float(percentages[c]) for c, _l, _x, _t in SEGMENTS},
            "n_txome_assigned": total,
        })
    return {"entries": entries, "source": str(fates_path)}

def draw(prepared, title=None, figsize=None, label_threshold=6.0):
    import matplotlib.pyplot as plt
    from matplotlib.patches import Patch
    import panel_style as ps

    ps.apply_rcparams()
    entries = prepared["entries"]
    figsize = figsize or (3.4, 1.1 * len(entries) + 1.6)
    figure, axis = plt.subplots(figsize=figsize)
    y = np.arange(len(entries))[::-1]

    for yi, entry in zip(y, entries):
        left = 0.0
        for category, _label, colour, text_colour in SEGMENTS:
            width = entry["pct"][category]
            axis.barh(yi, width, left=left, color=colour, edgecolor="white",
                      linewidth=0.6, height=0.55)
            if width >= label_threshold:
                axis.text(left + width / 2, yi, "%.0f%%" % width, va="center",
                          ha="center", fontsize=ps.FONT_ANNOTATION, color=text_colour)
            left += width

    axis.set_yticks(list(y))
    axis.set_yticklabels([e["gene_name"] or e["transcript_id"] for e in entries],
                         fontsize=ps.FONT_TICK)
    axis.set_ylim(-0.6, len(entries) - 0.4)
    axis.set_xlim(0, 100)
    axis.set_xlabel("% of transcriptome-\nassigned reads", fontsize=ps.FONT_LABEL)
    axis.grid(axis="x", alpha=0.15)
    if title is None:
        title = entries[0]["sample"] if entries else ""
    axis.set_title(title, fontsize=ps.FONT_TITLE, loc="left", fontweight="normal")

    figure.tight_layout()
    legend = ps.legend_below(
        axis,
        handles=[Patch(color=colour, label=label) for _c, label, colour, _t in SEGMENTS],
        fontsize=ps.FONT_ANNOTATION - 1, handlelength=1.1)
    return figure, axis, [legend]

def main(argv=None):
    parser = argparse.ArgumentParser(
        description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument("--fates", required=True, type=Path)
    parser.add_argument("--sample", help="restrict a multi-sample table to one sample")
    parser.add_argument("--title")
    parser.add_argument("--figsize", nargs=2, type=float)
    parser.add_argument("--output", required=True, type=Path)
    parser.add_argument("--format", dest="formats", default="pdf")
    parser.add_argument("--force", action="store_true")
    args = parser.parse_args(argv)

    sys.path.insert(0, str(HERE))
    import panel_style as ps

    prepared = prepare(args.fates, args.sample)
    for entry in prepared["entries"]:
        print("[panel] %-20s %-8s assigned %6d  %s"
              % (entry["transcript_id"], entry["gene_name"], entry["n_txome_assigned"],
                 "  ".join("%s %.1f%%" % (c.split("_")[-1], entry["pct"][c])
                           for c, _l, _x, _t in SEGMENTS)))
        assert sum(entry["counts"].values()) == entry["n_txome_assigned"]

    figure, _axis, extra = draw(prepared, args.title,
                                tuple(args.figsize) if args.figsize else None)
    written = ps.save(figure, args.output, ps.resolve_formats(args.formats), args.force,
                      extra_artists=extra)
    for path in written:
        print("[panel] wrote %s" % path)
    return 0

if __name__ == "__main__":
    sys.exit(main())
