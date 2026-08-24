#!/usr/bin/env python3
"""Shared loading for the five Figure-5 panels: the cohort ordering and its labels."""
from __future__ import annotations

import sys
from pathlib import Path

import pandas as pd

HERE = Path(__file__).resolve().parent

TAXONOMY_REQUIRED = (
    "sample", "n_universe", "n_gU_tU", "n_gU_tM", "n_gM_tU", "n_gM_tM",
    "n_gU_tA", "n_gM_tA", "n_gA_tU", "n_gA_tM",
    "n_genome_unique", "n_genome_multi", "n_genome_absent",
    "n_txome_unique", "n_txome_multi", "n_txome_absent")

#: The plot box shared by Figure 5's cohort panels A-D: distance from the PAGE TOP plus a
#: HEIGHT, both in points. Both must match across panels or the 24 rows misalign at
#: assembly; panels export with `tight=False` so the offsets survive to the file.
STACK_PAGE_HEIGHT_IN = 8.0
STACK_AXES_TOP_OFFSET_PT = 40.0     # page top -> top of the plot box (C and D's title band)
STACK_AXES_HEIGHT_PT = 470.0

def stack_axes_fractions(page_height_in=STACK_PAGE_HEIGHT_IN):
    """`(bottom, top)` for `subplots_adjust` that reproduces the shared box.

    Pages differ in height, so the fractions differ to land the box at the same points.
    """
    page_pt = page_height_in * 72.0
    top = 1.0 - STACK_AXES_TOP_OFFSET_PT / page_pt
    bottom = 1.0 - (STACK_AXES_TOP_OFFSET_PT + STACK_AXES_HEIGHT_PT) / page_pt
    if bottom <= 0.0:
        raise SystemExit(
            "a %.2f in page cannot hold the shared plot box (%.0f pt from the top, %.0f pt "
            "tall) and anything below it" % (page_height_in, STACK_AXES_TOP_OFFSET_PT,
                                             STACK_AXES_HEIGHT_PT))
    return bottom, top

def load_taxonomy(path):
    """The read-ID taxonomy master, with the union partition derived and checked."""
    sys.path.insert(0, str(HERE))
    import panel_style as ps

    frame = pd.read_csv(path, sep="\t")
    ps.require_columns(frame, TAXONOMY_REQUIRED, str(path))
    frame = frame.copy()
    frame["both_genome_unique"] = frame["n_gU_tU"] + frame["n_gU_tM"]
    frame["both_genome_multi"] = frame["n_gM_tU"] + frame["n_gM_tM"]
    frame["genome_only_unique"] = frame["n_gU_tA"]
    frame["genome_only_multi"] = frame["n_gM_tA"]
    frame["txome_only"] = frame["n_gA_tU"] + frame["n_gA_tM"]

    total = (frame["both_genome_unique"] + frame["both_genome_multi"]
             + frame["genome_only_unique"] + frame["genome_only_multi"]
             + frame["txome_only"])
    if not (total == frame["n_universe"]).all():
        bad = frame.loc[total != frame["n_universe"], "sample"].tolist()
        raise SystemExit("the union partition does not sum to n_universe for: %s"
                         % ", ".join(bad))

    frame["genome_present"] = frame["n_genome_unique"] + frame["n_genome_multi"]
    frame["txome_present"] = frame["n_txome_unique"] + frame["n_txome_multi"]
    for column, absent in (("genome_present", "n_genome_absent"),
                           ("txome_present", "n_txome_absent")):
        if not (frame[column] == frame["n_universe"] - frame[absent]).all():
            raise SystemExit("%s disagrees with n_universe - %s" % (column, absent))
    frame["delta_reads"] = frame["genome_present"] - frame["txome_present"]
    return frame

def load_labels(samples_csv=None):
    """{sample: GSM} from the sample table; keys normalised (spaces -> underscores)."""
    if not samples_csv:
        return {}
    frame = pd.read_csv(samples_csv, usecols=["cell_line", "ribo_GSM"])
    return {str(c).replace(" ", "_"): str(g)
            for c, g in zip(frame["cell_line"], frame["ribo_GSM"])}

def sample_order(frame, sort_column="delta_reads"):
    """The shared cohort ordering: cell lines by ascending `delta_reads`."""
    return frame.sort_values(sort_column)["sample"].tolist()

def median_iqr(values):
    import numpy as np
    finite = values[np.isfinite(values)]
    if finite.size == 0:
        return ""
    q1, median, q3 = np.percentile(finite, [25, 50, 75])
    return "median %.1f%%\nIQR [%.1f, %.1f]" % (median, q1, q3)
