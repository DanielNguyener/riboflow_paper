#!/usr/bin/env python3
"""Figure 4 C and D -- one cohort share per cell line, drawn as a side panel.

C: protein-coding to pseudogene tie share of genome multimapper reads (`--tie-master`).
D: omitted alternative-exon overlap of genome-only unique reads (`--reach-master`):
   the Figure 5A gene-level test, cohort-wide -- a read counts once when an aligned block
   of its primary genomic alignment overlaps exonic sequence absent from the overlapped
   gene's selected transcript.
"""
from __future__ import annotations

import sys
from pathlib import Path

import numpy as np
import pandas as pd

HERE = Path(__file__).resolve().parent

#: Per panel: the master's flag, its required columns, the share, and the panel's dress.
PANELS = {
    "C": {"master": "tie_master",
          "required": ("sample", "pct_cross_pp_pc", "pct_cross_pc_pp"),
          "share": lambda f: f["pct_cross_pp_pc"] + f["pct_cross_pc_pp"],
          "stat": "cross tie", "colour": "#1a7d1a",
          "title": "protein-coding ↔\npseudogene tie",
          "xlabel": "% of multimapper\nreads (primary)"},
    "D": {"master": "reach_master",
          "required": ("sample", "n_omitted_exon_overlap", "n_gU_tA"),
          "share": lambda f: 100.0 * f["n_omitted_exon_overlap"] / f["n_gU_tA"],
          "stat": "omitted-exon overlap", "colour": "#7fb9da",
          "title": "omitted\nalternative exon",
          "xlabel": "% of genome-only\nunique reads"},
}

def prepare(letter, master, taxonomy, samples_csv=None):
    sys.path.insert(0, str(HERE))
    import cohort_common as common
    import panel_style as ps

    spec = PANELS[letter]
    frame = pd.read_csv(master, sep="\t")
    ps.require_columns(frame, spec["required"], str(master))
    frame = frame.copy()
    frame["share"] = spec["share"](frame)

    order = common.sample_order(common.load_taxonomy(taxonomy))
    lookup = frame.set_index("sample")["share"]
    values = np.array([lookup.get(s, np.nan) for s in order], dtype=float)
    labels_map = common.load_labels(samples_csv)
    return {"order": order, "values": values,
            "order_provenance": "shared cohort order, derived from %s" % taxonomy,
            "labels": [labels_map.get(s, s) for s in order]}
