#!/usr/bin/env python3
"""Shared loading for the two Figure-2 heatmaps: sample x read-length matrices."""
from __future__ import annotations

import sys
from pathlib import Path

import numpy as np
import pandas as pd

HERE = Path(__file__).resolve().parent

#: Grid size and margins in inches, SHARED by both Figure-2 panels so cells and pages match
#: 1:1 at assembly. Each panel leaves the other's gutter BLANK, so both must save with
#: `tight=False` or the crop takes the blank margins back off.
AXES_SIZE = (7.0, 7.8)
MARGINS = (1.25, 1.75, 1.25, 0.15)          # left, bottom, right, top

def _as_bool(series):
    return series.astype(str).str.lower().isin(("true", "1"))

def load_selected(qc_path, value_column="psite_offset"):
    """Long table of phase-1 read lengths and one per-length value, per sample."""
    sys.path.insert(0, str(HERE))
    import panel_style as ps

    frame = pd.read_csv(qc_path)
    ps.require_columns(frame, ("sample", "read_length", "in_phase1", value_column),
                       str(qc_path))
    frame = frame[_as_bool(frame["in_phase1"])].copy()
    frame["read_length"] = frame["read_length"].astype(int)
    return frame[["sample", "read_length", value_column]]

def to_matrix(frame, samples, lengths, value_column):
    """sample x read_length matrix; NaN where that length was not selected."""
    matrix = np.full((len(samples), len(lengths)), np.nan)
    sample_index = {s: i for i, s in enumerate(samples)}
    length_index = {L: j for j, L in enumerate(lengths)}
    for sample, length, value in zip(frame["sample"], frame["read_length"],
                                     frame[value_column]):
        matrix[sample_index[sample], length_index[length]] = value
    return matrix

def load_pair(genome_qc, txome_qc, value_column):
    """Both routes as aligned matrices, over the union of samples and read lengths."""
    genome = load_selected(genome_qc, value_column)
    txome = load_selected(txome_qc, value_column)
    samples = sorted(set(genome["sample"]) | set(txome["sample"]))
    lengths = sorted(set(genome["read_length"]) | set(txome["read_length"]))
    return (samples, lengths,
            to_matrix(genome, samples, lengths, value_column),
            to_matrix(txome, samples, lengths, value_column))

def gsm_labels(samples_csv, samples):
    """GSM accessions for the y axis, tolerating spaced or underscored cell-line names."""
    frame = pd.read_csv(samples_csv, usecols=["cell_line", "ribo_GSM"])
    mapping = {}
    for cell_line, gsm in zip(frame["cell_line"], frame["ribo_GSM"]):
        mapping[cell_line] = gsm
        mapping[str(cell_line).replace(" ", "_")] = gsm
    return [mapping.get(s, s) for s in samples]

def draw_grid(axis, matrix, cmap, norm):
    """Draw a sample x read-length grid as TRUE VECTOR QUADS; row 0 at the TOP.

    NOT `imshow` -- its embedded raster gets resampled by PDF/SVG viewers into mottled cells.
    """
    n_rows, n_cols = matrix.shape
    mesh = axis.pcolormesh(np.arange(n_cols + 1), np.arange(n_rows + 1),
                           np.ma.masked_invalid(matrix), cmap=cmap, norm=norm,
                           shading="flat", edgecolors="white", linewidth=0.8,
                           antialiased=False)
    axis.set_xlim(0, n_cols)
    axis.set_ylim(n_rows, 0)
    axis.set_aspect("auto")
    return mesh

def cell_centre(index):
    """pcolormesh cells span [i, i+1), so the centre -- for ticks and text -- is i + 0.5."""
    return np.asarray(index) + 0.5

#: Type for the two grids at their two ship scales: "large" = standalone panel,
#: "base" = journal page (PLOS 8-12 pt window).
def grid_type(scale="large"):
    sys.path.insert(0, str(HERE))
    import panel_style as ps
    if scale == "large":
        return {"label": ps.FONT_LABEL_LARGE, "tick": ps.FONT_TICK_LARGE,
                "annotation": ps.FONT_ANNOTATION_LARGE}
    if scale == "base":
        return {"label": ps.FONT_LABEL, "tick": ps.FONT_TICK, "annotation": ps.FONT_INSET}
    raise SystemExit("unknown type scale %r (large|base)" % scale)

def style_grid(axis, samples, lengths, labels, sizes=None, show_ylabels=True):
    """Ticks and labels at cell centres; sizes default to panel_style's *_LARGE scale."""
    sizes = sizes or grid_type("large")

    axis.set_xticks(cell_centre(range(len(lengths))))
    axis.set_xticklabels(lengths, fontsize=sizes["tick"])
    axis.set_yticks(cell_centre(range(len(samples))))
    axis.set_yticklabels(labels if show_ylabels else [""] * len(samples),
                         fontsize=sizes["tick"])
    if not show_ylabels:
        axis.tick_params(axis="y", length=0)
    axis.set_xlabel("read length (nt)", fontsize=sizes["label"])
    axis.tick_params(which="minor", length=0)
