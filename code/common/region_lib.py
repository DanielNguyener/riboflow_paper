#!/usr/bin/env python3
"""Transcript-region classification: the five-way ribopy scheme, applied to both routes."""
from __future__ import annotations

import re
import sys
from pathlib import Path

import numpy as np
import pandas as pd

# bam_inputs sits beside this file and provides `config` plus the BAM accessors
_HERE = Path(__file__).resolve().parent
for _entry in (str(_HERE), str(_HERE / "ribo_seq_qc")):
    if _entry not in sys.path:
        sys.path.insert(0, _entry)
import bam_inputs as fc

_UTR5_RE = re.compile(r"\|UTR5:(\d+)-(\d+)\|")
_CDS_RE = re.compile(r"\|CDS:(\d+)-(\d+)\|")

REGIONS = ["UTR5", "UTR5J", "CDS", "UTR3J", "UTR3"]

DEFAULT_LEFT_SPAN = 35
DEFAULT_RIGHT_SPAN = 10

def load_ribo_params(ribo_path) -> dict:
    """Read left_span / right_span / length_min / length_max out of a `.ribo` file.

    Intentionally uncalled: kept so a holder of the original `.ribo` can confirm the spans.
    """
    import h5py

    if ribo_path is None:
        raise ValueError(
            "load_ribo_params needs an explicit .ribo path. No .ribo file is shipped "
            "with this repository; the spans default to DEFAULT_LEFT_SPAN=%d / "
            "DEFAULT_RIGHT_SPAN=%d." % (DEFAULT_LEFT_SPAN, DEFAULT_RIGHT_SPAN))
    with h5py.File(ribo_path, "r") as h:
        p = {k: int(h.attrs[k]) for k in
             ("left_span", "right_span", "length_min", "length_max")}
    p["ribo_path"] = ribo_path
    return p

# ── region classifier (vectorised port of get_extended_boundaries) ──────────────
def classify(x, start_site, stop_site, left_span, right_span):
    """Vectorised region code (0=UTR5 1=UTR5J 2=CDS 3=UTR3J 4=UTR3) for 5′-ends `x`.

    Mirrors ribopy's half-open boundaries exactly (second coordinate excluded).
    """
    x = np.asarray(x)
    s = np.asarray(start_site)
    e = np.asarray(stop_site)
    u5j_lo = s - left_span
    cds_lo = s + right_span + 1
    cds_hi = e - left_span
    u3j_hi = e + right_span + 1

    code = np.full(x.shape, -1, dtype=np.int8)
    code[(x >= 0) & (x < u5j_lo)] = 0
    code[(x >= u5j_lo) & (x < cds_lo)] = 1
    code[(x >= cds_lo) & (x < cds_hi)] = 2
    code[(x >= cds_hi) & (x < u3j_hi)] = 3
    code[(x >= u3j_hi)] = 4
    return code

def tally(codes) -> np.ndarray:
    """5-element count vector [UTR5, UTR5J, CDS, UTR3J, UTR3] from classify() codes."""
    return np.array([int((codes == i).sum()) for i in range(5)], dtype=np.int64)

def build_genome_exon_table(allowed_base_ids: set | None = None):
    """Per-(transcript, region) exon table with 5′→3′ cumulative offsets → (exons_df, base2ver).

    NOT a universe: the three-region requirement reproduces RiboPy's `get_length_dist("CDS")`
    population for read-length selection only; Figure 4 uses a different transcript set.
    """
    _as_bool = fc._as_bool
    cds_df = fc.config.load_annotation()
    utr_df = fc.config.load_appris_utr()
    meta_df = fc.config.load_appris_meta()

    length_ok = set(meta_df.loc[~_as_bool(meta_df["length_filtered"]),
                                "transcript_id"])

    cds_part = cds_df[["Chromosome", "Start", "End", "Strand", "transcript_id"]].copy()
    cds_part["region"] = "cds"
    utr_part = utr_df[["Chromosome", "Start", "End", "Strand", "transcript_id"]].copy()
    utr_part["region"] = utr_df["utr_type"].values

    exons = pd.concat([cds_part, utr_part], ignore_index=True)
    exons = exons[exons["transcript_id"].isin(length_ok)].copy()

    base2ver = {tid.split(".", 1)[0]: tid for tid in length_ok}
    if allowed_base_ids is not None:
        keep_ver = {base2ver[b] for b in allowed_base_ids if b in base2ver}
        exons = exons[exons["transcript_id"].isin(keep_ver)].copy()

    exons["exon_len"] = exons["End"] - exons["Start"]
    exons["order_key"] = np.where(exons["Strand"] == "+",
                                  exons["Start"], -exons["Start"])
    exons = exons.sort_values(["transcript_id", "region", "order_key"]) \
                 .reset_index(drop=True)
    grp = exons.groupby(["transcript_id", "region"], sort=False)["exon_len"]
    exons["cum_offset"] = grp.cumsum() - exons["exon_len"]
    return exons, base2ver
