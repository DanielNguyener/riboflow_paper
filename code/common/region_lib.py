#!/usr/bin/env python3
"""Transcript-region classification: the five-way ribopy scheme, applied to both routes."""
from __future__ import annotations

import os
import re
import sys
from pathlib import Path

import numpy as np
import pandas as pd

# bam_inputs sits beside this file and gives `config`, the BAM accessors, the MAPQ
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

    OPTIONAL, and takes an explicit path: nothing in this repository calls it as part of
    an analysis. It exists so that someone holding the original `.ribo` can confirm the
    spans this code uses match the ones baked into that file. Use `DEFAULT_LEFT_SPAN` /
    `DEFAULT_RIGHT_SPAN` otherwise.
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
    """Vectorised region code for transcript-coord 5′-ends `x`.

    All args except the two spans are 1-D numpy arrays of equal length; returns an
    int array with 0=UTR5 1=UTR5J 2=CDS 3=UTR3J 4=UTR3 (indexes into REGIONS).
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

def txome_region_map(bam) -> dict:
    """refname -> (base_tid, start_site, stop_site, L) for every reference with a CDS.

    start_site = CDS:start - 1 (0-based), stop_site = CDS:end (0-based exclusive),
    L = reference length. Returns only refs whose header carries UTR5+CDS.
    """
    out = {}
    lengths = dict(zip(bam.references, bam.lengths))
    for ref in bam.references:
        mc = _CDS_RE.search(ref)
        m5 = _UTR5_RE.search(ref)
        if not mc or not m5:
            continue
        start_site = int(mc.group(1)) - 1
        stop_site = int(mc.group(2))
        base = ref.split(".", 1)[0].split("|", 1)[0]
        out[ref] = (base, start_site, stop_site, lengths[ref])
    return out

def build_genome_exon_table(allowed_base_ids: set | None = None):
    """Per-(transcript, region) exon table with 5′→3′ cumulative offsets.

    region ∈ {five, cds, three}; cum_offset = nt from that region's 5′ boundary.
    Restricted to APPRIS principal isoforms that have all three regions
    (`length_filtered == False`), optionally further restricted to
    `allowed_base_ids` (base ENST, no version).

    THIS IS NOT A UNIVERSE. The three-region requirement reproduces RiboPy's own
    `get_length_dist("CDS")` population, which is what the read-length selection must be
    computed over; `qc_core.genome_cds_core_intervals` is its only caller. Figure 4 counts
    into the CANONICAL CDS over every APPRIS transcript that has one -- a different region
    and a different transcript set, built in `ribo_rna_lib` from the annotation bundle.

    Returns (exons_df, base2ver) where base2ver maps base ENST → versioned id.
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
