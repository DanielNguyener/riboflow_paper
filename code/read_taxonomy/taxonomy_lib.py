#!/usr/bin/env python3
"""The multi-mapping read taxonomy."""
from __future__ import annotations

import os
import sys
from pathlib import Path

import pandas as pd
import pysam

_HERE = Path(__file__).resolve().parent
_COMMON = _HERE.parent / "common"
for _entry in (str(_HERE), str(_COMMON), str(_COMMON / "ribo_seq_qc")):
    if _entry not in sys.path:
        sys.path.insert(0, _entry)
import bam_inputs as fc

def sample_to_gsm(samples_csv=None):
    """sample name -> ribo_GSM, from the sample table.

    `cell_line` uses spaces where sample names use underscores, so the key is normalised.
    """
    if samples_csv is None:
        samples_csv = os.environ.get("RIBOFLOW_PAPER_SAMPLES_CSV")
    if samples_csv is None:
        samples_csv = (Path(__file__).resolve().parents[2] / "supporting_information"
                       / "S1_Table" / "samples.csv")
    samples_csv = Path(samples_csv)
    if not samples_csv.exists():
        raise SystemExit(
            "sample table not found: %s\nPass an explicit path, set "
            "RIBOFLOW_PAPER_SAMPLES_CSV, or keep the table at "
            "supporting_information/S1_Table/samples.csv" % samples_csv)
    frame = pd.read_csv(samples_csv)
    return {row["cell_line"].replace(" ", "_"): row["ribo_GSM"]
            for _, row in frame.iterrows()}

STATES = ("unique", "multi", "absent")

def status_sets(bam_path, kind):
    """(all_mapped_qnames, unique_qnames) over PRIMARY alignments of one BAM.

    kind="genome": unique = NH==1 (bam_inputs.is_unique_genome_read).
    kind="txome" : unique = MAPQ >= TXOME_MIN_MAPQ (bowtie2, no NH tag).
    """
    if kind not in ("genome", "txome"):
        raise ValueError(f"kind must be 'genome' or 'txome', got {kind!r}")
    all_q, uniq_q = set(), set()
    bam = pysam.AlignmentFile(str(bam_path), "rb")
    for r in bam.fetch(until_eof=True):
        if r.is_unmapped or r.is_secondary or r.is_supplementary:
            continue
        q = r.query_name
        all_q.add(q)
        if kind == "genome":
            uniq = fc.is_unique_genome_read(r)
        else:
            uniq = r.mapping_quality >= fc.txome_min_mapq()
        if uniq:
            uniq_q.add(q)
    bam.close()
    return all_q, uniq_q

def status_of(q, all_q, uniq_q):
    if q in uniq_q:
        return "unique"
    if q in all_q:
        return "multi"
    return "absent"

def classify_sample(sample, log=print):
    """Return (counts keyed (genome_status, txome_status), n_universe) for one sample.

    Memory-lean: only the four qname sets are held, never a fifth `universe` set.
    """
    g_all, g_uniq = status_sets(fc.genome_bam(sample), "genome")
    t_all, t_uniq = status_sets(fc.txome_bam(sample), "txome")
    if log:
        log(f"  [{sample}] genome mapped={len(g_all):,} unique={len(g_uniq):,} | "
            f"txome mapped={len(t_all):,} unique={len(t_uniq):,}")

    inter = len(g_all & t_all)
    smaller = min(len(g_all), len(t_all))
    frac = inter / smaller if smaller else 0.0
    assert frac > 0.05, (
        f"[{sample}] QNAME-namespace mismatch: genome n txome intersection {inter:,} "
        f"= {frac:.1%} of min({len(g_all):,},{len(t_all):,}) — suffix/dedup bug suspected")
    if frac < 0.60 and log:
        log(f"  [{sample}] WARNING: low qname overlap ({frac:.1%} of smaller route) — "
            f"asymmetric read recovery; expect a large `absent` fraction")

    counts = {(gs, ts): 0 for gs in STATES for ts in STATES if not (gs == "absent" and ts == "absent")}
    for q in g_all:
        gs = status_of(q, g_all, g_uniq)
        ts = status_of(q, t_all, t_uniq)
        counts[(gs, ts)] += 1
    for q in t_all:
        if q in g_all:
            continue
        ts = status_of(q, t_all, t_uniq)
        counts[("absent", ts)] += 1

    n_universe = len(g_all | t_all)
    assert sum(counts.values()) == n_universe, \
        f"[{sample}] cell sum {sum(counts.values()):,} != n_universe {n_universe:,}"
    return counts, n_universe
