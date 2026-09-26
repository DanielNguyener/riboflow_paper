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

GENOME_STATES = ("unique", "multi", "absent")
TXOME_STATES = ("present", "absent")
#: The five cells of the read-ID taxonomy: genome status x transcriptome presence, minus
#: the empty (absent, absent) cell. Transcriptome presence is a primary alignment in the
#: post-dedup BAM; RiboFlow_v2 already filtered that BAM at MAPQ >= 10, so no further
#: threshold is applied here (the MAPQ >= 42 rule belongs to the coverage/TE analyses).
CELLS = tuple((g, t) for g in GENOME_STATES for t in TXOME_STATES
              if not (g == "absent" and t == "absent"))
ABBR = {"unique": "U", "multi": "M", "present": "P", "absent": "A"}

def cell_key(genome_status, txome_status):
    """`gU_tP`, `gM_tA`, ...: the column stem of one taxonomy cell."""
    return "g%s_t%s" % (ABBR[genome_status], ABBR[txome_status])

def genome_status_sets(bam_path):
    """(all_qnames, unique_qnames) over the primary alignments of a genome BAM;
    unique = NH == 1 (`bam_inputs.is_unique_genome_read`)."""
    all_q, uniq_q = set(), set()
    bam = pysam.AlignmentFile(str(bam_path), "rb")
    for r in bam.fetch(until_eof=True):
        if r.is_unmapped or r.is_secondary or r.is_supplementary:
            continue
        all_q.add(r.query_name)
        if fc.is_unique_genome_read(r):
            uniq_q.add(r.query_name)
    bam.close()
    return all_q, uniq_q

def txome_present_qnames(bam_path):
    """Every read id with a primary alignment in the transcriptome BAM."""
    all_q = set()
    bam = pysam.AlignmentFile(str(bam_path), "rb")
    for r in bam.fetch(until_eof=True):
        if r.is_unmapped or r.is_secondary or r.is_supplementary:
            continue
        all_q.add(r.query_name)
    bam.close()
    return all_q

def classify_sample(sample, log=print):
    """Return (counts keyed (genome_status, txome_status), n_universe) for one sample."""
    g_all, g_uniq = genome_status_sets(fc.genome_bam(sample))
    t_all = txome_present_qnames(fc.txome_bam(sample))
    if log:
        log(f"  [{sample}] genome mapped={len(g_all):,} unique={len(g_uniq):,} | "
            f"txome mapped={len(t_all):,}")

    inter = len(g_all & t_all)
    smaller = min(len(g_all), len(t_all))
    frac = inter / smaller if smaller else 0.0
    assert frac > 0.05, (
        f"[{sample}] QNAME-namespace mismatch: genome n txome intersection {inter:,} "
        f"= {frac:.1%} of min({len(g_all):,},{len(t_all):,}) — suffix/dedup bug suspected")
    if frac < 0.60 and log:
        log(f"  [{sample}] WARNING: low qname overlap ({frac:.1%} of smaller route) — "
            f"asymmetric read recovery; expect a large `absent` fraction")

    counts = {cell: 0 for cell in CELLS}
    for q in g_all:
        gs = "unique" if q in g_uniq else "multi"
        counts[(gs, "present" if q in t_all else "absent")] += 1
    counts[("absent", "present")] = len(t_all - g_all)

    n_universe = len(g_all | t_all)
    assert sum(counts.values()) == n_universe, \
        f"[{sample}] cell sum {sum(counts.values()):,} != n_universe {n_universe:,}"
    return counts, n_universe
