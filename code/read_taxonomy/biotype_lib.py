#!/usr/bin/env python3
"""GENCODE gene_type of genome-only unique reads."""
from __future__ import annotations

import sys
from pathlib import Path

_HERE = Path(__file__).resolve().parent
_COMMON = _HERE.parent / "common"
for _entry in (str(_HERE), str(_COMMON), str(_COMMON / "ribo_seq_qc")):
    if _entry not in sys.path:
        sys.path.insert(0, _entry)
import reference_lib as cl
fc = cl.fc

OUTDIR = fc.output_root() / "read_taxonomy" / "genome_only_biotype"
CACHE_DIR = fc.output_root() / ".cache" / "read_taxonomy"
GENE_BODY_CACHE = CACHE_DIR / "gene_body.pkl"

def _rank_int(gt):
    if gt == "protein_coding":
        return 0
    if gt == "lncRNA":
        return 1
    if "pseudogene" in gt:
        return 2
    return 3

def gene_body_pr(rebuild=False):
    """PyRanges of per-gene genomic bodies (min exon start .. max exon end) with gene_type.

    Cached keyed by the annotation fingerprint, so a changed GTF rebuilds it.
    """
    import pyranges as pr

    def build():
        exons = cl.build_exon_gene_table()
        return exons.groupby("gene_id", sort=False).agg(
            Chromosome=("Chromosome", "first"),
            Start=("Start", "min"),
            End=("End", "max"),
            gene_type=("gene_type", "first"),
        ).reset_index()

    if rebuild and GENE_BODY_CACHE.exists():
        GENE_BODY_CACHE.unlink()
    frame = fc.config.cached_frame(
        GENE_BODY_CACHE,
        fc.config.annotation_fingerprint(["gene_body_pr/1"]),
        build)
    return pr.PyRanges(frame)
