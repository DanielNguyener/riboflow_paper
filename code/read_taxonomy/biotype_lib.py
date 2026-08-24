#!/usr/bin/env python3
"""GENCODE gene_type of genome-only unique reads."""
from __future__ import annotations

import sys
from pathlib import Path

import numpy as np
import pandas as pd

_HERE = Path(__file__).resolve().parent
_COMMON = _HERE.parent / "common"
for _entry in (str(_HERE), str(_COMMON), str(_COMMON / "ribo_seq_qc")):
    if _entry not in sys.path:
        sys.path.insert(0, _entry)
import concordance_lib as cl
fc = cl.fc

OUTDIR = fc.output_root() / "read_taxonomy" / "genome_only_biotype"
CACHE_DIR = fc.output_root() / ".cache" / "read_taxonomy"
GENE_BODY_CACHE = CACHE_DIR / "gene_body.pkl"

def read_genome_unique_absent(bam_path, exclude_qnames):
    """qname -> (chrom, pos5) for primary unique (NH==1) genome reads NOT in `exclude_qnames`.

    Same uniqueness rule and 5'-most-base convention as concordance_lib.read_genome_unique.
    """
    import pysam
    out = {}
    bam = pysam.AlignmentFile(str(bam_path), "rb")
    for r in bam.fetch(until_eof=True):
        if r.is_unmapped or r.is_secondary or r.is_supplementary:
            continue
        q = r.query_name
        if q in exclude_qnames:
            continue
        uniq = fc.is_unique_genome_read(r)
        if not uniq:
            continue
        blocks = r.get_blocks()
        if not blocks:
            continue
        pos5 = blocks[0][0] if not r.is_reverse else blocks[-1][1] - 1
        out[q] = (r.reference_name, pos5)
    bam.close()
    return out

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

def classify_5p(qnames, genome_dict, exon_pr, gene_pr):
    """qname -> biotype (gene_type | 'intronic' | 'intergenic') for each read's 5'-end base."""
    import pyranges as pr
    if not qnames:
        return pd.Series(dtype=object, name="biotype")

    rows = [(q, genome_dict[q][0], int(genome_dict[q][1])) for q in qnames]
    reads = pd.DataFrame(rows, columns=["qname", "Chromosome", "Start"])
    reads["End"] = reads["Start"] + 1
    reads_pr = pr.PyRanges(reads)

    biotype = {}
    ej = reads_pr.join(exon_pr, strandedness=False, how=None).df
    if not ej.empty:
        ej = ej[["qname", "gene_type"]].copy()
        ej["rank"] = ej["gene_type"].map(_rank_int)
        ej = ej.sort_values(["qname", "rank", "gene_type"]).drop_duplicates("qname")
        biotype = dict(zip(ej["qname"], ej["gene_type"]))

    remaining = reads[~reads["qname"].isin(biotype)]
    if len(remaining):
        rem_pr = pr.PyRanges(remaining)
        gj = rem_pr.join(gene_pr, strandedness=False, how=None).df
        intronic = set(gj["qname"]) if not gj.empty else set()
        rem_qn = remaining["qname"]
        rem_bt = np.where(rem_qn.isin(intronic).to_numpy(), "intronic", "intergenic")
        biotype.update(dict(zip(rem_qn.to_numpy(), rem_bt)))

    return pd.Series(biotype, name="biotype")
