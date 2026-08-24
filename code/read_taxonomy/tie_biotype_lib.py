#!/usr/bin/env python3
"""Shared genome-multimapper reads whose primary is score-tied with a secondary at
protein_coding / processed_pseudogene loci."""
from __future__ import annotations

import sys
from collections import defaultdict
from pathlib import Path

import numpy as np
import pandas as pd

_HERE = Path(__file__).resolve().parent
sys.path.insert(0, str(_HERE))
import biotype_align_lib as bal
cl = bal.cl
fc = bal.fc
tl = bal.tl
bl = bal.bl

OUTDIR = bal.OUTDIR
PC = "protein_coding"
PP = "processed_pseudogene"
_MISSING_AS = -(10 ** 9)

def read_genome_multi_records_flagged(bam_path, target_qnames):
    """qname -> [(chrom, pos5, AS, is_secondary)] for every reported genome locus of each
    NH>1-primary read in `target_qnames`. Full BAM pass; keeps the primary/secondary flag
    and per-record AS needed for the tie test."""
    import pysam
    out = defaultdict(list)
    primary_multi = set()
    bam = pysam.AlignmentFile(str(bam_path), "rb")
    for r in bam.fetch(until_eof=True):
        if r.is_unmapped or r.is_supplementary:
            continue
        q = r.query_name
        if q not in target_qnames:
            continue
        try:
            nh = r.get_tag("NH")
        except KeyError:
            nh = None
        if not r.is_secondary and nh is not None and nh > 1:
            primary_multi.add(q)
        if nh is None or nh <= 1:
            continue
        blocks = r.get_blocks()
        if not blocks:
            continue
        pos5 = blocks[0][0] if not r.is_reverse else blocks[-1][1] - 1
        AS = int(r.get_tag("AS")) if r.has_tag("AS") else _MISSING_AS
        out[q].append((r.reference_name, pos5, AS, bool(r.is_secondary)))
    bam.close()
    return {q: recs for q, recs in out.items() if q in primary_multi}

def _classify_loci(records_by_qname, exon_pr, gene_pr):
    """Flat DataFrame of every locus with its 5'-base biotype.
    Columns: qname, AS, is_secondary, biotype. One row per locus."""
    import pyranges as pr
    rows = []
    idx = 0
    for q, recs in records_by_qname.items():
        for (chrom, pos5, AS, is_sec) in recs:
            rows.append((idx, q, chrom, int(pos5), int(AS), bool(is_sec)))
            idx += 1
    if not rows:
        return pd.DataFrame(columns=["qname", "AS", "is_secondary", "biotype"])

    base = pd.DataFrame(rows, columns=["locus_idx", "qname", "Chromosome", "Start", "AS", "is_secondary"])
    loc = base[["locus_idx", "Chromosome", "Start"]].copy()
    loc["End"] = loc["Start"] + 1
    loc_pr = pr.PyRanges(loc)

    biotype = {}
    ej = loc_pr.join(exon_pr, strandedness=False, how=None).df
    if not ej.empty:
        ej = ej[["locus_idx", "gene_type"]].copy()
        ej["rank"] = ej["gene_type"].map(bl._rank_int)
        ej = ej.sort_values(["locus_idx", "rank", "gene_type"]).drop_duplicates("locus_idx")
        biotype = dict(zip(ej["locus_idx"], ej["gene_type"]))

    remaining = loc[~loc["locus_idx"].isin(biotype)]
    if len(remaining):
        rem_pr = pr.PyRanges(remaining)
        gj = rem_pr.join(gene_pr, strandedness=False, how=None).df
        intronic = set(gj["locus_idx"]) if not gj.empty else set()
        rem = remaining["locus_idx"].to_numpy()
        rem_bt = np.where(np.isin(rem, list(intronic)), "intronic", "intergenic")
        biotype.update(dict(zip(rem, rem_bt)))

    base["biotype"] = base["locus_idx"].map(biotype)
    return base[["qname", "AS", "is_secondary", "biotype"]]

def categorize_reads(records_by_qname, exon_pr, gene_pr):
    """qname -> one of the four tie categories, or None for a read qualifying for none.

    The four categories are mutually exclusive by construction, so one label per read.
    """
    loci = _classify_loci(records_by_qname, exon_pr, gene_pr)
    if loci.empty:
        return pd.Series(dtype=object)

    prim = loci[~loci["is_secondary"]].drop_duplicates("qname").set_index("qname")
    prim_bt = prim["biotype"]
    prim_as = prim["AS"]

    sec = loci[loci["is_secondary"]].copy()
    sec["prim_as"] = sec["qname"].map(prim_as)
    tied = sec[sec["AS"] == sec["prim_as"]]
    has_pc = tied[tied["biotype"] == PC].groupby("qname").size().gt(0)
    has_pp = tied[tied["biotype"] == PP].groupby("qname").size().gt(0)

    df = pd.DataFrame({"prim_bt": prim_bt})
    df["has_pc"] = df.index.map(has_pc).fillna(False).astype(bool)
    df["has_pp"] = df.index.map(has_pp).fillna(False).astype(bool)

    is_pc = df["prim_bt"] == PC
    is_pp = df["prim_bt"] == PP
    label = pd.Series(None, index=df.index, dtype=object)
    label[is_pc & df["has_pc"] & ~df["has_pp"]] = "same_pc_pc"
    label[is_pp & df["has_pp"] & ~df["has_pc"]] = "same_pp_pp"
    label[is_pc & df["has_pp"]] = "cross_pc_pp"
    label[is_pp & df["has_pc"]] = "cross_pp_pc"
    return label

def categorize(records_by_qname, exon_pr, gene_pr):
    """(counts, n_reads): the four category counts plus n_qualifying over `records_by_qname`."""
    n_reads = len(records_by_qname)
    label = categorize_reads(records_by_qname, exon_pr, gene_pr)
    counts = {name: int((label == name).sum())
              for name in ("cross_pc_pp", "cross_pp_pc", "same_pc_pc", "same_pp_pp")}
    counts["n_qualifying"] = sum(counts.values())
    return counts, n_reads
