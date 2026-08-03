#!/usr/bin/env python3
"""GENCODE gene_type of genome-MULTIMAPPER ALIGNMENTS."""
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
import biotype_lib as bl
import mm_concordance_lib as mm
cl, fc, tl = bl.cl, bl.fc, bl.tl

OUTDIR = fc.output_root() / "read_taxonomy" / "multimap_biotype"

def population_records(sample, log=print):
    """qname -> [genome-locus records] for the dark-green population (gM_tU + gM_tM).

    Reads the txome dedup BAM for the present-qname set, then enumerates every reported
    genome locus of each read whose genome PRIMARY is a multimapper (NH>1) and present in
    that set. Returned dict's keys are exactly the dark-green reads; each value is the list
    of (chrom, strand, pos5, n_blocks, blk_min, blk_max, AS) tuples from `mm_concordance_lib`.
    """
    log(f"[{sample}] reading txome BAM (present qname set)...")
    t_all, _ = tl.status_sets(fc.txome_bam(sample), "txome")
    log(f"[{sample}] n_txome_present={len(t_all):,}; enumerating genome multimapper loci...")
    records_by_qname = mm.read_genome_multi_records(fc.genome_bam(sample), t_all)
    return records_by_qname

def read_genome_multi_primary(bam_path, target_qnames):
    """qname -> (chrom, pos5) for the PRIMARY genome record of each read whose primary is a
    multimapper (NH>1) and whose qname is present in `target_qnames` (the txome-present set)
    — i.e. exactly the gM_tU + gM_tM reads, ONE 5' base each (the primary locus only).

    Per-read, primary-only counterpart of population_records (which keeps every locus).
    Mirrors `biotype_lib.read_genome_unique_absent` but selects genome-multi primaries and
    INCLUDES only txome-present qnames. Genome-multi is decided by the primary's NH>1 (same
    rule as `read_genome_multi_records` / `taxonomy_lib`), and the 5'-most-base
    convention matches Parts 14/15/16.
    """
    import pysam
    out = {}
    bam = pysam.AlignmentFile(str(bam_path), "rb")
    for r in bam.fetch(until_eof=True):
        if r.is_unmapped or r.is_secondary or r.is_supplementary:
            continue
        q = r.query_name
        if q not in target_qnames:
            continue
        try:
            nh = r.get_tag("NH")
        except KeyError:
            nh = None
        if nh is None or nh <= 1:
            continue
        blocks = r.get_blocks()
        if not blocks:
            continue
        pos5 = blocks[0][0] if not r.is_reverse else blocks[-1][1] - 1
        out[q] = (r.reference_name, pos5)
    bam.close()
    return out

def classify_alignments(records_by_qname, exon_pr, gene_pr):
    """Classify EVERY genome locus (alignment) by GENCODE gene_type / intronic / intergenic.

    Generalizes `biotype_lib.classify_5p` to key on each alignment instead of each read: one row
    per locus (aln_idx), classified by the 5' base. Returns (counts, n_reads, n_alignments)
    where counts is {biotype -> n_alignments} and sum(counts.values()) == n_alignments.
    """
    import pyranges as pr

    rows = []
    aln_idx = 0
    for recs in records_by_qname.values():
        for (chrom, _strand, pos5, *_rest) in recs:
            rows.append((aln_idx, chrom, int(pos5)))
            aln_idx += 1
    n_reads = len(records_by_qname)
    n_alignments = len(rows)
    if n_alignments == 0:
        return {}, n_reads, 0

    reads = pd.DataFrame(rows, columns=["aln_idx", "Chromosome", "Start"])
    reads["End"] = reads["Start"] + 1
    reads_pr = pr.PyRanges(reads)

    biotype = {}
    ej = reads_pr.join(exon_pr, strandedness=False, how=None).df
    if not ej.empty:
        ej = ej[["aln_idx", "gene_type"]].copy()
        ej["rank"] = ej["gene_type"].map(bl._rank_int)
        ej = ej.sort_values(["aln_idx", "rank", "gene_type"]).drop_duplicates("aln_idx")
        biotype = dict(zip(ej["aln_idx"], ej["gene_type"]))

    remaining = reads[~reads["aln_idx"].isin(biotype)]
    if len(remaining):
        rem_pr = pr.PyRanges(remaining)
        gj = rem_pr.join(gene_pr, strandedness=False, how=None).df
        intronic = set(gj["aln_idx"]) if not gj.empty else set()
        rem = remaining["aln_idx"].to_numpy()
        rem_bt = np.where(np.isin(rem, list(intronic)), "intronic", "intergenic")
        biotype.update(dict(zip(rem, rem_bt)))

    counts = pd.Series(list(biotype.values())).value_counts().to_dict()
    assert sum(counts.values()) == n_alignments, (sum(counts.values()), n_alignments)
    return counts, n_reads, n_alignments
