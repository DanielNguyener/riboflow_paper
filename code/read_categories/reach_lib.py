#!/usr/bin/env python3
"""Genome-anchored reach classification."""
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
import bam_inputs as fc

OUTDIR = fc.output_root() / "read_categories"

#: How a genome-unique read that is ABSENT from the transcriptome BAM relates to the
#: selected transcript of its gene.
REACH_CATEGORIES = [
    "representable_not_present_in_dedup_bam",
    "splice_junction_absent",
    "nonselected_isoform_exon",
    "protein_coding_gene_omitted",
    "pseudogene",
    "non_protein_coding_gene",
    "intronic",
    "intergenic",
    "other_unclassified",
]

def gene_to_transcript_map(table):
    """gene_id -> transcript_id (1:1 — one selected APPRIS transcript per gene)."""
    return {v["gene_id"]: tid for tid, v in table.items()}

def omitted_pc_genes(exon_gene_df, selected_genes):
    """protein_coding genes in the GTF with NO selected (APPRIS) transcript."""
    pc = exon_gene_df[exon_gene_df["gene_type"] == "protein_coding"]
    all_pc = set(pc["gene_id"].unique())
    return all_pc - selected_genes

# ── direct overlap with omitted exonic sequence (Figure 4D) ─────────────────
# The gene-level test of Figure 5A (`read_categories/gene_read_partition_lib.
# alt_exon_overlap`) applied cohort-wide: a gU_tA read counts once per library when an
# aligned block of its primary genomic alignment overlaps exonic sequence of ANY gene
# that is absent from that gene's selected transcript. Strand-agnostic, indifferent to
# junctions and to whether the rest of the alignment fits the selected transcript.

from intervals import merge as _merge, subtract as _subtract  # noqa: E402

def build_omitted_exon_index(exon_gene_df, table):
    """{chrom: (starts, ends)}: the union, over every gene with a selected transcript, of
    the gene's exonic sequence its selected transcript omits.

    Per gene, exactly Figure 5A's `omitted_exonic_sequence`: the gene's exons are merged
    over EVERY annotated transcript (all `exon_gene_df` rows sharing the version-stripped
    gene_id, chromosome taken from the first row), and the selected transcript's exon
    intervals are subtracted. The per-gene results are unioned per chromosome, so a read
    that qualifies at several genes is still one read.
    """
    frame = exon_gene_df
    base = frame["gene_id"].astype(str).str.split(".").str[0]
    gene_index = {}
    for gene, rows in frame.groupby(base, sort=False):
        gene_index[gene] = (str(rows["Chromosome"].iat[0]),
                            _merge(zip(rows["Start"].tolist(), rows["End"].tolist())))
    by_chrom = {}
    for tid, info in table.items():
        entry = gene_index.get(str(info["gene_id"]).split(".", 1)[0])
        if entry is None:
            continue
        chrom, exons = entry
        selected = list(zip(info["g_start"].tolist(), info["g_end"].tolist()))
        by_chrom.setdefault(chrom, []).extend(_subtract(exons, selected))
    out = {}
    for chrom, intervals in by_chrom.items():
        merged = _merge(intervals)
        out[chrom] = (np.array([i[0] for i in merged], dtype=np.int64),
                      np.array([i[1] for i in merged], dtype=np.int64))
    return out

def omitted_exon_overlap_qnames(qnames, genome_blocks, omitted_by_chrom):
    """The qnames whose primary's aligned blocks overlap any omitted exonic interval."""
    import bisect

    hits = set()
    for q in qnames:
        rec = genome_blocks.get(q)
        if rec is None:
            continue
        chrom, _strand, blocks = rec
        entry = omitted_by_chrom.get(str(chrom))
        if entry is None:
            continue
        starts, ends = entry
        for b_start, b_end in blocks:
            # merged disjoint intervals sorted by start: the only candidate for a
            # `start < b_end` overlap is the last interval starting before b_end.
            j = bisect.bisect_left(starts, b_end)
            if j > 0 and ends[j - 1] > b_start:
                hits.add(q)
                break
    return hits

def _find_exon_idx(bs, be, g_start, g_end):
    """Index of the exon fully containing block [bs,be), or -1 if none."""
    hits = np.nonzero((g_start <= bs) & (be <= g_end))[0]
    return int(hits[0]) if len(hits) else -1

def representability(chrom, strand, blocks, t):
    """Classify a read (that overlaps selected gene `t`'s locus) against t's own
    selected-transcript exon structure. Returns one of:
      "representable", "nonselected_isoform_exon", "splice_junction_absent"
    """
    if chrom != t["chrom"]:
        return "nonselected_isoform_exon"
    g_start, g_end = t["g_start"], t["g_end"]
    idxs = [_find_exon_idx(bs, be, g_start, g_end) for bs, be in blocks]
    if any(i < 0 for i in idxs) or strand != t["strand"]:
        return "nonselected_isoform_exon"
    if len(idxs) == 1:
        return "representable"
    order = idxs if t["strand"] == "+" else idxs[::-1]
    if order == list(range(order[0], order[0] + len(order))):
        return "representable"
    return "splice_junction_absent"

def classify_gU_tA(qnames, genome_blocks, exon_gene_pr, exon_gene_df,
                   all_gene_body_pr, transcript_table, gene2tid, omitted_genes):
    """Return a Series qname -> category, for the gU_tA population."""
    import pyranges as pr

    rows = []
    for q in qnames:
        rec = genome_blocks.get(q)
        if rec is None:
            continue
        chrom, strand, blocks = rec
        rows.append((q, chrom, min(b[0] for b in blocks), max(b[1] for b in blocks)))
    reads_df = pd.DataFrame(rows, columns=["qname", "Chromosome", "Start", "End"])

    label = pd.Series("other_unclassified", index=reads_df["qname"], dtype=object)

    reads_pr = pr.PyRanges(reads_df.reset_index(drop=True))
    joined = reads_pr.join(exon_gene_pr, strandedness=False, how=None).df
    has_exon_hit = set(joined["qname"].unique()) if not joined.empty else set()

    no_exon = reads_df[~reads_df["qname"].isin(has_exon_hit)]
    if len(no_exon):
        no_exon_pr = pr.PyRanges(no_exon[["Chromosome", "Start", "End", "qname"]].reset_index(drop=True))
        body_joined = no_exon_pr.join(all_gene_body_pr, strandedness=False, how=None).df
        in_body = set(body_joined["qname"].unique()) if not body_joined.empty else set()
        for q in no_exon["qname"]:
            label[q] = "intronic" if q in in_body else "intergenic"

    if not joined.empty:
        def pick(sub):
            gene_ids = sub["gene_id"].tolist()
            types = sub["gene_type"].tolist()
            for gid, gt in zip(gene_ids, types):
                if gt == "protein_coding" and gid in gene2tid:
                    return ("selected_pc", gid)
            for gid, gt in zip(gene_ids, types):
                if gt == "protein_coding":
                    return ("omitted_pc", gid)
            for gt in types:
                if "pseudogene" in gt:
                    return ("pseudogene", None)
            return ("other_biotype", None)

        picks = joined.groupby("qname").apply(pick, include_groups=False)
        for q, (kind, gid) in picks.items():
            if kind == "omitted_pc":
                label[q] = "protein_coding_gene_omitted"
            elif kind == "pseudogene":
                label[q] = "pseudogene"
            elif kind == "other_biotype":
                label[q] = "non_protein_coding_gene"
            elif kind == "selected_pc":
                rec = genome_blocks.get(q)
                if rec is None:
                    continue
                chrom, strand, blocks = rec
                t = transcript_table[gene2tid[gid]]
                verdict = representability(chrom, strand, blocks, t)
                if verdict == "representable":
                    label[q] = "representable_not_present_in_dedup_bam"
                else:
                    label[q] = verdict

    return label
