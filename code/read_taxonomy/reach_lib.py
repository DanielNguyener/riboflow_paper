#!/usr/bin/env python3
"""Genome-anchored reach classification."""
from __future__ import annotations

import sys
from pathlib import Path

import numpy as np
import pandas as pd
import pysam

_HERE = Path(__file__).resolve().parent
_COMMON = _HERE.parent / "common"
for _entry in (str(_HERE), str(_COMMON), str(_COMMON / "ribo_seq_qc")):
    if _entry not in sys.path:
        sys.path.insert(0, _entry)
import bam_inputs as fc
import concordance_lib as cl

OUTDIR = fc.output_root() / "read_taxonomy" / "reach"

REACH_CATEGORIES = [
    "shared_unique_concordant",
    "shared_unique_discordant",
    "genome_unique_transcriptome_multimapped",
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

CATEGORY_LABEL = {
    "shared_unique_concordant": "concordant (both unique, same locus)",
    "shared_unique_discordant": "discordant (both unique, different locus)",
    "genome_unique_transcriptome_multimapped": "genome unique, transcriptome multimapped",
    "representable_not_present_in_dedup_bam": "representable, not in dedup'd txome BAM",
    "splice_junction_absent": "splice junction absent from selected isoform",
    "nonselected_isoform_exon": "exon of a nonselected isoform",
    "protein_coding_gene_omitted": "protein-coding gene omitted from transcriptome",
    "pseudogene": "pseudogene",
    "non_protein_coding_gene": "non-protein-coding gene",
    "intronic": "intronic",
    "intergenic": "intergenic",
    "other_unclassified": "other/unclassified",
}

UNREACHABLE_CATEGORIES = [
    "splice_junction_absent", "nonselected_isoform_exon", "protein_coding_gene_omitted",
    "pseudogene", "non_protein_coding_gene", "intronic", "intergenic", "other_unclassified",
]

def genome_status_sets(bam_path):
    """(all_q, uniq_q) over primary genome alignments — identical rule to
    `taxonomy_lib.status_sets(kind="genome")`."""
    all_q, uniq_q = set(), set()
    bam = pysam.AlignmentFile(str(bam_path), "rb")
    for r in bam.fetch(until_eof=True):
        if r.is_unmapped or r.is_secondary or r.is_supplementary:
            continue
        q = r.query_name
        all_q.add(q)
        uniq = fc.is_unique_genome_read(r)
        if uniq:
            uniq_q.add(q)
    bam.close()
    return all_q, uniq_q

def txome_all_qnames(bam_path):
    """All primary txome qnames (any MAPQ) — used to define gU_tA = genome-unique
    minus this set (matches `taxonomy_lib`'s txome `all_q`)."""
    all_q = set()
    bam = pysam.AlignmentFile(str(bam_path), "rb")
    for r in bam.fetch(until_eof=True):
        if r.is_unmapped or r.is_secondary or r.is_supplementary:
            continue
        all_q.add(r.query_name)
    bam.close()
    return all_q

def read_genome_blocks(bam_path, qnames):
    """qname -> (chrom, strand, blocks) for the given qname set only (gU_tA is small,
    ~1M reads/sample — cheap to keep full per-block coordinates, unlike the
    all-genome-unique pass in concordance_lib.read_genome_unique)."""
    out = {}
    bam = pysam.AlignmentFile(str(bam_path), "rb")
    for r in bam.fetch(until_eof=True):
        if r.is_unmapped or r.is_secondary or r.is_supplementary:
            continue
        if r.query_name not in qnames:
            continue
        blocks = r.get_blocks()
        if not blocks:
            continue
        strand = "-" if r.is_reverse else "+"
        out[r.query_name] = (r.reference_name, strand, blocks)
    bam.close()
    return out

def gene_to_transcript_map(table):
    """gene_id -> transcript_id (1:1 — one selected APPRIS transcript per gene)."""
    return {v["gene_id"]: tid for tid, v in table.items()}

def omitted_pc_genes(exon_gene_df, selected_genes):
    """protein_coding genes in the GTF with NO selected (APPRIS) transcript."""
    pc = exon_gene_df[exon_gene_df["gene_type"] == "protein_coding"]
    all_pc = set(pc["gene_id"].unique())
    return all_pc - selected_genes

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
