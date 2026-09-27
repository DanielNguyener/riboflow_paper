#!/usr/bin/env python3
"""The transcriptome reference in genome coordinates: the APPRIS transcript table
(exon layout, cumulative positions) and the GENCODE exon-gene table, cached under
results/.cache/, plus the one-pass reader of a transcriptome BAM's primaries."""
from __future__ import annotations

import gzip
import re
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

CACHE_DIR = fc.output_root() / ".cache" / "read_categories"
EXON_GENE_CACHE = CACHE_DIR / "exon_gene_table.pkl"
TRANSCRIPT_TABLE_CACHE = CACHE_DIR / "transcript_coord_table.pkl"

def read_txome_primary(bam_path, base2ver):
    """One pass over every primary alignment of a transcriptome BAM (RiboFlow_v2 already
    filtered it at MAPQ >= 10; nothing more is applied) -> (present, all_qnames):
    {qname: transcript_id} for the APPRIS-resolving primaries, and every primary qname."""
    present = {}
    all_qnames = set()
    bam = pysam.AlignmentFile(str(bam_path), "rb")
    for r in bam.fetch(until_eof=True):
        if r.is_unmapped or r.is_secondary or r.is_supplementary:
            continue
        all_qnames.add(r.query_name)
        base = r.reference_name.split(".", 1)[0].split("|", 1)[0]
        tid = base2ver.get(base)
        if tid is not None:
            present[r.query_name] = tid
    bam.close()
    return present, all_qnames

def _merge_contiguous(exons):
    """Merge genomically contiguous CDS/UTR rows (one exon split at the start/stop codon).

    Without this, the CDS/UTR seam counts as a splice junction and inflates splice_discordant.
    """
    out = []
    for tid, g in exons.groupby("transcript_id", sort=False):
        g = g.sort_values("Start").reset_index(drop=True)
        starts = g["Start"].to_numpy().tolist()
        ends = g["End"].to_numpy().tolist()
        chrom, strand = g["Chromosome"].iat[0], g["Strand"].iat[0]
        m_start, m_end = [starts[0]], [ends[0]]
        for s, e in zip(starts[1:], ends[1:]):
            if s <= m_end[-1]:
                m_end[-1] = max(m_end[-1], e)
            else:
                m_start.append(s); m_end.append(e)
        out.append(pd.DataFrame({
            "Chromosome": chrom, "Start": m_start, "End": m_end, "Strand": strand,
            "transcript_id": tid,
        }))
    return pd.concat(out, ignore_index=True)

def build_transcript_table(rebuild=False):
    """Per-APPRIS-transcript whole-transcript (UTR5+CDS+UTR3) coordinate table.

    Returns {"table": {tid: {...}}, "base2ver": {base_ENST: versioned_ENST}}.
    """
    if rebuild and TRANSCRIPT_TABLE_CACHE.exists():
        TRANSCRIPT_TABLE_CACHE.unlink()
    return fc.config.cached_frame(
        TRANSCRIPT_TABLE_CACHE,
        fc.config.annotation_fingerprint(["build_transcript_table/1"]),
        _build_transcript_table)

def _build_transcript_table():
    cds_df = fc.config.load_annotation()
    utr_df = fc.config.load_appris_utr()
    meta_df = fc.config.load_appris_meta()
    body_df = fc.config.load_gene_bodies()

    gene_id_map = cds_df.drop_duplicates("transcript_id").set_index("transcript_id")["gene_id"]
    body_lookup = body_df.set_index("transcript_id")

    cds_part = cds_df[["Chromosome", "Start", "End", "Strand", "transcript_id"]].copy()
    utr_part = utr_df[["Chromosome", "Start", "End", "Strand", "transcript_id"]].copy()
    exons = pd.concat([cds_part, utr_part], ignore_index=True)
    exons = _merge_contiguous(exons)
    exons["exon_len"] = exons["End"] - exons["Start"]
    exons["order_key"] = np.where(exons["Strand"] == "+", exons["Start"], -exons["Start"])
    exons = exons.sort_values(["transcript_id", "order_key"]).reset_index(drop=True)
    exons["cum_start"] = exons.groupby("transcript_id", sort=False)["exon_len"].cumsum() - exons["exon_len"]

    table = {}
    for tid, g in exons.groupby("transcript_id", sort=False):
        chrom = g["Chromosome"].iat[0]
        strand = g["Strand"].iat[0]
        if tid in body_lookup.index:
            body_start = int(body_lookup.loc[tid, "Start"])
            body_end = int(body_lookup.loc[tid, "End"])
        else:
            body_start, body_end = int(g["Start"].min()), int(g["End"].max())
        table[tid] = dict(
            gene_id=gene_id_map.get(tid, ""),
            chrom=chrom, strand=strand,
            body_start=body_start, body_end=body_end,
            total_len=int(g["exon_len"].sum()),
            cum_start=g["cum_start"].to_numpy(),
            g_start=g["Start"].to_numpy(),
            g_end=g["End"].to_numpy(),
        )

    base2ver = {tid.split(".", 1)[0]: tid for tid in meta_df["transcript_id"]}
    return {"table": table, "base2ver": base2ver}

def build_exon_gene_table(rebuild=False):
    """ALL GTF exon intervals (every gene_type) with gene_id + gene_type kept.

    Strand is dropped because buckets 5-7 ask about ANY other gene's exon, strandless.
    """
    if rebuild and EXON_GENE_CACHE.exists():
        EXON_GENE_CACHE.unlink()
    return fc.config.cached_frame(
        EXON_GENE_CACHE,
        fc.config.annotation_fingerprint(["build_exon_gene_table/1"]),
        _build_exon_gene_table)

def _build_exon_gene_table():
    gtf = fc.config.gtf_path()
    gid_re = re.compile(r'gene_id "([^"]+)"')
    gtype_re = re.compile(r'gene_type "([^"]+)"')
    rows = []
    _open = gzip.open if gtf.endswith(".gz") else open
    with _open(gtf, "rt") as fh:
        for line in fh:
            if line[0] == "#":
                continue
            f = line.rstrip("\n").split("\t")
            if len(f) < 9 or f[2] != "exon":
                continue
            m_gid = gid_re.search(f[8])
            m_gt = gtype_re.search(f[8])
            rows.append({
                "Chromosome": f[0], "Start": int(f[3]) - 1, "End": int(f[4]),
                "gene_id": m_gid.group(1) if m_gid else "",
                "gene_type": m_gt.group(1) if m_gt else "",
            })
    return pd.DataFrame(rows)

def load_exon_gene_pr(rebuild=False):
    import pyranges as pr
    df = build_exon_gene_table(rebuild=rebuild)
    return pr.PyRanges(df.reset_index(drop=True))


# ── gene bodies with biotype  the tie tests annotation side ───────────────────

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
        exons = build_exon_gene_table()
        return exons.groupby("gene_id", sort=False).agg(
            Chromosome=("Chromosome", "first"),
            Start=("Start", "min"),
            End=("End", "max"),
            gene_type=("gene_type", "first"),
        ).reset_index()

    cache = fc.output_root() / ".cache" / "read_categories" / "gene_body.pkl"
    if rebuild and cache.exists():
        cache.unlink()
    frame = fc.config.cached_frame(
        cache,
        fc.config.annotation_fingerprint(["gene_body_pr/1"]),
        build)
    return pr.PyRanges(frame)
