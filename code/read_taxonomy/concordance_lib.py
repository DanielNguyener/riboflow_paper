#!/usr/bin/env python3
"""Read-ID-level genome<->transcriptome alignment concordance."""
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

OUTDIR = fc.output_root() / ".cache" / "read_taxonomy" / "alignment_concordance"
CACHE_DIR = fc.output_root() / ".cache" / "read_taxonomy"
EXON_GENE_CACHE = CACHE_DIR / "exon_gene_table.pkl"
TRANSCRIPT_TABLE_CACHE = CACHE_DIR / "transcript_coord_table.pkl"

COORD_TOL = 0

CATEGORIES = [
    "concordant",
    "same_gene_coord_discordant",
    "different_pc_gene",
    "pseudogene",
    "other_biotype",
    "intron_or_intergenic",
    "strand_discordant",
    "splice_discordant",
]

CATEGORY_LABEL = {
    "concordant": "concordant at expected gene and locus",
    "same_gene_coord_discordant": "same gene, coordinate-discordant",
    "different_pc_gene": "different protein-coding gene",
    "pseudogene": "assigned to a pseudogene",
    "other_biotype": "another biotype",
    "intron_or_intergenic": "intronic or intergenic",
    "strand_discordant": "strand-discordant",
    "splice_discordant": "splice-structure-discordant",
}

def read_txome_primary(bam_path, base2ver):
    """One pass over every primary txome alignment (any MAPQ) -> (present, unique_qnames,
    all_qnames): {qname: (tid, tx_pos, tx_len)} for APPRIS-resolving primaries, the
    MAPQ>=TXOME_MIN_MAPQ subset, and every primary qname (resolved or not).
    """
    present = {}
    unique_qnames = set()
    all_qnames = set()
    bam = pysam.AlignmentFile(str(bam_path), "rb")
    for r in bam.fetch(until_eof=True):
        if r.is_unmapped or r.is_secondary or r.is_supplementary:
            continue
        all_qnames.add(r.query_name)
        base = r.reference_name.split(".", 1)[0].split("|", 1)[0]
        tid = base2ver.get(base)
        if tid is None:
            continue
        present[r.query_name] = (tid, r.reference_start, r.reference_length)
        if r.mapping_quality >= fc.txome_min_mapq():
            unique_qnames.add(r.query_name)
    bam.close()
    return present, unique_qnames, all_qnames

def read_genome_unique(bam_path):
    """qname -> (chrom, strand, pos5, n_blocks, blk_min, blk_max) for primary, unique genome reads.

    Unique = NH==1; pos5 is the 5'-most genomic base (leftmost fwd / rightmost-1 rev).
    """
    out = {}
    bam = pysam.AlignmentFile(str(bam_path), "rb")
    for r in bam.fetch(until_eof=True):
        if r.is_unmapped or r.is_secondary or r.is_supplementary:
            continue
        uniq = fc.is_unique_genome_read(r)
        if not uniq:
            continue
        blocks = r.get_blocks()
        if not blocks:
            continue
        strand = "-" if r.is_reverse else "+"
        pos5 = blocks[0][0] if strand == "+" else blocks[-1][1] - 1
        out[r.query_name] = (r.reference_name, strand, pos5, len(blocks),
                              blocks[0][0], blocks[-1][1])
    bam.close()
    return out

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

    Strand is dropped deliberately: buckets 5-7 ask about ANY other gene's exon, strandless.
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

def classify_all(shared_qnames, txome_dict, genome_dict, transcript_payload, exon_gene_pr):
    """Return a DataFrame (qname, tid, label) for every shared read."""
    table = transcript_payload["table"]

    rows = [(q,) + txome_dict[q] for q in shared_qnames]
    df = pd.DataFrame(rows, columns=["qname", "tid", "tx_pos", "tx_len"])
    grows = [(q,) + genome_dict[q] for q in df["qname"]]
    gdf = pd.DataFrame(grows, columns=["qname", "g_chrom", "g_strand", "g_pos5",
                                        "g_nblocks", "g_min", "g_max"])
    df = df.merge(gdf, on="qname", how="left")

    label = np.full(len(df), "intron_or_intergenic", dtype=object)
    outside_mask = np.ones(len(df), dtype=bool)

    for tid, idx in df.groupby("tid", sort=False).groups.items():
        t = table.get(tid)
        if t is None:
            continue
        pos = df.index.get_indexer(idx)
        sub = df.loc[idx]

        cum_start = t["cum_start"]
        g_start = t["g_start"]
        g_end = t["g_end"]
        n_exons = len(cum_start)

        tx_pos = sub["tx_pos"].to_numpy()
        tx_len = sub["tx_len"].to_numpy()
        idx_start = np.clip(np.searchsorted(cum_start, tx_pos, side="right") - 1, 0, n_exons - 1)
        idx_end = np.clip(np.searchsorted(cum_start, tx_pos + tx_len - 1, side="right") - 1, 0, n_exons - 1)
        expected_njunc = idx_end - idx_start
        offset = tx_pos - cum_start[idx_start]
        if t["strand"] == "+":
            expected_pos5 = g_start[idx_start] + offset
        else:
            expected_pos5 = g_end[idx_start] - 1 - offset

        obs_chrom = sub["g_chrom"].to_numpy()
        obs_strand = sub["g_strand"].to_numpy()
        obs_pos5 = sub["g_pos5"].to_numpy(dtype=np.int64)
        obs_min = sub["g_min"].to_numpy()
        obs_max = sub["g_max"].to_numpy()
        obs_njunc = sub["g_nblocks"].to_numpy() - 1

        same_chrom = obs_chrom == t["chrom"]
        within_body = same_chrom & (obs_min < t["body_end"]) & (obs_max > t["body_start"])

        lab = np.full(len(sub), "", dtype=object)
        strand_bad = within_body & (obs_strand != t["strand"])
        lab[strand_bad] = "strand_discordant"
        splice_bad = within_body & ~strand_bad & (obs_njunc != expected_njunc)
        lab[splice_bad] = "splice_discordant"
        coord_bad = (within_body & ~strand_bad & ~splice_bad &
                     (np.abs(obs_pos5 - expected_pos5.astype(np.int64)) > COORD_TOL))
        lab[coord_bad] = "same_gene_coord_discordant"
        concordant = within_body & ~strand_bad & ~splice_bad & ~coord_bad
        lab[concordant] = "concordant"

        label[pos[within_body]] = lab[within_body]
        outside_mask[pos[within_body]] = False

    df["label"] = label

    if outside_mask.any():
        out_labels = _classify_outside(df.loc[outside_mask], exon_gene_pr)
        df.loc[outside_mask, "label"] = df.loc[outside_mask, "qname"].map(out_labels).fillna(
            "intron_or_intergenic").to_numpy()

    return df[["qname", "tid", "label"]]

def _classify_outside(sub_df, exon_gene_pr):
    """Priority: different protein-coding gene > pseudogene > other biotype > (no hit)."""
    import pyranges as pr

    reads_df = sub_df[["qname", "g_chrom", "g_min", "g_max"]].rename(
        columns={"g_chrom": "Chromosome", "g_min": "Start", "g_max": "End"})
    reads_pr = pr.PyRanges(reads_df.reset_index(drop=True))
    joined = reads_pr.join(exon_gene_pr, strandedness=False, how=None).df
    if joined.empty or "gene_type" not in joined.columns:
        return pd.Series(dtype=object)

    def pick(types):
        if any(t == "protein_coding" for t in types):
            return "different_pc_gene"
        if any("pseudogene" in t for t in types):
            return "pseudogene"
        return "other_biotype"

    return joined.groupby("qname")["gene_type"].apply(lambda s: pick(set(s)))

def load_exon_gene_pr(rebuild=False):
    import pyranges as pr
    df = build_exon_gene_table(rebuild=rebuild)
    return pr.PyRanges(df.reset_index(drop=True))
