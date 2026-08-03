#!/usr/bin/env python3
"""Sample inputs shared by the BAM-reading analyses: BAM discovery, the CDS-exon coordinate table, and the offset-controlled read-length set."""
from __future__ import annotations

import os
import re
import sys
from pathlib import Path

import numpy as np
import pandas as pd

_HERE = Path(__file__).resolve().parent
REPO = _HERE.parents[1]

for _entry in (str(_HERE), str(_HERE / "ribo_seq_qc")):
    if _entry not in sys.path:
        sys.path.insert(0, _entry)
import config

GENOME_BAM_TEMPLATE = "{s}/genome/alignment_ribo/merged/{s}.post_dedup.bam"
TXOME_BAM_TEMPLATE = (
    "{s}/transcriptome/alignment_ribo/merged/{s}.transcriptome.post_dedup.bam")
RNA_GENOME_BAM_TEMPLATE = "{s}/rnaseq/genome/alignment_ribo/merged/{s}.rnaseq.post_dedup.bam"
RNA_TXOME_BAM_TEMPLATE = ("{s}/rnaseq/transcriptome/alignment_ribo/merged/"
                          "{s}.rnaseq.transcriptome.post_dedup.bam")

#: Transcriptome BAMs carry NO NH tag (verified). Bowtie2 MAPQ maxes at 42 for a confident
#: stripped, so MAPQ == 42 is the operational "uniquely mapped" set (~94 % of reads).
DEFAULT_TXOME_MIN_MAPQ = 42

_CDS_HEADER = re.compile(r"\|CDS:(\d+)-(\d+)\|")

class InputError(RuntimeError):
    pass

def bams_root(required=True):
    """The RiboFlow output tree. No default: an absent tree is an error, not a guess."""
    value = os.environ.get("RIBOFLOW_PAPER_BAMS")
    if not value:
        if not required:
            return None
        raise InputError(
            "no BAM tree configured. Set RIBOFLOW_PAPER_BAMS to a RiboFlow-genome output "
            "directory, or pass --bams to code/make_tables.py, which sets it for you.\n"
            "Expected layout, per sample:\n"
            "    <root>/" + GENOME_BAM_TEMPLATE + "\n"
            "    <root>/" + TXOME_BAM_TEMPLATE)
    root = Path(value)
    if required and not root.is_dir():
        raise InputError("RIBOFLOW_PAPER_BAMS is not a directory: %s" % root)
    return root

def output_root() -> Path:
    """Where analyses write. Defaults to `results/` inside this repository."""
    return Path(os.environ.get("RIBOFLOW_PAPER_OUT", REPO / "results"))

def txome_min_mapq() -> int:
    return int(os.environ.get("RIBOFLOW_PAPER_TXOME_MIN_MAPQ", DEFAULT_TXOME_MIN_MAPQ))

# ── the one uniqueness policy ────────────────────────────────────────────────
# is about multimappers and deliberately does NOT use these -- it classifies with its own
# GENOME is `NH == 1`, read from the tag, not inferred from MAPQ. STAR happens to encode
# uniqueness as MAPQ 255 and the two agree exactly on these BAMs (0 disagreements over
# 4M primaries), but MAPQ is an aligner convention and `NH` is the count itself. A genome
# BAM with no `NH` is an input error rather than a silent fallback to a MAPQ heuristic:
# guessing the uniqueness rule is how a multimapper quietly becomes a unique read.

def is_unique_genome_read(read) -> bool:
    """A primary, uniquely-mapping genome alignment: `NH == 1`.

    Raises `InputError` if the record has no `NH` tag. The declared STAR inputs all carry
    one; a BAM that does not is a different input than this pipeline documents.
    """
    if read.is_unmapped or read.is_secondary or read.is_supplementary:
        return False
    try:
        return int(read.get_tag("NH")) == 1
    except KeyError:
        raise InputError(
            "genome alignment %r carries no NH tag, so uniqueness cannot be determined. "
            "This pipeline defines a unique genome read as NH == 1 and will not fall back "
            "to a MAPQ heuristic. Re-align with STAR (which always emits NH) or supply the "
            "BAM this cohort documents." % (read.query_name,))

def is_unique_txome_read(read) -> bool:
    """A primary transcriptome alignment at or above the Bowtie2 confident-unique ceiling."""
    if read.is_unmapped or read.is_secondary or read.is_supplementary:
        return False
    return read.mapping_quality >= txome_min_mapq()

def _template(name, default):
    return os.environ.get(name, default)

def genome_bam(sample: str) -> Path:
    return bams_root() / _template(
        "RIBOFLOW_PAPER_GENOME_BAM_TPL", GENOME_BAM_TEMPLATE).format(s=sample)

def txome_bam(sample: str) -> Path:
    return bams_root() / _template(
        "RIBOFLOW_PAPER_TXOME_BAM_TPL", TXOME_BAM_TEMPLATE).format(s=sample)

def rna_genome_bam(sample: str) -> Path:
    return bams_root() / _template(
        "RIBOFLOW_PAPER_RNA_GENOME_BAM_TPL", RNA_GENOME_BAM_TEMPLATE).format(s=sample)

def rna_txome_bam(sample: str) -> Path:
    return bams_root() / _template(
        "RIBOFLOW_PAPER_RNA_TXOME_BAM_TPL", RNA_TXOME_BAM_TEMPLATE).format(s=sample)

def discover_samples() -> list:
    """Samples with BOTH a genome and a transcriptome ribo BAM present."""
    root = bams_root()
    return [d.name for d in sorted(root.iterdir())
            if d.is_dir() and genome_bam(d.name).exists() and txome_bam(d.name).exists()]

def require_bams(samples) -> None:
    """Fail before any compute starts, naming every missing BAM at once."""
    missing = []
    for sample in samples:
        for label, path in (("genome", genome_bam(sample)), ("txome", txome_bam(sample))):
            if not path.exists():
                missing.append("  %-10s %-8s %s" % (sample, label, path))
    if missing:
        raise InputError("these BAMs do not exist:\n" + "\n".join(missing))

def build_cds_table() -> dict:
    """Per-transcript CDS exon table plus lookups.

    Returns a dict with:
      cds          CDS exons 5'->3' with exon_index, exon_len, cds_cum_start/end
                   (CDS-relative nt from the start codon)
      tx_cumstarts {tid: array of cds_cum_start}
      cds_total    {tid: genomic CDS length in nt, stop codon excluded}
      n_exons      {tid: number of CDS exons}
      gene_name / gene_id / strand / chrom   {tid: value}
      junctions    {tid: array of internal CDS-relative exon boundaries}
    """
    cds = config.load_annotation()[
        ["Chromosome", "Start", "End", "Strand", "Phase",
         "transcript_id", "gene_id", "gene_name"]
    ].copy()
    cds["exon_len"] = cds["End"] - cds["Start"]
    cds["order_key"] = np.where(cds["Strand"] == "+", cds["Start"], -cds["Start"])
    cds = cds.sort_values(["transcript_id", "order_key"]).reset_index(drop=True)
    grp = cds.groupby("transcript_id", sort=False)
    cds["exon_index"] = grp.cumcount()
    cds["cds_cum_start"] = grp["exon_len"].cumsum() - cds["exon_len"]
    cds["cds_cum_end"] = cds["cds_cum_start"] + cds["exon_len"]

    tx_cumstarts = {tid: sub["cds_cum_start"].to_numpy()
                    for tid, sub in cds.groupby("transcript_id", sort=False)}
    # internal boundaries = cds_cum_start of every exon after the first: the splice
    junctions = {tid: arr[1:].astype(np.int64) for tid, arr in tx_cumstarts.items()}

    return {
        "cds": cds,
        "tx_cumstarts": tx_cumstarts,
        "cds_total": grp["exon_len"].sum().astype(int).to_dict(),
        "n_exons": grp.size().to_dict(),
        "gene_name": grp["gene_name"].first().to_dict(),
        "gene_id": grp["gene_id"].first().to_dict(),
        "strand": grp["Strand"].first().to_dict(),
        "chrom": grp["Chromosome"].first().to_dict(),
        "junctions": junctions,
    }

def _as_bool(series):
    return series.map(lambda x: str(x).strip().lower() in ("true", "1"))

def genome_qc_path() -> Path:
    return Path(config.out_dir()) / "tables" / "readlen_window_qc.csv"

def txome_qc_path() -> Path:
    return Path(config.tx_out_dir()) / "tables" / "readlen_window_qc.csv"
