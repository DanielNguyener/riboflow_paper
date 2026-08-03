#!/usr/bin/env python3
"""CDS-assigned read counts per APPRIS transcript, on all four BAMs of a sample."""
from __future__ import annotations

import sys
from pathlib import Path

import numpy as np
import pandas as pd
import pyranges as pr
import pysam
from scipy.stats import pearsonr, spearmanr

_HERE = Path(__file__).resolve().parent
_COMMON = _HERE.parent / "common"
_COVERAGE = _HERE.parent / "coverage"
for _entry in (str(_HERE), str(_COMMON), str(_COMMON / "ribo_seq_qc"), str(_COVERAGE)):
    if _entry not in sys.path:
        sys.path.insert(0, _entry)
import bam_inputs as fc

rna_genome_bam = fc.rna_genome_bam
rna_txome_bam = fc.rna_txome_bam

class RiboRnaError(RuntimeError):
    pass

def load_bundle(gtf, appris, regions=None, cache_path=None, left_span=35, right_span=10):
    """The shared coordinate bundle -- the SAME one the coverage builder uses.

    Reusing `annotation_cache` rather than re-deriving the annotation is the point: a
    second parse could disagree with the coverage figure about what a CDS is.
    """
    import annotation_cache as ac

    return ac.load_or_build(cache_path, gtf, appris, regions, left_span, right_span)

def transcript_ids_in_reference(txome_bam) -> set:
    """Versioned transcript ids named by a transcriptome BAM's `@SQ` records.

    The reference name is split at the first `|` and nothing else is read out of it. The
    CDS comes from the annotation; the reference only says which transcript this is.
    """
    handle = pysam.AlignmentFile(str(txome_bam), "rb")
    try:
        names = list(handle.references)
    finally:
        handle.close()
    ids = [name.split("|", 1)[0] for name in names]
    if len(set(ids)) != len(ids):
        seen, duplicates = set(), set()
        for tid in ids:
            (duplicates if tid in seen else seen).add(tid)
        raise RiboRnaError(
            "the transcriptome reference names %d transcript(s) more than once (e.g. %s). "
            "A transcript id must select exactly one reference sequence."
            % (len(duplicates), sorted(duplicates)[:3]))
    return set(ids)

def build_universe(bundle, txome_bam):
    """The transcripts this comparison can count on BOTH routes, and why any are missing.

    Returns `(universe, report)`. `universe` is the ordered list of VERSIONED transcript
    ids present in the annotation bundle and named by the transcriptome reference; every
    one has exactly one canonical CDS. `report` records the size and each side's
    shortfall, so a universe that quietly shrinks is visible rather than absorbed into the
    correlation.

    No transcript is excluded for lacking a UTR. A CDS that starts at transcript
    coordinate 0 is a transcript with no annotated 5'UTR, not a defective one, and it has
    the same countable CDS as any other.
    """
    annotated = list(bundle["transcripts"]["transcript_id"])
    if len(set(annotated)) != len(annotated):
        raise RiboRnaError("the annotation bundle repeats a transcript id")
    in_reference = transcript_ids_in_reference(txome_bam)

    universe = sorted(t for t in annotated if t in in_reference)
    report = {
        "n_universe": len(universe),
        "n_annotated": len(annotated),
        "n_annotated_not_in_reference": len(set(annotated) - in_reference),
        "n_reference_not_annotated": len(in_reference - set(annotated)),
    }

    genes = bundle["transcripts"].set_index("transcript_id").loc[universe, "gene_id"]
    if genes.nunique() != len(universe):
        repeated = genes[genes.duplicated(keep=False)]
        raise RiboRnaError(
            "APPRIS must select at most one transcript per gene, but %d transcripts share "
            "%d gene ids (e.g. %s). Excluding cross-gene overlaps assumes one transcript "
            "per gene and would be wrong otherwise."
            % (len(repeated), repeated.nunique(), sorted(repeated.index)[:4]))
    return universe, report

def transcript_cds_spans(bundle, universe) -> dict:
    """transcript_id -> (start, end), the canonical CDS in TRANSCRIPT coordinates."""
    regions = bundle["regions"]
    cds = regions[regions["label"] == "CDS"]
    counts = cds["transcript_id"].value_counts()
    repeated = counts[counts > 1]
    if len(repeated):
        raise RiboRnaError(
            "%d transcripts carry more than one canonical CDS region (e.g. %s)"
            % (len(repeated), list(repeated.index[:3])))
    spans = dict(zip(cds["transcript_id"], zip(cds["start"], cds["end"])))
    missing = [t for t in universe if t not in spans]
    if missing:
        raise RiboRnaError("%d transcripts in the universe have no canonical CDS (e.g. %s)"
                           % (len(missing), missing[:3]))
    return {t: (int(spans[t][0]), int(spans[t][1])) for t in universe}

def genome_cds_intervals(bundle, universe):
    """The same canonical CDS as exon-aware GENOMIC segments, as a PyRanges.

    Asserts that each transcript's genomic CDS length equals its transcript-coordinate CDS
    length. The two are built by different code paths from the same annotation; if they
    ever disagreed, a genome count and a transcriptome count of the same transcript would
    be measuring different intervals and the route comparison would be meaningless.
    """
    wanted = set(universe)
    table = bundle["cds_table"]
    table = table[table["transcript_id"].isin(wanted)]

    genomic_length = table.groupby("transcript_id")["exon_len"].sum()
    spans = transcript_cds_spans(bundle, universe)
    mismatched = [t for t in universe
                  if int(genomic_length.get(t, -1)) != spans[t][1] - spans[t][0]]
    if mismatched:
        raise RiboRnaError(
            "%d transcripts have a genomic CDS length that differs from their "
            "transcript-coordinate CDS length (e.g. %s). The two routes would then be "
            "counting different intervals." % (len(mismatched), mismatched[:3]))

    frame = table[["Chromosome", "Start", "End", "Strand", "transcript_id"]].copy()
    return pr.PyRanges(frame)

def assert_single_end(bam_path, limit=200000):
    """Refuse a paired BAM.

    The counting rule is one count per retained alignment. That is correct for the
    single-end RNA-seq this cohort uses and WRONG for paired data, where it would count a
    fragment twice. There is no mate-aware branch here on purpose: a paired library needs
    a deliberate decision about which end represents the fragment, and silently guessing
    is how a doubled count reaches a figure.
    """
    handle = pysam.AlignmentFile(str(bam_path), "rb")
    try:
        for i, read in enumerate(handle.fetch(until_eof=True)):
            if read.is_paired:
                raise RiboRnaError(
                    "%s contains paired-end alignments (first at read %r). This pipeline "
                    "counts one read per alignment, which double-counts a paired "
                    "fragment. The cohort's RNA-seq is single-end; a paired library needs "
                    "an explicit fragment rule, not this one."
                    % (Path(bam_path).name, read.query_name))
            if limit and i >= limit:
                break
    finally:
        handle.close()

def _five_prime(read):
    """The reference-aligned 5' nucleotide: `reference_start` forward, `end - 1` reverse."""
    return (read.reference_end - 1) if read.is_reverse else read.reference_start

def count_genome_cds(bam_path, cds_pr, stranded: bool, read_lengths=None):
    """CDS-assigned counts from a GENOME BAM, keyed by versioned transcript id.

    `read_lengths` restricts to a set of selected read lengths (ribo) or is None (rna).

    `stranded=True` (ribo) joins same-strand: ribo footprints are reliably sense.
    `stranded=False` (rna) joins strand-agnostically, because this cohort's RNA-seq library
    strandedness is MIXED across studies (HEK293 is reverse-stranded). A same-strand rule
    would silently zero the antisense-protocol libraries and nothing downstream would show
    it.

    Returns `(counts, n_assigned, n_ambiguous_excluded, n_retained)`, where `n_retained` is
    every read that passed the filters regardless of where it landed -- the denominator of
    the CDS-capture fraction.
    """
    chrom, pos5, strand = [], [], []
    handle = pysam.AlignmentFile(str(bam_path), "rb")
    try:
        for read in handle.fetch(until_eof=True):
            if not fc.is_unique_genome_read(read):
                continue
            if read_lengths is not None and read.query_length not in read_lengths:
                continue
            chrom.append(read.reference_name)
            pos5.append(_five_prime(read))
            strand.append("-" if read.is_reverse else "+")
    finally:
        handle.close()

    n_retained = len(chrom)
    if n_retained == 0:
        return {}, 0, 0, 0

    positions = np.asarray(pos5, dtype=np.int64)
    reads = pr.PyRanges(pd.DataFrame({
        "Chromosome": chrom, "Start": positions, "End": positions + 1,
        "Strand": strand, "read_idx": np.arange(n_retained, dtype=np.int64),
    }))
    frame = reads.join(cds_pr, strandedness="same" if stranded else False).df
    if frame.empty:
        return {}, 0, 0, n_retained
    counts, n_assigned, n_ambiguous = _resolve_overlaps(frame)
    return counts, n_assigned, n_ambiguous, n_retained

def _resolve_overlaps(frame):
    """(counts, n_assigned, n_ambiguous_excluded) from a read x CDS-exon join.

    A read that lands in two CDS exons of the SAME transcript is one assignment, so
    (read, transcript) pairs are collapsed first. A read still left with more than one
    transcript is in a cross-gene CDS overlap: it is dropped and tallied. Both steps are
    set operations, so the result cannot depend on the order rows came out of the join.
    """
    pairs = frame[["read_idx", "transcript_id"]].drop_duplicates()
    per_read = pairs.groupby("read_idx", sort=False)["transcript_id"].size()
    ambiguous = set(per_read.index[per_read > 1])
    n_ambiguous = len(ambiguous)
    if n_ambiguous:
        pairs = pairs[~pairs["read_idx"].isin(ambiguous)]
    counts = pairs["transcript_id"].value_counts()
    return {str(k): int(v) for k, v in counts.items()}, int(counts.sum()), n_ambiguous

def count_txome_cds(bam_path, cds_spans, read_lengths=None):
    """CDS-assigned counts from a TRANSCRIPTOME BAM, keyed by versioned transcript id.

    One reference means one transcript, so no ambiguity is possible and none is reported.
    The 5' nucleotide is taken the same way as on the genome route: bowtie2 `--norc` makes
    ribo alignments forward, but the RNA-seq transcriptome BAMs do contain reverse
    alignments, and those are placed by `reference_end - 1` rather than assumed forward.

    Returns `(counts, n_assigned, n_retained)`.
    """
    counts = {}
    n_assigned = 0
    n_retained = 0
    handle = pysam.AlignmentFile(str(bam_path), "rb")
    try:
        tid_of_ref = [name.split("|", 1)[0] for name in handle.references]
        span_of_ref = [cds_spans.get(t) for t in tid_of_ref]
        for read in handle.fetch(until_eof=True):
            if not fc.is_unique_txome_read(read):
                continue
            if read_lengths is not None and read.query_length not in read_lengths:
                continue
            n_retained += 1
            span = span_of_ref[read.reference_id]
            if span is None:
                continue
            position = _five_prime(read)
            if span[0] <= position < span[1]:
                tid = tid_of_ref[read.reference_id]
                counts[tid] = counts.get(tid, 0) + 1
                n_assigned += 1
    finally:
        handle.close()
    return counts, n_assigned, n_retained

def correlate(ribo_counts: dict, rna_counts: dict, universe: list) -> dict:
    """Ribo-vs-RNA agreement over the whole universe, zeros kept.

    Two statistics, both on RAW counts:

      spearman_rho       rank-based, so invariant to any monotonic per-assay transform;
      pearson_log2_raw   Pearson on `log2(count + 1)`, un-normalised.

    There is no CPM here. Dividing each route by its own library size (or by its own
    within-universe sum) rescales away precisely the quantity the route comparison is
    about -- how much of the library each reference can place in a CDS -- so a normalised
    correlation answers a different question than the one asked. That capture fraction is
    reported separately, as `ribo_cds_frac` / `rna_cds_frac`.
    """
    rb = np.array([ribo_counts.get(t, 0) for t in universe], dtype=np.float64)
    rn = np.array([rna_counts.get(t, 0) for t in universe], dtype=np.float64)
    if rb.sum() == 0 or rn.sum() == 0:
        return {"n_transcripts": len(universe),
                "pearson_log2_raw": float("nan"), "spearman_rho": float("nan")}
    log_rb, log_rn = np.log2(rb + 1.0), np.log2(rn + 1.0)
    pearson = (float(pearsonr(log_rb, log_rn)[0])
               if log_rb.std() > 0 and log_rn.std() > 0 else float("nan"))
    rho = (float(spearmanr(rb, rn).correlation)
           if rb.std() > 0 and rn.std() > 0 else float("nan"))
    return {"n_transcripts": len(universe),
            "pearson_log2_raw": pearson, "spearman_rho": rho}
