#!/usr/bin/env python3
"""Build one sample's shared-coordinate coverage HDF5 from its two ribo BAMs."""
from __future__ import annotations

import argparse
import json
import re
import sys
import time
from pathlib import Path

import numpy as np
import pandas as pd

sys.path.insert(0, str(Path(__file__).resolve().parent.parent / "common"))
import bam_inputs                      # the one uniqueness policy

_CDS_HEADER = re.compile(r"\|CDS:(\d+)-(\d+)\|")

class BuildError(RuntimeError):
    pass

def log(message):
    print("[coverage] %s" % message, flush=True)

def _exon_pyranges(exons, transcripts):
    """The full exon map as PyRanges, carrying what placement needs."""
    import pyranges as pr

    strand = transcripts["strand"].to_numpy()[exons["transcript_index"].to_numpy()]
    return pr.PyRanges(pd.DataFrame({
        "Chromosome": exons["chrom"].to_numpy(),
        "Start": exons["g_start"].to_numpy(),
        "End": exons["g_end"].to_numpy(),
        "Strand": strand,
        "tx_index": exons["transcript_index"].to_numpy(),
        "ex_tx_start": exons["tx_start"].to_numpy(),
        "ex_g_start": exons["g_start"].to_numpy(),
        "ex_g_end": exons["g_end"].to_numpy(),
    }))

def _cds_pyranges(cds_table):
    import pyranges as pr
    return pr.PyRanges(cds_table[["Chromosome", "Start", "End", "Strand",
                                  "transcript_id", "cds_cum_start"]])

def read_genome_psites(bam_path, offsets):
    """Stream the genome BAM -> per-read (chrom, psite, strand).

    Uniquely-mapping primaries (`bam_inputs.is_unique_genome_read`, NH == 1) whose read
    length is in `offsets`.
    Placement is CIGAR-aware (`psite_placement.place`), so a junction-spanning read's
    P-site lands on a base the read covers instead of inside the intron.
    """
    import pysam
    import psite_placement

    chroms, positions, strands = [], [], []
    bam = pysam.AlignmentFile(str(bam_path), "rb")
    try:
        for read in bam.fetch(until_eof=True):
            if not bam_inputs.is_unique_genome_read(read):
                continue
            offset = offsets.get(read.query_length)
            if offset is None:
                continue
            position = psite_placement.place(read, offset)
            if position is None:
                continue
            chroms.append(read.reference_name)
            positions.append(position)
            strands.append("-" if read.is_reverse else "+")
    finally:
        bam.close()
    return chroms, np.asarray(positions, dtype=np.int64), strands

def read_genome_blocks(bam_path, lengths):
    """Stream the genome BAM -> aligned blocks. `get_blocks()` splits on N, so introns are
    already excluded and each block is contiguous on the reference."""
    import pysam

    chroms, starts, ends, strands, read_ids = [], [], [], [], []
    index = 0
    bam = pysam.AlignmentFile(str(bam_path), "rb")
    try:
        for read in bam.fetch(until_eof=True):
            if not bam_inputs.is_unique_genome_read(read):
                continue
            if read.query_length not in lengths:
                continue
            strand = "-" if read.is_reverse else "+"
            for block_start, block_end in read.get_blocks():
                chroms.append(read.reference_name)
                starts.append(block_start)
                ends.append(block_end)
                strands.append(strand)
                read_ids.append(index)
            index += 1
    finally:
        bam.close()
    return (chroms, np.asarray(starts, dtype=np.int64), np.asarray(ends, dtype=np.int64),
            strands, np.asarray(read_ids, dtype=np.int64), index)

def _stage1_psite_assignment(chroms, positions, strands, cds_pr, cds_total_by_id, tx_index_of_id):
    """The P-site rule: FIRST CDS-exon overlap, clipped to [0, cds_total)."""
    import pyranges as pr

    reads = pr.PyRanges(pd.DataFrame({
        "Chromosome": chroms, "Start": positions, "End": positions + 1,
        "Strand": strands, "read_idx": np.arange(len(chroms), dtype=np.int64)}))
    joined = reads.join(cds_pr, strandedness="same").df
    if joined.empty:
        return np.empty(0, dtype=np.int64), np.empty(0, dtype=np.int64)
    joined = joined.drop_duplicates("read_idx", keep="first")
    plus = (joined["Strand"] == "+").to_numpy()
    within = np.where(plus,
                      joined["Start"].to_numpy() - joined["Start_b"].to_numpy(),
                      joined["End_b"].to_numpy() - 1 - joined["Start"].to_numpy())
    cds_rel = joined["cds_cum_start"].to_numpy() + within
    total = joined["transcript_id"].map(cds_total_by_id).to_numpy()
    keep = (cds_rel >= 0) & (cds_rel < total)
    joined = joined[keep]
    return (joined["read_idx"].to_numpy(),
            joined["transcript_id"].map(tx_index_of_id).to_numpy(dtype=np.int64))

def project_genome_psites(chroms, positions, strands, exon_pr, cds_pr,
                          cds_total_by_id, tx_index_of_id, coverage_offset):
    """Genome P-sites -> absolute indices into the concatenated coverage array.

    Returns (indices, stats). Stage 1 sees the CDS exon table alone; stage 2 picks up the
    reads it left, using the full exon map. Splitting them is what keeps CDS values
    independent of the UTR rows: adding UTR rows to the stage-1 join would change what
    `keep="first"` and `idxmax` select, and so could move a published CDS value.
    """
    import pyranges as pr

    n_reads = len(chroms)
    stage1_reads, stage1_tx = _stage1_psite_assignment(
        chroms, positions, strands, cds_pr, cds_total_by_id, tx_index_of_id)

    assigned = np.full(n_reads, -1, dtype=np.int64)
    assigned[stage1_reads] = stage1_tx

    reads = pr.PyRanges(pd.DataFrame({
        "Chromosome": chroms, "Start": positions, "End": positions + 1,
        "Strand": strands, "read_idx": np.arange(n_reads, dtype=np.int64)}))
    joined = reads.join(exon_pr, strandedness="same").df
    if joined.empty:
        return np.empty(0, dtype=np.int64), {
            "n_assigned_stage1": 0, "n_assigned_stage2_utr": 0, "n_unassigned": n_reads}

    read_idx = joined["read_idx"].to_numpy()
    stage2_mask = assigned[read_idx] == -1
    if stage2_mask.any():
        leftovers = joined[stage2_mask].drop_duplicates("read_idx", keep="first")
        assigned[leftovers["read_idx"].to_numpy()] = leftovers["tx_index"].to_numpy()
    n_stage2 = int((assigned != -1).sum() - len(stage1_reads))

    keep = assigned[read_idx] == joined["tx_index"].to_numpy()
    rows = joined[keep]
    if rows.empty:
        return np.empty(0, dtype=np.int64), {
            "n_assigned_stage1": int(len(stage1_reads)),
            "n_assigned_stage2_utr": n_stage2,
            "n_unassigned": int((assigned == -1).sum())}
    rows = rows.drop_duplicates("read_idx", keep="first")

    plus = (rows["Strand"] == "+").to_numpy()
    position = rows["Start"].to_numpy()
    offset_in_exon = np.where(plus,
                              position - rows["ex_g_start"].to_numpy(),
                              rows["ex_g_end"].to_numpy() - 1 - position)
    tx_position = rows["ex_tx_start"].to_numpy() + offset_in_exon
    indices = coverage_offset[rows["tx_index"].to_numpy()] + tx_position
    return indices, {
        "n_assigned_stage1": int(len(stage1_reads)),
        "n_assigned_stage2_utr": n_stage2,
        "n_unassigned": int((assigned == -1).sum()),
    }

def project_genome_footprints(blocks, exon_pr, cds_pr, tx_index_of_id, coverage_offset):
    """Genome footprint blocks -> (start, end) index ranges in the coverage array.

    Stage 1 is the max-CDS-overlap rule against CDS exons alone; stage 2 uses the full exon map for
    reads stage 1 did not claim. Each surviving (block, exon) overlap becomes one interval.
    """
    import pyranges as pr

    chroms, starts, ends, strands, read_ids, n_reads = blocks
    block_pr = pr.PyRanges(pd.DataFrame({
        "Chromosome": chroms, "Start": starts, "End": ends,
        "Strand": strands, "read_idx": read_ids}))

    cds_joined = block_pr.join(cds_pr, strandedness="same").df
    assigned = np.full(n_reads, -1, dtype=np.int64)
    n_stage1 = 0
    if not cds_joined.empty:
        overlap = (np.minimum(cds_joined["End"].to_numpy(), cds_joined["End_b"].to_numpy())
                   - np.maximum(cds_joined["Start"].to_numpy(),
                                cds_joined["Start_b"].to_numpy()))
        cds_joined = cds_joined.assign(olen=overlap)
        cds_joined = cds_joined[cds_joined["olen"] > 0]
        totals = (cds_joined.groupby(["read_idx", "transcript_id"], sort=False)["olen"]
                  .sum().reset_index())
        winner = totals.loc[totals.groupby("read_idx", sort=False)["olen"].idxmax()]
        assigned[winner["read_idx"].to_numpy()] = \
            winner["transcript_id"].map(tx_index_of_id).to_numpy(dtype=np.int64)
        n_stage1 = int(len(winner))

    joined = block_pr.join(exon_pr, strandedness="same").df
    if joined.empty:
        return (np.empty(0, dtype=np.int64), np.empty(0, dtype=np.int64),
                {"n_assigned_stage1": n_stage1, "n_assigned_stage2_utr": 0,
                 "n_unassigned": n_reads - n_stage1})

    overlap_start = np.maximum(joined["Start"].to_numpy(), joined["ex_g_start"].to_numpy())
    overlap_end = np.minimum(joined["End"].to_numpy(), joined["ex_g_end"].to_numpy())
    joined = joined.assign(ostart=overlap_start, oend=overlap_end,
                           olen=overlap_end - overlap_start)
    joined = joined[joined["olen"] > 0]

    read_idx = joined["read_idx"].to_numpy()
    stage2_mask = assigned[read_idx] == -1
    if stage2_mask.any():
        leftovers = joined[stage2_mask]
        totals = (leftovers.groupby(["read_idx", "tx_index"], sort=False)["olen"]
                  .sum().reset_index())
        winner = totals.loc[totals.groupby("read_idx", sort=False)["olen"].idxmax()]
        assigned[winner["read_idx"].to_numpy()] = winner["tx_index"].to_numpy()
    n_stage2 = int((assigned != -1).sum()) - n_stage1

    keep = assigned[joined["read_idx"].to_numpy()] == joined["tx_index"].to_numpy()
    rows = joined[keep]
    if rows.empty:
        return (np.empty(0, dtype=np.int64), np.empty(0, dtype=np.int64),
                {"n_assigned_stage1": n_stage1, "n_assigned_stage2_utr": n_stage2,
                 "n_unassigned": int((assigned == -1).sum())})

    plus = (rows["Strand"] == "+").to_numpy()
    ostart = rows["ostart"].to_numpy()
    oend = rows["oend"].to_numpy()
    ex_tx_start = rows["ex_tx_start"].to_numpy()
    ex_g_start = rows["ex_g_start"].to_numpy()
    ex_g_end = rows["ex_g_end"].to_numpy()
    # on '-' the transcript runs the other way, so the overlap's 5' end is its genomic END
    tx_start = np.where(plus, ex_tx_start + (ostart - ex_g_start),
                        ex_tx_start + (ex_g_end - oend))
    base = coverage_offset[rows["tx_index"].to_numpy()]
    return (base + tx_start, base + tx_start + (oend - ostart),
            {"n_assigned_stage1": n_stage1, "n_assigned_stage2_utr": n_stage2,
             "n_unassigned": int((assigned == -1).sum())})

def txome_reference_map(bam_path, tx_index_of_base, transcript_len):
    """{reference name: (transcript index, reference length)} and a length cross-check.

    The transcriptome reference IS the shared coordinate, so its `@SQ` length must equal
    `transcript_len`. A mismatch means the two routes are not on one ruler.
    """
    import pysam

    mapping, mismatches = {}, []
    bam = pysam.AlignmentFile(str(bam_path), "rb")
    try:
        lengths = dict(zip(bam.references, bam.lengths))
    finally:
        bam.close()
    for name, length in lengths.items():
        if not _CDS_HEADER.search(name):
            continue
        base = name.split("|", 1)[0].split(".", 1)[0]
        index = tx_index_of_base.get(base)
        if index is None:
            continue
        if int(length) != int(transcript_len[index]):
            mismatches.append((name.split("|", 1)[0], int(length),
                               int(transcript_len[index])))
        mapping[name] = index
    if mismatches:
        detail = "\n".join("    %s  @SQ %d  coordinate %d" % row for row in mismatches[:10])
        raise BuildError(
            "%d transcriptome reference(s) have an @SQ length differing from the shared "
            "coordinate. The two routes would not be on the same ruler.\n%s"
            % (len(mismatches), detail))
    return mapping

def read_txome_signals(bam_path, offsets, reference_map, coverage_offset,
                       transcript_len):
    """ONE pass over the transcriptome BAM -> both P-site indices and footprint ranges.

    The two signals pass exactly the same reads: uniquely-mapping primaries
    (`bam_inputs.is_unique_txome_read`, MAPQ >= 42), and a
    read length in the selected window (`offsets` supplies both the window and the offset).
    Reading the BAM twice to apply the same filter twice is pure I/O, and the outputs are
    compact int64 index arrays -- for a human library roughly 54 MB of P-site indices and
    107 MB of footprint ranges, against the 282 MB coverage array they feed.

    P-sites are placed at `reference_start + offset`, the forward transcript coordinate. The
    published transcriptome alignments are Bowtie2 `--norc`, so every read is on the
    transcript's forward strand. The CIGAR-aware walk is the GENOME route's rule, where a
    junction-spanning read would otherwise land in an intron; a transcript reference has no
    introns.

    Returns (psite_indices, footprint_starts, footprint_ends), all absolute indices into the
    concatenated coverage array.
    """
    import pysam

    psites, starts, ends = [], [], []
    bam = pysam.AlignmentFile(str(bam_path), "rb")
    try:
        for read in bam.fetch(until_eof=True):
            if not bam_inputs.is_unique_txome_read(read):
                continue
            offset = offsets.get(read.query_length)
            if offset is None:
                continue
            index = reference_map.get(read.reference_name)
            if index is None:
                continue
            base = coverage_offset[index]
            length = int(transcript_len[index])

            position = read.reference_start + offset
            if 0 <= position < length:
                psites.append(base + position)

            begin = max(0, read.reference_start)
            finish = min(length, read.reference_end)
            if finish > begin:
                starts.append(base + begin)
                ends.append(base + finish)
    finally:
        bam.close()
    return (np.asarray(psites, dtype=np.int64),
            np.asarray(starts, dtype=np.int64),
            np.asarray(ends, dtype=np.int64))

INT32_MAX = int(np.iinfo(np.int32).max)

def accumulate_points(indices, n_positions):
    """Point counts -> int32[n_positions].

    Accumulates directly into int32 rather than via `np.bincount`, whose int64 result would
    double the peak footprint on a 70.5-million-position coordinate. Overflow is impossible
    to reach here (the bound is the read count) but is checked anyway.
    """
    if indices.size and int(indices.size) > INT32_MAX:
        raise BuildError("more points (%d) than int32 can count" % indices.size)
    counts = np.zeros(n_positions, dtype=np.int32)
    if indices.size:
        np.add.at(counts, indices, 1)
    return counts

def accumulate_intervals(starts, ends, n_positions):
    """Half-open intervals -> depth, via a difference array.

    ONE global difference array is correct because every interval lies inside a single
    transcript's span, so the running sum returns to zero before the next transcript
    begins. That is asserted, not assumed: if any interval end escaped its transcript the
    array would not balance and depth would leak forward.

    Overflow is bounded BEFORE the cumulative sum, using the total positive increment as an
    upper bound on any achievable depth. That keeps the sum itself in int32 -- an int64
    cumsum would need 564 MB where 282 MB does.
    """
    if starts.size:
        if starts.size != ends.size:
            raise BuildError("got %d interval starts but %d ends"
                             % (starts.size, ends.size))
        if int(starts.min()) < 0 or int(ends.max()) > n_positions:
            raise BuildError(
                "an interval falls outside the coordinate [0, %d]: min start %d, max end "
                "%d. An end past the coordinate would leak depth into another transcript."
                % (n_positions, int(starts.min()), int(ends.max())))
        if int((ends < starts).sum()):
            raise BuildError("%d interval(s) have end < start"
                             % int((ends < starts).sum()))

    diff = np.zeros(n_positions + 1, dtype=np.int32)
    if starts.size:
        np.add.at(diff, starts, 1)
        np.add.at(diff, ends, -1)

    total = int(diff.sum(dtype=np.int64))
    if total != 0:
        raise BuildError(
            "the difference array does not balance (sum %d); some interval end is outside "
            "its transcript, which would leak depth into the next transcript" % total)
    upper_bound = int(diff[diff > 0].sum(dtype=np.int64)) if diff.size else 0
    if upper_bound > INT32_MAX:
        raise BuildError(
            "footprint depth could reach %d, which overflows int32 (max %d). Coverage is "
            "stored as int32; this needs a schema change, not a silent wrap."
            % (upper_bound, INT32_MAX))

    depth = np.cumsum(diff[:-1], dtype=np.int32)
    if depth.size and int(depth.min()) < 0:
        raise BuildError("footprint depth went negative -- intervals are malformed")
    return depth

#: The obvious vectorised forms allocate a full-length int64 array to produce ~19,736

def per_transcript_sums(values, coverage_offset, transcript_len):
    """Sum a full-coordinate array within each transcript's span."""
    return np.fromiter(
        (values[offset:offset + length].sum(dtype=np.int64)
         for offset, length in zip(coverage_offset, transcript_len)),
        dtype=np.int64, count=len(coverage_offset))

def region_slice_sums(values, coverage_offset, starts, ends):
    """Sum a full-coordinate array over one [start, end) window per transcript.

    `starts`/`ends` are transcript-relative; a window with end <= start (a CDS shorter than
    twice the trim, for instance) contributes 0.
    """
    return np.fromiter(
        (values[offset + max(int(start), 0):offset + max(int(end), int(start))]
         .sum(dtype=np.int64)
         for offset, start, end in zip(coverage_offset, starts, ends)),
        dtype=np.int64, count=len(coverage_offset))

def sha256_file(path):
    import hashlib
    digest = hashlib.sha256()
    with open(path, "rb") as handle:
        for chunk in iter(lambda: handle.read(1 << 20), b""):
            digest.update(chunk)
    return digest.hexdigest()

def file_identity(path, with_digest=True, record_path=False):
    """Identify an input by name, size and content digest.

    The full filesystem path is recorded ONLY under `--record-input-paths`. A coverage
    file is meant to be shareable, and an absolute path names the machine that built it,
    the directory layout of a private project, and often the operator. The `sha256` is
    what actually identifies the input; the path is convenience.
    """
    path = Path(path)
    record = {"name": path.name, "bytes": path.stat().st_size}
    if record_path:
        record["path"] = str(path)
    if with_digest:
        record["sha256"] = sha256_file(path)
    return record

def bam_identity(path, hash_bams=False, record_path=False):
    """Size + index digest by default. The index is derived from the BAM's block structure
    and alignment positions, so it cannot survive a content change; hashing many GB on the
    fast path would cost more than the check is worth. `--hash-bams` adds the full digest."""
    path = Path(path)
    record = {"name": path.name, "bytes": path.stat().st_size}
    if record_path:
        record["path"] = str(path)
    for suffix in (".bai", ".csi"):
        index = Path(str(path) + suffix)
        if index.exists():
            record["index"] = index.name
            record["index_sha256"] = sha256_file(index)
            break
    if hash_bams:
        record["sha256"] = sha256_file(path)
    return record

def build(config):
    """Build one sample's coverage file. Returns the finished path and a run report."""
    import coverage_schema
    import psite_placement
    import transcript_coords
    import transcript_regions

    started = time.time()
    report = {"sample": config.sample, "steps": {}}

    import annotation_cache

    bundle, reused = annotation_cache.load_or_build(
        getattr(config, "annotation_cache", None), config.gtf, config.appris,
        config.regions, config.left_span, config.right_span)
    report["annotation_cache_reused"] = reused

    headers = bundle["headers"]
    coords = bundle["coords"]
    cds_table = bundle["cds_table"]
    transcripts, exons = bundle["transcripts"], bundle["exons"]
    n_positions = bundle["n_positions"]
    region_summary = bundle["region_summary"]
    stop_ids = bundle["stop_ids"]
    index_of_id, index_of_base = bundle["index_of_id"], bundle["index_of_base"]
    regions, ribo_bins = bundle["regions"], bundle["ribo_bins"]

    transcripts = transcripts.copy()
    transcripts["region_offset"] = _first_row_offset(regions["transcript_index"].to_numpy(),
                                                     len(transcripts))
    transcripts["bin_offset"] = _first_row_offset(
        ribo_bins["transcript_index"].to_numpy() if len(ribo_bins)
        else np.empty(0, dtype=np.int64), len(transcripts))
    transcripts["in_transcriptome_reference"] = True
    transcripts["has_annotated_stop"] = [tid in stop_ids
                                         for tid in transcripts["transcript_id"]]
    transcripts["cds_divisible_by_3"] = (transcripts["cds_len_gtf"] % 3 == 0)
    transcripts["length_filtered"] = _length_filtered(headers, transcripts, regions)

    coverage_offset = transcripts["coverage_offset"].to_numpy()
    transcript_len = transcripts["transcript_len"].to_numpy()
    cds_total_by_id = dict(zip(transcripts["transcript_id"], transcripts["cds_len_gtf"]))

    log("offsets: reading the two QC masters")
    genome_offsets = psite_placement.load_offsets(config.qc_genome, config.sample)
    txome_offsets = psite_placement.load_offsets(config.qc_txome, config.sample)
    log("  genome %s" % genome_offsets)
    log("  txome  %s" % txome_offsets)

    exon_pr = _exon_pyranges(exons, transcripts)
    cds_pr = _cds_pyranges(cds_table)
    reference_map = txome_reference_map(config.txome_bam, index_of_base, transcript_len)
    log("  transcriptome references matched to the coordinate: %d" % len(reference_map))

    provenance = _provenance(config, coords, cds_table, region_summary,
                             genome_offsets, txome_offsets)
    bin_attrs = transcript_regions.ribo_bin_provenance(
        config.left_span, config.right_span, config.span_source)

    out_path = Path(config.output) / ("%s.shared_coverage.h5" % config.sample)
    writer = coverage_schema.CoverageWriter(
        out_path, sample=config.sample, transcripts=transcripts, exons=exons,
        regions=regions, ribo_bins=ribo_bins,
        offsets={"genome": genome_offsets, "transcriptome": txome_offsets},
        provenance=provenance, paper_cds_trim=config.trim,
        reference_name=config.reference_name, chunk=config.chunk,
        gzip_level=config.gzip_level, shuffle=config.shuffle,
        ribo_bin_attrs=bin_attrs, assay=getattr(config, "assay", "ribo"))

    try:
        cds_starts, cds_ends = _cds_windows(regions, len(transcripts))

        log("genome P-sites: streaming the BAM")
        chroms, positions, strands = read_genome_psites(
            config.genome_bam, genome_offsets)
        log("  %d reads placed; projecting" % len(chroms))
        indices, stats = project_genome_psites(
            chroms, positions, strands, exon_pr, cds_pr, cds_total_by_id,
            index_of_id, coverage_offset)
        del chroms, positions, strands
        report["steps"]["genome_psite"] = stats
        values = accumulate_points(indices, n_positions)
        del indices
        writer.write_signal("genome_psite", values)
        writer.set_transcript_counts(
            "n_genome_psite_events", per_transcript_sums(values, coverage_offset, transcript_len))
        writer.set_transcript_counts(
            "hist_cds_genome_psite_key",
            region_slice_sums(values, coverage_offset, cds_starts, cds_ends) > 0)
        del values

        log("transcriptome: streaming the BAM once for both signals")
        indices, txome_fp_starts, txome_fp_ends = read_txome_signals(
            config.txome_bam, txome_offsets,
            reference_map, coverage_offset, transcript_len)
        report["steps"]["txome_psite"] = {"n_placed": int(indices.size)}
        values = accumulate_points(indices, n_positions)
        del indices
        writer.write_signal("txome_psite", values)
        writer.set_transcript_counts(
            "n_txome_psite_events", per_transcript_sums(values, coverage_offset, transcript_len))
        writer.set_transcript_counts(
            "hist_cds_txome_psite_key",
            region_slice_sums(values, coverage_offset, cds_starts, cds_ends) > 0)
        del values

        trimmed_starts = cds_starts + config.trim
        trimmed_ends = cds_ends - config.trim

        log("genome footprints: streaming the BAM")
        blocks = read_genome_blocks(config.genome_bam, set(genome_offsets))
        log("  %d reads, %d aligned blocks; projecting" % (blocks[5], len(blocks[0])))
        starts, ends, stats = project_genome_footprints(
            blocks, exon_pr, cds_pr, index_of_id, coverage_offset)
        del blocks
        report["steps"]["genome_footprint"] = stats
        values = accumulate_intervals(starts, ends, n_positions)
        del starts, ends
        writer.write_signal("genome_footprint", values)
        writer.set_transcript_counts(
            "n_genome_footprint_bases",
            per_transcript_sums(values, coverage_offset, transcript_len))
        writer.set_transcript_counts(
            "hist_cds_genome_footprint_key",
            region_slice_sums(values, coverage_offset, trimmed_starts, trimmed_ends) > 0)
        del values

        # already read, in the same pass as the transcriptome P-sites
        starts, ends = txome_fp_starts, txome_fp_ends
        del txome_fp_starts, txome_fp_ends
        report["steps"]["txome_footprint"] = {"n_intervals": int(starts.size)}
        values = accumulate_intervals(starts, ends, n_positions)
        del starts, ends
        writer.write_signal("txome_footprint", values)
        writer.set_transcript_counts(
            "n_txome_footprint_bases",
            per_transcript_sums(values, coverage_offset, transcript_len))
        writer.set_transcript_counts(
            "hist_cds_txome_footprint_key",
            region_slice_sums(values, coverage_offset, trimmed_starts, trimmed_ends) > 0)
        del values

        final = writer.finalize()
    except BaseException:
        writer.abort()
        raise

    report["output"] = str(final)
    report["bytes"] = final.stat().st_size
    report["elapsed_seconds"] = round(time.time() - started, 1)
    report["n_transcripts"] = int(len(transcripts))
    report["n_positions"] = int(n_positions)
    log("wrote %s (%.1f MB) in %.1f min"
        % (final, report["bytes"] / 1e6, report["elapsed_seconds"] / 60.0))
    return final, report

def _first_row_offset(transcript_index, n_transcripts):
    """For each transcript, the index of its first row in a table sorted by transcript."""
    offsets = np.zeros(n_transcripts, dtype=np.int64)
    if transcript_index.size:
        counts = np.bincount(transcript_index, minlength=n_transcripts)
        offsets[1:] = np.cumsum(counts)[:-1]
    return offsets

def _cds_windows(regions, n_transcripts):
    """Per-transcript normalized CDS [start, end). Absent -> (0, 0)."""
    starts = np.zeros(n_transcripts, dtype=np.int64)
    ends = np.zeros(n_transcripts, dtype=np.int64)
    cds = regions[regions["label"] == "CDS"]
    starts[cds["transcript_index"].to_numpy()] = cds["start"].to_numpy()
    ends[cds["transcript_index"].to_numpy()] = cds["end"].to_numpy()
    return starts, ends

MIN_FIVE_UTR, MIN_CDS, MIN_THREE_UTR = 30, 150, 30

def _length_filtered(headers, transcripts, regions):
    """The 30/150/30 UTR/CDS minima, computed exactly as `ribo_seq_qc/config.py` does.

    Two details make this NOT the same as applying the thresholds to the normalized
    regions, and getting either one wrong silently redefines a documented flag:

      1. The lengths come from the RAW HEADER, so the CDS still INCLUDES the stop codon
         and the 3'UTR does not.
      2. A missing header UTR5 or UTR3 field filters the transcript outright. After
         normalization 495 transcripts gain a 3'UTR consisting solely of the relocated
         stop codon; those must still count as HAVING NO 3'UTR here.

    Carried as a flag, never a filter: nothing drops a transcript for it. `False` holds
    for 17,071 of the 19,736, which is the number to check a changed definition against.

    It asks whether a transcript's UTRs are long enough for a metagene window, which is a
    QC convenience and says nothing about whether its CDS can be counted -- so it is NOT a
    universe. Figure 4 counts all 19,736.
    """
    del transcripts, regions                     # deliberately unused: the header decides
    flags = []
    for tid in sorted(headers):
        header = headers[tid]
        utr5 = header["hdr_utr5"]
        utr3 = header["hdr_utr3"]
        cds_start, cds_end = header["hdr_cds"]
        cds_len = cds_end - cds_start + 1
        flags.append(
            (utr5 is None or (utr5[1] - utr5[0] + 1) < MIN_FIVE_UTR)
            or cds_len < MIN_CDS
            or (utr3 is None or (utr3[1] - utr3[0] + 1) < MIN_THREE_UTR))
    return np.asarray(flags, dtype=bool)

def _provenance(config, coords, cds_table, region_summary, genome_offsets, txome_offsets):
    import coverage_schema
    import h5py
    import psite_placement
    import pysam
    import scipy

    keep_paths = bool(getattr(config, "record_input_paths", False))
    inputs = {
        "gtf": file_identity(config.gtf, record_path=keep_paths),
        "appris_lengths": file_identity(config.appris, record_path=keep_paths),
        "genome_bam": bam_identity(config.genome_bam, config.hash_bams, keep_paths),
        "transcriptome_bam": bam_identity(config.txome_bam, config.hash_bams, keep_paths),
        "qc_genome": file_identity(config.qc_genome, record_path=keep_paths),
        "qc_transcriptome": file_identity(config.qc_txome, record_path=keep_paths),
    }
    if config.regions:
        inputs["actual_regions_bed"] = file_identity(config.regions, record_path=keep_paths)
    return {
        "schema": coverage_schema.SCHEMA,
        "sample": config.sample,
        "assay": getattr(config, "assay", "ribo"),
        "routes": list(coverage_schema.ROUTES),
        "generation": coverage_schema.invocation(record_paths=keep_paths),
        "code_version": coverage_schema.code_version(),
        "inputs": inputs,
        "parameters": {
            "paper_cds_trim": config.trim,
            "genome_uniqueness": "NH==1",
            "txome_uniqueness": "MAPQ>=%d" % bam_inputs.txome_min_mapq(),
            "left_span": config.left_span,
            "right_span": config.right_span,
            "psite_placement": psite_placement.PSITE_PLACEMENT,
            "appris_principal_ranks_consumed": False,
        },
        "assignment_policies": {
            "psite": {"rule": "first_exon_overlap", "stage1_exon_set": "cds_only",
                      "tie_break": "historical_join_order",
                      "utr_fallback": "first_exon_overlap"},
            "footprint": {"rule": "max_exon_overlap", "stage1_exon_set": "cds_only",
                          "tie_break": "historical_join_order",
                          "utr_fallback": "max_exon_overlap"},
        },
        "regions": region_summary,
        "counts": {"n_transcripts": int(len(coords["transcripts"])),
                   "n_exons": int(len(coords["exons"])),
                   "n_cds_exons": int(len(cds_table)),
                   "n_positions": int(coords["n_positions"])},
        "offsets": {"genome": {str(k): v for k, v in sorted(genome_offsets.items())},
                    "transcriptome": {str(k): v for k, v in sorted(txome_offsets.items())}},
        "software": {"numpy": np.__version__, "pandas": pd.__version__,
                     "pysam": pysam.__version__, "scipy": scipy.__version__,
                     "h5py": h5py.__version__,
                     "python": "%d.%d.%d" % sys.version_info[:3]},
    }

def _build_parser():
    parser = argparse.ArgumentParser(
        description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument("--sample", required=True)
    parser.add_argument("--genome-bam", required=True, type=Path)
    parser.add_argument("--transcriptome-bam", required=True, type=Path, dest="txome_bam")
    parser.add_argument("--gtf", required=True, type=Path)
    parser.add_argument("--appris", required=True, type=Path)
    parser.add_argument("--regions", type=Path,
                        help="appris_human_v2_actual_regions.bed -- a cross-check")
    parser.add_argument("--annotation-cache", type=Path, default=None,
                        help="reuse (or create) the sample-independent annotation bundle "
                             "here instead of reparsing the GTF. The cohort driver builds "
                             "it once and passes it to every sample.")
    parser.add_argument("--qc-genome", required=True, type=Path)
    parser.add_argument("--qc-txome", required=True, type=Path)
    parser.add_argument("--output", type=Path, default=Path("results/coverage"))
    parser.add_argument("--trim", type=int, default=15)
    parser.add_argument("--left-span", type=int, default=35)
    parser.add_argument("--right-span", type=int, default=10)
    parser.add_argument("--assay", default="ribo", choices=("ribo", "rna"),
                        help="recorded in the file's provenance; both BAMs must be the "
                             "same assay, since the two routes are compared to each other")
    parser.add_argument("--reference-name", default="appris_human_v2_selected")
    parser.add_argument("--chunk", type=int, default=1 << 16)
    parser.add_argument("--gzip-level", type=int, default=9)
    parser.add_argument("--no-shuffle", dest="shuffle", action="store_false", default=True)
    parser.add_argument("--hash-bams", action="store_true")
    parser.add_argument("--record-input-paths", action="store_true",
                        help="store full filesystem paths in the provenance. Off by "
                             "default: an absolute path names the machine that built the "
                             "file, and the sha256 is what identifies the input.")
    return parser

REQUIRED_INPUTS = (("--genome-bam", "genome_bam"), ("--transcriptome-bam", "txome_bam"),
                   ("--gtf", "gtf"), ("--appris", "appris"),
                   ("--qc-genome", "qc_genome"), ("--qc-txome", "qc_txome"))

def check_inputs(args):
    """Fail before any compute starts, naming every missing input at once."""
    missing = ["  %-22s %s" % (flag, getattr(args, attr))
               for flag, attr in REQUIRED_INPUTS if not Path(getattr(args, attr)).exists()]
    if args.regions and not Path(args.regions).exists():
        missing.append("  %-22s %s" % ("--regions", args.regions))
    if missing:
        raise SystemExit("these required inputs do not exist:\n" + "\n".join(missing))

def main(argv=None):
    args = _build_parser().parse_args(argv)
    args.span_source = ("cli" if (args.left_span, args.right_span) != (35, 10)
                        else "cli_default")
    if args.trim % 3:
        raise SystemExit("--trim must be a multiple of 3 to keep the CDS slice in frame; "
                         "got %d" % args.trim)
    check_inputs(args)
    _final, report = build(args)
    return 0

if __name__ == "__main__":
    sys.path.insert(0, str(Path(__file__).resolve().parent))
    sys.exit(main())
