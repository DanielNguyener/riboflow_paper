#!/usr/bin/env python3
"""Build one sample's shared-coordinate coverage HDF5 from its two ribo BAMs."""
from __future__ import annotations

import argparse
import re
import sys
import time
from pathlib import Path

import numpy as np
import pandas as pd

sys.path.insert(0, str(Path(__file__).resolve().parent.parent / "common"))
import bam_inputs                      # the one uniqueness policy lives here
from inputs import make_log, require_existing, sha256_of

_CDS_HEADER = re.compile(r"\|CDS:(\d+)-(\d+)\|")
REFERENCE_NAME = "appris_human_v2_selected"

#: junction window spans that key the annotation cache  schema 3 stores no junction
#: bins so these only fingerprint the cache  no stored number changes
LEFT_SPAN, RIGHT_SPAN = 35, 10


class BuildError(RuntimeError):
    pass

log = make_log("coverage")

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

def read_genome_signals(bam_path, offsets, on_record=None, counts_sink=None,
                        collect_coverage=True):
    """Stream the genome BAM once for everything the ribo pass needs from it.

    Coverage keeps unique primaries with a read length in `offsets`. P-sites are
    CIGAR-aware (never inside an intron) and undefined placements are dropped from the
    P-site signal only; footprints are the aligned blocks (`get_blocks()` splits on N,
    introns excluded).

    `on_record` (categories) is called first for EVERY fetched record -- it needs
    secondaries, so no filter runs before it. `counts_sink` (Figure 3), three lists,
    receives (chrom, raw 5' end, strand) for exactly the coverage-filtered reads: no
    offset is applied, and the key set of `offsets` IS the selected-length window.
    """
    import pysam
    if collect_coverage:
        import psite_placement

    p_chroms, p_positions, p_strands = [], [], []
    f_chroms, f_starts, f_ends, f_strands, f_read_ids = [], [], [], [], []
    index = 0
    collect = collect_coverage or counts_sink is not None
    bam = pysam.AlignmentFile(str(bam_path), "rb")
    try:
        for read in bam.fetch(until_eof=True):
            if on_record is not None:
                on_record(read)
            if not collect:
                continue
            if not bam_inputs.is_unique_genome_read(read):
                continue
            offset = offsets.get(read.query_length)
            if offset is None:
                continue
            strand = "-" if read.is_reverse else "+"

            if counts_sink is not None:
                counts_sink[0].append(read.reference_name)
                counts_sink[1].append(read.reference_end - 1 if read.is_reverse
                                      else read.reference_start)
                counts_sink[2].append(strand)
            if not collect_coverage:
                continue

            position = psite_placement.place(read, offset)
            if position is not None:
                p_chroms.append(read.reference_name)
                p_positions.append(position)
                p_strands.append(strand)

            for block_start, block_end in read.get_blocks():
                f_chroms.append(read.reference_name)
                f_starts.append(block_start)
                f_ends.append(block_end)
                f_strands.append(strand)
                f_read_ids.append(index)
            index += 1
    finally:
        bam.close()
    psites = (p_chroms, np.asarray(p_positions, dtype=np.int64), p_strands)
    blocks = (f_chroms, np.asarray(f_starts, dtype=np.int64),
              np.asarray(f_ends, dtype=np.int64), f_strands,
              np.asarray(f_read_ids, dtype=np.int64), index)
    return psites, blocks


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
    """Genome P-sites -> absolute indices into the concatenated coverage array; (indices, stats).

    Stage 1 (CDS exons only) must stay separate from stage 2 (full exon map): merging them
    would change what `keep="first"`/`idxmax` select and could move a published CDS value.
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

    Stage 1: max-CDS-overlap on CDS exons alone; stage 2: full exon map for unclaimed reads.
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
    # on "-" transcript runs the other way  overlap 5' end is its genomic END
    tx_start = np.where(plus, ex_tx_start + (ostart - ex_g_start),
                        ex_tx_start + (ex_g_end - oend))
    base = coverage_offset[rows["tx_index"].to_numpy()]
    return (base + tx_start, base + tx_start + (oend - ostart),
            {"n_assigned_stage1": n_stage1, "n_assigned_stage2_utr": n_stage2,
             "n_unassigned": int((assigned == -1).sum())})

def txome_reference_map(bam_path, tx_index_of_base, transcript_len):
    """{reference name: transcript index}, cross-checking @SQ length == transcript_len
    (the transcriptome reference IS the shared coordinate)."""
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
                       transcript_len, on_record=None, counts_tally=None,
                       collect_coverage=True):
    """One transcriptome-BAM pass -> (psite_indices, footprint_starts, footprint_ends).

    P-site = reference_start + offset (Bowtie2 --norc: all reads forward; no introns, so no
    CIGAR walk needed); both signals pass exactly the same reads (MAPQ >= 42, window length).

    `on_record` (categories) is called first for EVERY fetched record. `counts_tally`
    (Figure 3, a `ribo_rna_lib.TxomeCdsTally`) receives exactly the reads that pass the
    uniqueness and window filters -- BEFORE the reference-map filter, which only coverage
    applies.
    """
    import pysam

    psites, starts, ends = [], [], []
    collect = collect_coverage or counts_tally is not None
    bam = pysam.AlignmentFile(str(bam_path), "rb")
    try:
        for read in bam.fetch(until_eof=True):
            if on_record is not None:
                on_record(read)
            if not collect:
                continue
            if not bam_inputs.is_unique_txome_read(read):
                continue
            offset = offsets.get(read.query_length)
            if offset is None:
                continue
            if counts_tally is not None:
                counts_tally.add(read)
            if not collect_coverage:
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
    """Point counts -> int32[n_positions] directly (an int64 bincount would double peak memory)."""
    if indices.size and int(indices.size) > INT32_MAX:
        raise BuildError("more points (%d) than int32 can count" % indices.size)
    counts = np.zeros(n_positions, dtype=np.int32)
    if indices.size:
        np.add.at(counts, indices, 1)
    return counts

def accumulate_intervals(starts, ends, n_positions):
    """Half-open intervals -> depth via one global difference array.

    Correct only because every interval stays inside one transcript's span (asserted below);
    int32 overflow is bounded BEFORE the cumsum so the sum itself stays int32.
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

def file_identity(path, record_path=False):
    """Identify an input by name, size and content digest.

    The full path is recorded only under `--record-input-paths` (shareable files, no machine names).
    """
    path = Path(path)
    record = {"name": path.name, "bytes": path.stat().st_size}
    if record_path:
        record["path"] = str(path)
    record["sha256"] = sha256_of(path)
    return record

def bam_identity(path, hash_bams=False, record_path=False):
    """Size + index digest by default (the index cannot survive a content change);
    `--hash-bams` adds the full digest."""
    path = Path(path)
    record = {"name": path.name, "bytes": path.stat().st_size}
    if record_path:
        record["path"] = str(path)
    for suffix in (".bai", ".csi"):
        index = Path(str(path) + suffix)
        if index.exists():
            record["index"] = index.name
            record["index_sha256"] = sha256_of(index)
            break
    if hash_bams:
        record["sha256"] = sha256_of(path)
    return record

PRODUCTS = ("coverage", "counts", "categories")

def _requested_products(config):
    """The products this run builds; `--only coverage` is the historical behavior."""
    value = getattr(config, "only", None) or ("coverage",)
    if isinstance(value, str):
        value = tuple(part.strip() for part in value.split(",") if part.strip())
    unknown = sorted(set(value) - set(PRODUCTS))
    if unknown:
        raise BuildError("unknown product(s) %s; valid: %s"
                         % (", ".join(unknown), ", ".join(PRODUCTS)))
    return set(value)

def _ribo_rna_lib():
    sys.path.insert(0, str(Path(__file__).resolve().parent.parent / "ribo_rna"))
    import ribo_rna_lib
    return ribo_rna_lib

def _library_scan():
    sys.path.insert(0, str(Path(__file__).resolve().parent.parent / "read_categories"))
    import library_scan
    return library_scan

def build(config):
    """Build one sample's requested products from ONE pass over each ribo BAM.

    coverage (Fig 2, the HDF5), counts (Fig 3, the staged ribo count column pair) and
    categories (Fig 4, the three staged rows) share the two traversals; each product
    applies its own per-record filter. Returns the coverage path (None when coverage was
    not requested) and a run report.
    """
    import coverage_schema
    import psite_placement
    started = time.time()
    products = _requested_products(config)
    do_cov = "coverage" in products
    do_counts = "counts" in products
    do_cats = "categories" in products
    report = {"sample": config.sample, "steps": {}}

    bundle = None
    if do_cov or do_counts:
        import annotation_cache

        bundle, reused = annotation_cache.load_or_build(
            config.annotation_cache, config.gtf, config.appris,
            config.regions, LEFT_SPAN, RIGHT_SPAN)
        report["annotation_cache_reused"] = reused

        log("offsets: reading the two QC masters")
        genome_offsets = psite_placement.load_offsets(config.qc_genome, config.sample)
        txome_offsets = psite_placement.load_offsets(config.qc_txome, config.sample)
        log("  genome %s" % genome_offsets)
        log("  txome  %s" % txome_offsets)
    else:
        genome_offsets = txome_offsets = None

    reference_map = coverage_offset = transcript_len = writer = None
    if do_cov:
        headers = bundle["headers"]
        coords = bundle["coords"]
        cds_table = bundle["cds_table"]
        transcripts, exons = bundle["transcripts"], bundle["exons"]
        n_positions = bundle["n_positions"]
        region_summary = bundle["region_summary"]
        index_of_id, index_of_base = bundle["index_of_id"], bundle["index_of_base"]
        regions = bundle["regions"]

        transcripts = transcripts.copy()
        cds_starts, cds_ends = _cds_windows(regions, len(transcripts))
        transcripts["cds_start"] = cds_starts
        transcripts["cds_end"] = cds_ends
        del headers

        coverage_offset = transcripts["coverage_offset"].to_numpy()
        transcript_len = transcripts["transcript_len"].to_numpy()
        cds_total_by_id = dict(zip(transcripts["transcript_id"],
                                   transcripts["cds_len_gtf"]))

        exon_pr = _exon_pyranges(exons, transcripts)
        cds_pr = _cds_pyranges(cds_table)
        reference_map = txome_reference_map(config.txome_bam, index_of_base,
                                            transcript_len)
        log("  transcriptome references matched to the coordinate: %d"
            % len(reference_map))

        provenance = _provenance(config, coords, cds_table, region_summary,
                                 genome_offsets, txome_offsets)

        out_path = Path(config.output) / ("%s.shared_coverage.h5" % config.sample)
        writer = coverage_schema.CoverageWriter(
            out_path, sample=config.sample,
            transcripts=transcripts[list(coverage_schema.TRANSCRIPT_COLUMNS)],
            provenance=provenance, paper_cds_trim=config.trim,
            chunk=getattr(config, "chunk", 1 << 16),
            gzip_level=getattr(config, "gzip_level", 9),
            assay=getattr(config, "assay", "ribo"))

    counts_tally = counts_sink = universe = rrl = None
    if do_counts:
        import pysam

        counts_staging = getattr(config, "counts_staging", None)
        if not counts_staging:
            raise BuildError("counts need a staging directory (--counts-staging)")
        rrl = _ribo_rna_lib()
        for path in (config.genome_bam, config.txome_bam):
            rrl.assert_single_end(path)
        universe, universe_report = rrl.build_universe(bundle, config.txome_bam)
        log("counts: shared reference set %d APPRIS transcripts"
            % universe_report["n_universe"])
        cds_spans = rrl.transcript_cds_spans(bundle, universe)
        counts_cds_pr = rrl.genome_cds_intervals(bundle, universe)
        handle = pysam.AlignmentFile(str(config.txome_bam), "rb")
        try:
            counts_tally = rrl.TxomeCdsTally(handle.references, cds_spans)
        finally:
            handle.close()
        counts_sink = ([], [], [])

    cat_state = ls = None
    if do_cats:
        ls = _library_scan()
        cat_state = ls.new_state()

    # transcriptome BAM goes FIRST  categories genome collector tests membership in the
    # transcriptome presence set  set must be complete by then  h5 WRITE order below
    # is unchanged
    try:
        log("transcriptome: streaming the BAM once")
        tx_psites, txome_fp_starts, txome_fp_ends = read_txome_signals(
            config.txome_bam, txome_offsets, reference_map, coverage_offset,
            transcript_len,
            on_record=(lambda r: ls.collect_txome_record(cat_state, r))
            if do_cats else None,
            counts_tally=counts_tally, collect_coverage=do_cov)

        log("genome: streaming the BAM once")
        (chroms, positions, strands), genome_blocks = read_genome_signals(
            config.genome_bam, genome_offsets,
            on_record=(lambda r: ls.collect_genome_record(cat_state, r))
            if do_cats else None,
            counts_sink=counts_sink, collect_coverage=do_cov)

        if do_cov:
            log("  %d reads placed; projecting" % len(chroms))
            indices, stats = project_genome_psites(
                chroms, positions, strands, exon_pr, cds_pr, cds_total_by_id,
                index_of_id, coverage_offset)
            del chroms, positions, strands
            report["steps"]["genome_psite"] = stats
            values = accumulate_points(indices, n_positions)
            del indices
            writer.write_signal("genome_psite", values)
            del values

            report["steps"]["txome_psite"] = {"n_placed": int(tx_psites.size)}
            values = accumulate_points(tx_psites, n_positions)
            del tx_psites
            writer.write_signal("txome_psite", values)
            del values

            # already read  same pass as the genome P sites
            blocks = genome_blocks
            del genome_blocks
            log("  %d reads, %d aligned blocks; projecting" % (blocks[5], len(blocks[0])))
            starts, ends, stats = project_genome_footprints(
                blocks, exon_pr, cds_pr, index_of_id, coverage_offset)
            del blocks
            report["steps"]["genome_footprint"] = stats
            values = accumulate_intervals(starts, ends, n_positions)
            del starts, ends
            writer.write_signal("genome_footprint", values)
            del values

            # already read  same pass as the transcriptome P sites
            starts, ends = txome_fp_starts, txome_fp_ends
            del txome_fp_starts, txome_fp_ends
            report["steps"]["txome_footprint"] = {"n_intervals": int(starts.size)}
            values = accumulate_intervals(starts, ends, n_positions)
            del starts, ends
            writer.write_signal("txome_footprint", values)
            del values

            final = writer.finalize()
        else:
            final = None
    except BaseException:
        if writer is not None:
            writer.abort()
        raise

    if do_counts:
        counts, n_assigned, n_ambiguous = rrl.count_genome_points(
            counts_sink[0], counts_sink[1], counts_sink[2], counts_cds_pr,
            stranded=True)
        frame = pd.DataFrame({"transcript_id": universe})
        frame["genome_ribo_reads"] = np.array(
            [counts.get(t, 0) for t in universe], dtype=int)
        frame["txome_ribo_reads"] = np.array(
            [counts_tally.counts.get(t, 0) for t in universe], dtype=int)
        del counts_sink, counts
        staging = Path(config.counts_staging)
        staging.mkdir(parents=True, exist_ok=True)
        counts_path = staging / ("%s.tsv" % config.sample)
        frame.to_csv(counts_path, sep="\t", index=False, lineterminator="\n")
        report["steps"]["counts"] = {
            "genome_assigned": int(n_assigned), "genome_ambiguous": int(n_ambiguous),
            "txome_assigned": int(counts_tally.n_assigned)}
        log("counts: staged %s (genome %d, txome %d CDS-assigned reads)"
            % (counts_path, n_assigned, counts_tally.n_assigned))

    if do_cats:
        # row builders load annotation through the RIBOFLOW_PAPER_* environment  same as
        # the standalone scan always did  fill it in for standalone runs
        import os
        os.environ.setdefault("RIBOFLOW_PAPER_GTF", str(config.gtf))
        os.environ.setdefault("RIBOFLOW_PAPER_APPRIS", str(config.appris))
        ls.finish_state(cat_state)
        ls.stage_rows(config.sample, cat_state, log=log)

    report["elapsed_seconds"] = round(time.time() - started, 1)
    if do_cov:
        report["output"] = str(final)
        report["bytes"] = final.stat().st_size
        report["n_transcripts"] = int(len(transcripts))
        report["n_positions"] = int(n_positions)
        log("wrote %s (%.1f MB) in %.1f min"
            % (final, report["bytes"] / 1e6, report["elapsed_seconds"] / 60.0))
    return final, report

def _cds_windows(regions, n_transcripts):
    """Per-transcript normalized CDS [start, end), stop codon excluded. Absent -> (-1, -1)."""
    import coverage_schema
    starts = np.full(n_transcripts, coverage_schema.NO_CDS, dtype=np.int64)
    ends = np.full(n_transcripts, coverage_schema.NO_CDS, dtype=np.int64)
    cds = regions[regions["label"] == "CDS"]
    starts[cds["transcript_index"].to_numpy()] = cds["start"].to_numpy()
    ends[cds["transcript_index"].to_numpy()] = cds["end"].to_numpy()
    return starts, ends

def _provenance(config, coords, cds_table, region_summary, genome_offsets, txome_offsets):
    import coverage_schema
    import h5py
    import psite_placement
    import pysam
    import scipy

    keep_paths = bool(getattr(config, "record_input_paths", False))
    hash_bams = bool(getattr(config, "hash_bams", False))
    inputs = {
        "gtf": file_identity(config.gtf, record_path=keep_paths),
        "appris_lengths": file_identity(config.appris, record_path=keep_paths),
        "genome_bam": bam_identity(config.genome_bam, hash_bams, keep_paths),
        "transcriptome_bam": bam_identity(config.txome_bam, hash_bams, keep_paths),
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
            "psite_placement": psite_placement.PSITE_PLACEMENT,
            "stop_codon_assignment": "utr3",
            "exon_source": "gencode_exon_features",
            "reference_name": REFERENCE_NAME,
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
    parser.add_argument("--regions", type=Path)
    parser.add_argument("--annotation-cache", type=Path, default=None)
    parser.add_argument("--qc-genome", type=Path)
    parser.add_argument("--qc-txome", type=Path)
    parser.add_argument("--only", default="coverage")
    parser.add_argument("--counts-staging", type=Path, default=None, dest="counts_staging")
    parser.add_argument("--output", type=Path, default=Path("results/coverage"))
    parser.add_argument("--trim", type=int, default=15)
    return parser

REQUIRED_INPUTS = (("--genome-bam", "genome_bam"), ("--transcriptome-bam", "txome_bam"),
                   ("--gtf", "gtf"), ("--appris", "appris"))
QC_INPUTS = (("--qc-genome", "qc_genome"), ("--qc-txome", "qc_txome"))

def check_inputs(args):
    """Fail before any compute starts, naming every missing input at once.

    The QC masters are inputs to coverage and counts only; a categories-only run never
    reads them.
    """
    products = _requested_products(args)
    required = list(REQUIRED_INPUTS)
    if products & {"coverage", "counts"}:
        required += list(QC_INPUTS)
    pairs = [(flag, getattr(args, attr)) for flag, attr in required]
    if args.regions:
        pairs.append(("--regions", args.regions))
    require_existing(pairs)
    if "counts" in products and not getattr(args, "counts_staging", None):
        raise SystemExit("--counts-staging is required with --only ...counts...")

def main(argv=None):
    args = _build_parser().parse_args(argv)
    if args.trim % 3:
        raise SystemExit("--trim must be a multiple of 3 to keep the CDS slice in frame; "
                         "got %d" % args.trim)
    check_inputs(args)
    build(args)
    return 0

if __name__ == "__main__":
    sys.path.insert(0, str(Path(__file__).resolve().parent))
    sys.exit(main())
