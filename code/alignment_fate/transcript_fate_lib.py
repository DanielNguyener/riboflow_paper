#!/usr/bin/env python3
"""Per-transcript alignment fates: where a transcript's transcriptome reads went in the genome."""
from __future__ import annotations

import sys
from pathlib import Path

import numpy as np
import pandas as pd

REPO = Path(__file__).resolve().parents[2]
READ_TAXONOMY = REPO / "code" / "read_taxonomy"

PROCESSED_PSEUDOGENE = "processed_pseudogene"
CATEGORIES = ("genome_unique", "genome_multi_pseudogene_tie", "other")
CATEGORY_LABELS = {
    "genome_unique": "genome-unique (kept)",
    "genome_multi_pseudogene_tie": "genome-multi: pseudogene tie (excluded)",
    "other": "other (genome-multi w/o tie, or unaligned)",
}
WIDE_COLUMNS = [
    "sample", "tid", "n_txome_assigned", "n_txome_all_primary", "n_shared_both_bams",
    "n_genome_unique", "n_genome_multi", "n_top_on_target",
    "n_top_on_target_and_pp_tie", "n_top_on_target_and_any_pseudo_tie",
    "pct_pp_tie_of_assigned", "pct_pp_tie_of_shared", "pct_pp_tie_of_genome_multi"]
TIDY_COLUMNS = ["sample", "transcript_id", "gene_id", "gene_name", "category",
                "n_reads", "pct_of_txome_assigned"]

PANEL_TRANSCRIPTS = (("ENST00000396861.5", "GAPDH"),
                     ("ENST00000361682.11", "COMT"))

class FateError(RuntimeError):
    pass

def load_libraries():
    """Import the read-taxonomy libraries by path, without disturbing sys.path."""
    saved = list(sys.path)
    try:
        for directory in (str(READ_TAXONOMY), str(REPO / "code" / "common"),
                          str(REPO / "code" / "common" / "ribo_seq_qc")):
            if directory not in sys.path:
                sys.path.insert(0, directory)
        import concordance_lib
        import mm_concordance_lib
        return concordance_lib, mm_concordance_lib
    finally:
        sys.path[:] = saved + [p for p in sys.path if p not in saved]

def resolve_transcripts(table, gene_ids=(), transcript_ids=(), coverage=None):
    """Gene and/or transcript IDs -> versioned transcript IDs present in `table`.

    Refuses to guess when a gene maps to more than one candidate, listing them. When a
    coverage HDF5 is supplied it is used for resolution (it carries gene IDs directly);
    otherwise the APPRIS transcript table is searched.
    """
    resolved = []
    for tid in transcript_ids:
        if tid in table:
            resolved.append(tid)
            continue
        base = tid.split(".", 1)[0]
        hits = [t for t in table if t.split(".", 1)[0] == base]
        if len(hits) == 1:
            resolved.append(hits[0])
        elif not hits:
            raise FateError("transcript %r is not in the APPRIS transcript table" % tid)
        else:
            raise FateError("transcript %r is ambiguous: %s" % (tid, ", ".join(hits)))

    for gene_id in gene_ids:
        if coverage is not None:
            index = coverage.resolve_gene(gene_id)
            resolved.append(coverage.transcript_info(index)["transcript_id"])
            continue
        base = gene_id.split(".", 1)[0]
        hits = [tid for tid, entry in table.items()
                if str(entry.get("gene_id", "")).split(".", 1)[0] == base]
        if len(hits) == 1:
            resolved.append(hits[0])
        elif not hits:
            raise FateError("gene %r has no transcript in the APPRIS table" % gene_id)
        else:
            raise FateError(
                "gene %r maps to %d transcripts and no unique one resolves it. Pass "
                "--transcript-id to choose; this is not guessed.\n%s"
                % (gene_id, len(hits), "\n".join("    %s" % h for h in sorted(hits))))

    seen, ordered = set(), []
    for tid in resolved:
        if tid not in seen:
            seen.add(tid)
            ordered.append(tid)
    return ordered

def collect_txome_populations(concordance_lib, txome_bam, base2ver, tids):
    """One transcriptome pass -> per transcript, its confident-unique read population.

    Returns (populations, all_primary_counts) where populations[tid] maps qname ->
    (tid, tx_pos, tx_len) for reads whose PRIMARY resolves there at MAPQ >= 42.
    """
    present, unique_qnames, _all_qnames = concordance_lib.read_txome_primary(
        txome_bam, base2ver)
    wanted = set(tids)
    populations = {tid: {} for tid in tids}
    all_primary = {tid: 0 for tid in tids}
    for qname, value in present.items():
        tid = value[0]
        if tid not in wanted:
            continue
        all_primary[tid] += 1
        if qname in unique_qnames:
            populations[tid][qname] = value
    return populations, all_primary

def collect_genome_records(genome_bam, target_qnames):
    """One genome pass -> primary presence, uniqueness, and every reported locus.

    One pass, not two: presence, uniqueness and every reported locus all come from the
    same scan of the genome BAM. Uniqueness is the
    repository's one genome rule -- the PRIMARY record's `NH == 1`
    (`bam_inputs.is_unique_genome_read`). There is no MAPQ fallback: this panel's whole
    subject is multimapping, so a BAM that cannot report `NH` cannot answer its question
    and is rejected rather than approximated.
    """
    import pysam
    from collections import defaultdict

    sys.path.insert(0, str(Path(__file__).resolve().parent.parent / "common"))
    import bam_inputs

    present = set()
    unique = set()
    primary_multi = set()
    records = defaultdict(list)

    bam = pysam.AlignmentFile(str(genome_bam), "rb")
    try:
        for read in bam.fetch(until_eof=True):
            if read.is_unmapped or read.is_supplementary:
                continue
            qname = read.query_name
            if qname not in target_qnames:
                continue
            if not read.is_secondary:
                present.add(qname)
                if bam_inputs.is_unique_genome_read(read):
                    unique.add(qname)
            nh = int(read.get_tag("NH"))
            if not read.is_secondary and nh > 1:
                primary_multi.add(qname)
            if nh <= 1:
                continue
            blocks = read.get_blocks()
            if not blocks:
                continue
            try:
                score = read.get_tag("AS")
            except KeyError:
                score = None
            pos5 = read.reference_end - 1 if read.is_reverse else read.reference_start
            records[qname].append((
                read.reference_name, "-" if read.is_reverse else "+", pos5,
                len(blocks), min(b[0] for b in blocks), max(b[1] for b in blocks), score))
    finally:
        bam.close()
    return present, unique, {q: v for q, v in records.items() if q in primary_multi}

def _project_match(concordance_lib, frame, transcript):
    """Is each genome locus concordant with `transcript`? Verbatim Part-14 projection."""
    cum_start = transcript["cum_start"]
    g_start = transcript["g_start"]
    g_end = transcript["g_end"]
    n_exons = len(cum_start)
    tx_pos = frame["tx_pos"].to_numpy()
    tx_len = frame["tx_len"].to_numpy()
    idx_start = np.clip(np.searchsorted(cum_start, tx_pos, side="right") - 1, 0, n_exons - 1)
    idx_end = np.clip(np.searchsorted(cum_start, tx_pos + tx_len - 1, side="right") - 1,
                      0, n_exons - 1)
    expected_junctions = idx_end - idx_start
    offset = tx_pos - cum_start[idx_start]
    expected_pos5 = (g_start[idx_start] + offset if transcript["strand"] == "+"
                     else g_end[idx_start] - 1 - offset)

    same_chrom = frame["g_chrom"].to_numpy() == transcript["chrom"]
    within_body = (same_chrom
                   & (frame["g_min"].to_numpy() < transcript["body_end"])
                   & (frame["g_max"].to_numpy() > transcript["body_start"]))
    strand_ok = frame["g_strand"].to_numpy() == transcript["strand"]
    splice_ok = (frame["g_nblocks"].to_numpy() - 1) == expected_junctions
    coord_ok = np.abs(frame["g_pos5"].to_numpy(dtype=np.int64)
                      - expected_pos5.astype(np.int64)) <= concordance_lib.COORD_TOL
    return within_body & strand_ok & splice_ok & coord_ok

def classify_transcript(concordance_lib, tid, transcript, population, all_primary,
                        genome_present, genome_unique, records, exon_pr, sample):
    """The wide summary row for one transcript, plus its three-way partition."""
    import pyranges as pr

    n_txome_assigned = len(population)
    if n_txome_assigned == 0:
        raise FateError(
            "%s: no transcriptome-assigned reads for %s. Either the sample has no coverage "
            "there, or the transcript id is not in this reference." % (sample, tid))

    shared = {q for q in population if q in genome_present}
    n_shared = len(shared)
    n_unique = len(shared & genome_unique)
    transcript_records = {q: v for q, v in records.items() if q in population}
    n_multi = len(transcript_records)

    rows = []
    for qname, locus_list in transcript_records.items():
        _tid, tx_pos, tx_len = population[qname]
        for chrom, strand, pos5, n_blocks, g_min, g_max, score in locus_list:
            rows.append((qname, tx_pos, tx_len, chrom, strand, pos5, n_blocks,
                         g_min, g_max, score))

    n_target_at_max = n_pp_tie = n_any_pseudo_tie = 0
    if rows:
        frame = pd.DataFrame(rows, columns=[
            "qname", "tx_pos", "tx_len", "g_chrom", "g_strand", "g_pos5", "g_nblocks",
            "g_min", "g_max", "AS"])
        frame["match_target"] = _project_match(concordance_lib, frame, transcript)
        frame["max_as"] = frame.groupby("qname", sort=False)["AS"].transform("max")
        frame["at_max"] = frame["AS"] == frame["max_as"]

        locus = frame[["g_chrom", "g_pos5"]].reset_index().rename(columns={"index": "row"})
        locus_pr = pr.PyRanges(pd.DataFrame({
            "Chromosome": locus["g_chrom"], "Start": locus["g_pos5"],
            "End": locus["g_pos5"] + 1, "row": locus["row"]}))
        joined = locus_pr.join(exon_pr, strandedness=False).df
        if joined.empty:
            pp_rows, pseudo_rows = set(), set()
        else:
            pp_rows = set(joined.loc[joined["gene_type"] == PROCESSED_PSEUDOGENE, "row"])
            pseudo_rows = set(joined.loc[
                joined["gene_type"].str.contains("pseudogene", na=False), "row"])
        frame["is_pp"] = frame.index.isin(pp_rows)
        frame["is_pseudo"] = frame.index.isin(pseudo_rows)

        by_read = frame.groupby("qname", sort=False)
        target_at_max = (frame["match_target"] & frame["at_max"]).groupby(
            frame["qname"], sort=False).any()
        pp_at_max = (frame["is_pp"] & frame["at_max"]).groupby(
            frame["qname"], sort=False).any()
        pseudo_at_max = (frame["is_pseudo"] & frame["at_max"]).groupby(
            frame["qname"], sort=False).any()
        del by_read
        n_target_at_max = int(target_at_max.sum())
        n_pp_tie = int((target_at_max & pp_at_max).sum())
        n_any_pseudo_tie = int((target_at_max & pseudo_at_max).sum())

    n_other = n_txome_assigned - n_unique - n_pp_tie
    if n_other < 0:
        raise FateError(
            "%s %s: genome-unique (%d) + pseudogene-tie (%d) exceeds the "
            "transcriptome-assigned total (%d)"
            % (sample, tid, n_unique, n_pp_tie, n_txome_assigned))
    partition = {"genome_unique": n_unique,
                 "genome_multi_pseudogene_tie": n_pp_tie,
                 "other": n_other}
    if sum(partition.values()) != n_txome_assigned:
        raise FateError("%s %s: the three categories sum to %d, not %d"
                        % (sample, tid, sum(partition.values()), n_txome_assigned))

    wide = {
        "sample": sample, "tid": tid,
        "n_txome_assigned": n_txome_assigned,
        "n_txome_all_primary": all_primary,
        "n_shared_both_bams": n_shared,
        "n_genome_unique": n_unique,
        "n_genome_multi": n_multi,
        "n_top_on_target": n_target_at_max,
        "n_top_on_target_and_pp_tie": n_pp_tie,
        "n_top_on_target_and_any_pseudo_tie": n_any_pseudo_tie,
        "pct_pp_tie_of_assigned": round(100.0 * n_pp_tie / n_txome_assigned, 4),
        "pct_pp_tie_of_shared": round(100.0 * n_pp_tie / n_shared, 4) if n_shared else 0.0,
        "pct_pp_tie_of_genome_multi": round(100.0 * n_pp_tie / n_multi, 4) if n_multi else 0.0,
    }
    return wide, partition

def compute_fates(sample, genome_bam, txome_bam, gene_ids=(), transcript_ids=(),
                  coverage=None):
    """The whole computation: resolve targets, two BAM passes, one row set per transcript."""
    concordance_lib, _mm = load_libraries()
    payload = concordance_lib.build_transcript_table()
    table, base2ver = payload["table"], payload["base2ver"]

    tids = resolve_transcripts(table, gene_ids, transcript_ids, coverage)
    if not tids:
        raise FateError("no transcripts requested")
    for tid in tids:
        if tid not in table:
            raise FateError("%s is not in the APPRIS transcript table" % tid)

    populations, all_primary = collect_txome_populations(
        concordance_lib, txome_bam, base2ver, tids)
    targets = set()
    for population in populations.values():
        targets |= set(population)
    genome_present, genome_unique, records = collect_genome_records(genome_bam, targets)
    exon_pr = concordance_lib.load_exon_gene_pr()

    display_names, gene_ids_out = {}, {}
    for tid in tids:
        name = table[tid].get("gene_name", "")
        if isinstance(name, float) or name is None:
            name = ""
        display_names[tid] = str(name)
        gene_ids_out[tid] = str(table[tid].get("gene_id", "") or "")
    if coverage is not None:
        for tid in tids:
            try:
                index = coverage.index_of_transcript(tid)
            except Exception:
                continue
            info = coverage.transcript_info(index)
            display_names[tid] = info["gene_name"]
            gene_ids_out[tid] = info["gene_id"]

    wide_rows, tidy_rows = [], []
    for tid in tids:
        wide, partition = classify_transcript(
            concordance_lib, tid, table[tid], populations[tid], all_primary[tid],
            genome_present, genome_unique, records, exon_pr, sample)
        wide_rows.append(wide)
        for category in CATEGORIES:
            tidy_rows.append({
                "sample": sample,
                "transcript_id": tid,
                "gene_id": gene_ids_out[tid],
                "gene_name": display_names[tid],
                "category": category,
                "n_reads": partition[category],
                "pct_of_txome_assigned": round(
                    100.0 * partition[category] / wide["n_txome_assigned"], 6),
            })
    return pd.DataFrame(wide_rows)[WIDE_COLUMNS], pd.DataFrame(tidy_rows)[TIDY_COLUMNS]
