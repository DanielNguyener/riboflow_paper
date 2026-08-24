#!/usr/bin/env python3
"""The shared transcript coordinate: complete GTF exons, ordered 5' to 3'."""
from __future__ import annotations

import gzip
import hashlib
import re
from pathlib import Path

import numpy as np
import pandas as pd

_RE_TID = re.compile(r'transcript_id "([^"]+)"')
_RE_GID = re.compile(r'gene_id "([^"]+)"')
_RE_GNAME = re.compile(r'gene_name "([^"]+)"')
_FEATURES = ("exon", "CDS")

CDS_TABLE_COLUMNS = (
    "Chromosome", "Start", "End", "Strand", "Phase", "transcript_id", "gene_id",
    "gene_name", "exon_len", "order_key", "exon_index", "cds_cum_start", "cds_cum_end")

TRANSCRIPT_COLUMNS = (
    "transcript_id", "gene_id", "transcript_name", "gene_name", "chrom", "strand",
    "transcript_len", "n_exons", "cds_len_gtf", "n_cds_exons",
    "coverage_offset", "exon_offset",
)
EXON_COLUMNS = (
    "transcript_index", "exon_index", "chrom", "g_start", "g_end", "tx_start", "tx_end",
)

class CoordinateError(RuntimeError):
    """Raised when the exon map cannot be built or fails an invariant."""

def sha256_file(path: Path) -> str:
    digest = hashlib.sha256()
    with open(path, "rb") as fh:
        for chunk in iter(lambda: fh.read(1 << 20), b""):
            digest.update(chunk)
    return digest.hexdigest()

def parse_gtf_features(gtf: Path, wanted: set) -> dict:
    """One GTF pass -> {tid: {"exon": [...], "CDS": [...]}}, coordinates 0-based half-open.

    List order is GTF file order, NOT 5'->3'; `build_transcript_coords` does the strand-aware sort.
    """
    if not gtf.exists():
        raise CoordinateError("GTF not found: %s" % gtf)
    collected = {}
    opener = gzip.open if str(gtf).endswith(".gz") else open
    with opener(gtf, "rt") as handle:
        for line in handle:
            if line[0] == "#":
                continue
            fields = line.split("\t", 9)
            if len(fields) < 9:
                continue
            feature = fields[2]
            if feature not in _FEATURES:
                continue
            match = _RE_TID.search(fields[8])
            if match is None:
                continue
            tid = match.group(1)
            if tid not in wanted:
                continue
            entry = collected.setdefault(tid, {"exon": [], "CDS": []})
            if feature == "CDS":
                entry["CDS"].append((
                    fields[0], int(fields[3]) - 1, int(fields[4]), fields[6],
                    int(fields[7]) if fields[7] != "." else 0,
                    _RE_GID.search(fields[8]).group(1) if _RE_GID.search(fields[8]) else "",
                    _RE_GNAME.search(fields[8]).group(1)
                    if _RE_GNAME.search(fields[8]) else ""))
            else:
                entry["exon"].append(
                    (fields[0], int(fields[3]) - 1, int(fields[4]), fields[6]))
    return collected

def build_cds_exon_table(features: dict) -> pd.DataFrame:
    """The CDS-exon table (CDS pieces, not whole exons), one row per piece in a FIXED order.

    Row order is load-bearing: stage-1 tie-breaks read the first matching row, so the sort
    must match `bam_inputs.build_cds_table()["cds"]` exactly (asserted equal).
    """
    rows = []
    for tid in sorted(features):
        for chrom, start, end, strand, phase, gene_id, gene_name in features[tid]["CDS"]:
            rows.append((chrom, start, end, strand, phase, tid, gene_id, gene_name))
    cds = pd.DataFrame(rows, columns=[
        "Chromosome", "Start", "End", "Strand", "Phase", "transcript_id",
        "gene_id", "gene_name"])
    cds["exon_len"] = cds["End"] - cds["Start"]
    cds["order_key"] = np.where(cds["Strand"] == "+", cds["Start"], -cds["Start"])
    cds = cds.sort_values(["transcript_id", "order_key"]).reset_index(drop=True)
    grouped = cds.groupby("transcript_id", sort=False)
    cds["exon_index"] = grouped.cumcount()
    cds["cds_cum_start"] = grouped["exon_len"].cumsum() - cds["exon_len"]
    cds["cds_cum_end"] = cds["cds_cum_start"] + cds["exon_len"]
    return cds[list(CDS_TABLE_COLUMNS)]

def _order_exons(tid: str, exons: list) -> tuple:
    """Sort exons 5'->3' and validate they are same-chrom, same-strand, non-overlapping."""
    chroms = {e[0] for e in exons}
    strands = {e[3] for e in exons}
    if len(chroms) != 1:
        raise CoordinateError(
            "%s has exons on %d chromosomes (%s); a transcript coordinate is undefined"
            % (tid, len(chroms), ", ".join(sorted(chroms))))
    if len(strands) != 1:
        raise CoordinateError(
            "%s has exons on both strands (%s)" % (tid, ", ".join(sorted(strands))))
    chrom, strand = chroms.pop(), strands.pop()
    if strand not in ("+", "-"):
        raise CoordinateError("%s has strand %r, expected '+' or '-'" % (tid, strand))

    ordered = sorted(exons, key=lambda e: e[1], reverse=(strand == "-"))
    genomic = sorted((e[1], e[2]) for e in ordered)
    for (a_start, a_end), (b_start, b_end) in zip(genomic, genomic[1:]):
        if a_end > b_start:
            raise CoordinateError(
                "%s has overlapping exons [%d, %d) and [%d, %d)"
                % (tid, a_start, a_end, b_start, b_end))
    for start, end in genomic:
        if start >= end:
            raise CoordinateError("%s has an empty exon [%d, %d)" % (tid, start, end))
    return chrom, strand, ordered

def build_transcript_coords(features: dict, headers: dict) -> dict:
    """Build the transcripts and exons tables, in sorted-transcript_id storage order.

    Storage order is load-bearing for the pooled-Pearson reconstruction downstream.
    """
    missing = sorted(set(headers) - set(features))
    if missing:
        raise CoordinateError(
            "%d selected transcript(s) have no GTF `exon` feature, e.g. %s. The GTF and the "
            "transcriptome reference describe different annotations."
            % (len(missing), ", ".join(missing[:5])))

    transcript_rows, exon_rows = [], []
    coverage_offset = 0
    mismatches = []

    for index, tid in enumerate(sorted(headers)):
        header = headers[tid]
        entry = features[tid]
        chrom, strand, ordered = _order_exons(tid, entry["exon"])

        exon_offset = len(exon_rows)
        spliced = 0
        for exon_index, (_chrom, g_start, g_end) in enumerate(
                (e[0], e[1], e[2]) for e in ordered):
            length = g_end - g_start
            exon_rows.append((index, exon_index, chrom, g_start, g_end,
                              spliced, spliced + length))
            spliced += length

        if spliced != header["transcript_len"]:
            mismatches.append((tid, spliced, header["transcript_len"]))

        cds = entry["CDS"]
        transcript_rows.append((
            tid, header["gene_id"], header["transcript_name"], header["gene_name"],
            chrom, strand, header["transcript_len"], len(ordered),
            sum(row[2] - row[1] for row in cds), len(cds),
            coverage_offset, exon_offset,
        ))
        coverage_offset += header["transcript_len"]

    if mismatches:
        detail = "\n".join("    %s  spliced %d  reference %d  (diff %+d)"
                           % (t, s, r, s - r) for t, s, r in mismatches[:10])
        raise CoordinateError(
            "spliced exon length != transcriptome reference length for %d transcript(s). "
            "The genome and transcriptome routes would not share a coordinate.\n%s%s"
            % (len(mismatches), detail,
               "\n    ... and %d more" % (len(mismatches) - 10)
               if len(mismatches) > 10 else ""))

    transcripts = pd.DataFrame(transcript_rows, columns=list(TRANSCRIPT_COLUMNS))
    exons = pd.DataFrame(exon_rows, columns=list(EXON_COLUMNS))
    _validate(transcripts, exons)
    return {
        "transcripts": transcripts,
        "exons": exons,
        "n_positions": int(coverage_offset),
    }

def _validate(transcripts: pd.DataFrame, exons: pd.DataFrame) -> None:
    """Invariants that must hold before anything downstream trusts this table."""
    if transcripts["transcript_id"].duplicated().any():
        dup = transcripts.loc[transcripts["transcript_id"].duplicated(), "transcript_id"]
        raise CoordinateError("duplicate transcript_id: %s" % ", ".join(dup.head(5)))

    ids = transcripts["transcript_id"].tolist()
    if ids != sorted(ids):
        raise CoordinateError(
            "transcripts are not in sorted transcript_id order; storage order is "
            "load-bearing for the pooled-Pearson reconstruction")

    expected = np.concatenate([[0], np.cumsum(transcripts["transcript_len"].to_numpy())[:-1]])
    if not np.array_equal(transcripts["coverage_offset"].to_numpy(), expected):
        raise CoordinateError("coverage_offset is not the running sum of transcript_len")

    spliced = np.bincount(
        exons["transcript_index"].to_numpy(),
        weights=(exons["g_end"] - exons["g_start"]).to_numpy(),
        minlength=len(transcripts)).astype(np.int64)
    if not np.array_equal(spliced, transcripts["transcript_len"].to_numpy()):
        bad = int(np.argmax(spliced != transcripts["transcript_len"].to_numpy()))
        raise CoordinateError(
            "per-transcript spliced exon length != transcript_len, first at index %d (%s): "
            "%d vs %d" % (bad, transcripts.at[bad, "transcript_id"], spliced[bad],
                          transcripts.at[bad, "transcript_len"]))

    starts = exons["tx_start"].to_numpy()
    ends = exons["tx_end"].to_numpy()
    idx = exons["transcript_index"].to_numpy()
    if not np.all(ends > starts):
        raise CoordinateError("some exon has tx_end <= tx_start")

    is_first = np.r_[True, idx[1:] != idx[:-1]]
    if not np.all(starts[is_first] == 0):
        raise CoordinateError("the first exon of some transcript does not start at tx 0")

    same_transcript = idx[1:] == idx[:-1]
    if not np.all(starts[1:][same_transcript] == ends[:-1][same_transcript]):
        raise CoordinateError("exon tx_start/tx_end are not contiguous within a transcript")

    exon_index = exons["exon_index"].to_numpy()
    if not np.all(exon_index[is_first] == 0):
        raise CoordinateError("some transcript's first exon does not have exon_index 0")
    if not np.all(exon_index[1:][same_transcript] == exon_index[:-1][same_transcript] + 1):
        raise CoordinateError("exon_index is not consecutive within a transcript")

def transcript_exons(coords: dict, transcript_index: int) -> pd.DataFrame:
    transcripts, exons = coords["transcripts"], coords["exons"]
    start = int(transcripts.at[transcript_index, "exon_offset"])
    count = int(transcripts.at[transcript_index, "n_exons"])
    return exons.iloc[start:start + count]

def tx_to_genomic(coords: dict, transcript_index: int, positions) -> np.ndarray:
    """Transcript positions -> genomic positions (0-based), strand-aware; out-of-range yields -1."""
    strand = coords["transcripts"].at[transcript_index, "strand"]
    exons = transcript_exons(coords, transcript_index)
    tx_start = exons["tx_start"].to_numpy()
    tx_end = exons["tx_end"].to_numpy()
    g_start = exons["g_start"].to_numpy()
    g_end = exons["g_end"].to_numpy()

    pos = np.asarray(positions, dtype=np.int64)
    out = np.full(pos.shape, -1, dtype=np.int64)
    inside = (pos >= 0) & (pos < tx_end[-1])
    if not inside.any():
        return out
    which = np.searchsorted(tx_start, pos[inside], side="right") - 1
    offset = pos[inside] - tx_start[which]
    out[inside] = (g_start[which] + offset if strand == "+"
                   else g_end[which] - 1 - offset)
    return out

def genomic_to_tx(coords: dict, transcript_index: int, positions) -> np.ndarray:
    """Genomic positions -> transcript positions, strand-aware. -1 when not on an exon."""
    strand = coords["transcripts"].at[transcript_index, "strand"]
    exons = transcript_exons(coords, transcript_index)
    tx_start = exons["tx_start"].to_numpy()
    g_start = exons["g_start"].to_numpy()
    g_end = exons["g_end"].to_numpy()

    pos = np.asarray(positions, dtype=np.int64)
    out = np.full(pos.shape, -1, dtype=np.int64)
    for i in range(len(g_start)):
        hit = (pos >= g_start[i]) & (pos < g_end[i])
        if not hit.any():
            continue
        offset = (pos[hit] - g_start[i] if strand == "+" else g_end[i] - 1 - pos[hit])
        out[hit] = tx_start[i] + offset
    return out

