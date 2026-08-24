#!/usr/bin/env python3
"""Transcript region overlays: UTR5 / CDS / UTR3, and the derived ribopy bins."""
from __future__ import annotations

import gzip
import re
from pathlib import Path

UTR5, CDS, UTR3 = "UTR5", "CDS", "UTR3"
CANONICAL_LABELS = (UTR5, CDS, UTR3)

RIBO_BIN_LABELS = ("UTR5_OUTER", "START_WINDOW", "CDS_CORE", "STOP_WINDOW", "UTR3_OUTER")
RIBOPY_ALIASES = {
    "UTR5_OUTER": "UTR5",
    "START_WINDOW": "UTR5J",
    "CDS_CORE": "CDS",
    "STOP_WINDOW": "UTR3J",
    "UTR3_OUTER": "UTR3",
}
DEFAULT_LEFT_SPAN = 35
DEFAULT_RIGHT_SPAN = 10

STOP_CODON_NT = 3

_RE_UTR5 = re.compile(r"\|UTR5:(\d+)-(\d+)\|")
_RE_CDS = re.compile(r"\|CDS:(\d+)-(\d+)\|")
_RE_UTR3 = re.compile(r"\|UTR3:(\d+)-(\d+)\|")
_RE_TID = re.compile(r'transcript_id "([^"]+)"')

class RegionError(RuntimeError):
    """Raised when the region sources disagree or an invariant fails."""

def parse_reference_headers(appris_lengths: Path) -> dict:
    """Parse `appris_human_v2_transcript_lengths.tsv` -> {transcript_id: identity + raw regions}.

    Header intervals stay 1-based inclusive (start, end), or None when the field is absent.
    """
    out = {}
    with open(appris_lengths) as fh:
        for lineno, line in enumerate(fh, 1):
            line = line.rstrip("\n")
            if not line:
                continue
            try:
                name, length = line.split("\t")
            except ValueError:
                raise RegionError(
                    "%s:%d is not <name>TAB<length>: %.80r" % (appris_lengths, lineno, line))
            fields = name.split("|")
            if len(fields) < 7:
                raise RegionError(
                    "%s:%d header has %d pipe-delimited fields, expected at least 7 "
                    "(transcript_id|gene_id|havana_gene|havana_tx|tx_name|gene_name|length|...): "
                    "%.80r" % (appris_lengths, lineno, len(fields), name))
            tid = fields[0]
            m_cds = _RE_CDS.search(name)
            if m_cds is None:
                raise RegionError(
                    "%s:%d (%s) has no |CDS:start-end| field. Every transcript in this "
                    "reference must define a CDS." % (appris_lengths, lineno, tid))
            m5, m3 = _RE_UTR5.search(name), _RE_UTR3.search(name)
            out[tid] = {
                "transcript_id": tid,
                "gene_id": fields[1],
                "transcript_name": fields[4],
                "gene_name": fields[5],
                "transcript_len": int(length),
                "reference_name": name,
                "hdr_utr5": (int(m5.group(1)), int(m5.group(2))) if m5 else None,
                "hdr_cds": (int(m_cds.group(1)), int(m_cds.group(2))),
                "hdr_utr3": (int(m3.group(1)), int(m3.group(2))) if m3 else None,
            }
    if not out:
        raise RegionError("%s contained no transcripts" % appris_lengths)
    return out

def parse_stop_codon_transcripts(gtf: Path, wanted: set) -> set:
    """Transcript IDs with a GENCODE `stop_codon` feature -- the authoritative relocation test."""
    found = set()
    opener = gzip.open if str(gtf).endswith(".gz") else open
    with opener(gtf, "rt") as fh:
        for line in fh:
            if line[0] == "#":
                continue
            fields = line.split("\t", 9)
            if len(fields) < 9 or fields[2] != "stop_codon":
                continue
            m = _RE_TID.search(fields[8])
            if m and m.group(1) in wanted:
                found.add(m.group(1))
    return found

def parse_actual_regions_bed(bed_path: Path) -> dict:
    """Parse `appris_human_v2_actual_regions.bed` -> {transcript_id: {label: (start, end)}}.

    BED6, 0-based half-open, already in transcript coordinates.
    """
    out = {}
    with open(bed_path) as fh:
        for lineno, line in enumerate(fh, 1):
            line = line.rstrip("\n")
            if not line or line.startswith(("#", "track", "browser")):
                continue
            fields = line.split("\t")
            if len(fields) < 4:
                raise RegionError(
                    "%s:%d has %d columns, expected at least 4 (name, start, end, label)"
                    % (bed_path, lineno, len(fields)))
            tid = fields[0].split("|", 1)[0]
            label = fields[3]
            if label not in CANONICAL_LABELS:
                raise RegionError(
                    "%s:%d has unknown region label %r; expected one of %s"
                    % (bed_path, lineno, label, ", ".join(CANONICAL_LABELS)))
            out.setdefault(tid, {})[label] = (int(fields[1]), int(fields[2]))
    if not out:
        raise RegionError("%s contained no regions" % bed_path)
    return out

def derive_normalized_regions(header: dict, has_annotated_stop: bool) -> dict:
    """{label: (start, end)} 0-based half-open, tiling [0, L); an annotated stop codon
    moves into UTR3. A region that does not exist is ABSENT, never zero-length."""
    length = header["transcript_len"]
    cds_start_1, cds_end_1 = header["hdr_cds"]
    start, end = cds_start_1 - 1, cds_end_1

    if not 0 <= start < end <= length:
        raise RegionError(
            "%s: header CDS %d-%d is not inside [1, %d]"
            % (header["transcript_id"], cds_start_1, cds_end_1, length))

    if has_annotated_stop:
        if end - start < STOP_CODON_NT:
            raise RegionError(
                "%s: a stop codon is annotated but the header CDS is only %d nt, so 3 nt "
                "cannot be relocated" % (header["transcript_id"], end - start))
        end -= STOP_CODON_NT

    regions = {}
    if start > 0:
        regions[UTR5] = (0, start)
    regions[CDS] = (start, end)
    if end < length:
        regions[UTR3] = (end, length)
    return regions

def check_tiling(transcript_id: str, regions: dict, length: int) -> None:
    """Regions must partition [0, L) with no gap and no overlap."""
    spans = sorted(regions.values())
    if not spans:
        raise RegionError("%s: no regions at all" % transcript_id)
    if spans[0][0] != 0 or spans[-1][1] != length:
        raise RegionError(
            "%s: regions span [%d, %d), expected [0, %d)"
            % (transcript_id, spans[0][0], spans[-1][1], length))
    for (a_start, a_end), (b_start, b_end) in zip(spans, spans[1:]):
        if a_end != b_start:
            raise RegionError(
                "%s: regions are not contiguous -- [%d, %d) then [%d, %d)"
                % (transcript_id, a_start, a_end, b_start, b_end))
    for start, end in spans:
        if start >= end:
            raise RegionError("%s: empty or inverted region [%d, %d)"
                              % (transcript_id, start, end))

def build_regions(headers: dict, stop_codon_ids: set, bed: dict = None) -> tuple:
    """Per-transcript region overlays (raw + normalized) -> (rows, summary).

    A supplied `bed` is a cross-check, never a fallback: any disagreement raises.
    """
    rows = []
    n_relocated = n_not_relocated = 0
    n_bed_checked = n_bed_missing = 0
    disagreements = []

    for tid in sorted(headers):
        header = headers[tid]
        length = header["transcript_len"]
        has_stop = tid in stop_codon_ids
        normalized = derive_normalized_regions(header, has_stop)
        check_tiling(tid, normalized, length)
        n_relocated += has_stop
        n_not_relocated += not has_stop

        bed_regions = bed.get(tid) if bed is not None else None
        if bed is not None:
            if bed_regions is None:
                n_bed_missing += 1
            else:
                n_bed_checked += 1
                if bed_regions != normalized:
                    disagreements.append((tid, dict(sorted(normalized.items())),
                                          dict(sorted(bed_regions.items()))))

        raw_header = {UTR5: header["hdr_utr5"], CDS: header["hdr_cds"], UTR3: header["hdr_utr3"]}
        for label in CANONICAL_LABELS:
            if label not in normalized:
                continue
            start, end = normalized[label]
            hdr = raw_header[label]
            bed_iv = (bed_regions or {}).get(label)
            rows.append({
                "transcript_id": tid,
                "label": label,
                "raw_header_start_1based": hdr[0] if hdr else -1,
                "raw_header_end_1based": hdr[1] if hdr else -1,
                "raw_bed_start": bed_iv[0] if bed_iv else -1,
                "raw_bed_end": bed_iv[1] if bed_iv else -1,
                "start": start,
                "end": end,
                "source": "bed" if bed_iv == (start, end) else "derived",
            })

    if disagreements:
        detail = "\n".join(
            "    %s  derived %s  bed %s" % (tid, d, b) for tid, d, b in disagreements[:10])
        raise RegionError(
            "the derived regions disagree with %d transcript(s) in the supplied "
            "actual_regions.bed. This means the relocation rule does not describe this "
            "reference and must not be used for it.\n%s%s"
            % (len(disagreements), detail,
               "\n    ... and %d more" % (len(disagreements) - 10)
               if len(disagreements) > 10 else ""))

    summary = {
        "n_transcripts": len(headers),
        "n_region_rows": len(rows),
        "n_stop_relocated": n_relocated,
        "n_no_annotated_stop": n_not_relocated,
        "n_checked_against_bed": n_bed_checked,
        "n_absent_from_bed": n_bed_missing,
        "stop_codon_assignment": "utr3",
        "relocation_rule": "gtf_stop_codon_feature" if stop_codon_ids else "none",
    }
    return rows, summary

def heuristic_stop_codon_ids(headers: dict) -> set:
    """No-GTF fallback: header UTR3 present, or CDS length % 3 == 0.

    Matches the GTF stop_codon set on the shipped v2 reference, but IS a heuristic
    (divisibility alone misses 5'-incomplete annotations that carry a UTR3).
    """
    out = set()
    for tid, header in headers.items():
        cds_start, cds_end = header["hdr_cds"]
        if header["hdr_utr3"] is not None or (cds_end - cds_start + 1) % 3 == 0:
            out.add(tid)
    return out

def build_ribo_region_bins(headers: dict, left_span: int, right_span: int) -> list:
    """ribopy's five-way binning (port of `region_lib.classify` boundaries), derived, not annotation.

    Uses the STOP-INCLUSIVE header CDS end -- ribopy's convention, deliberately not harmonised
    with the canonical regions above; bins are clipped to [0, L) and empty bins omitted.
    """
    if left_span < 0 or right_span < 0:
        raise RegionError("left_span and right_span must be >= 0, got %d / %d"
                          % (left_span, right_span))
    rows = []
    for tid in sorted(headers):
        header = headers[tid]
        length = header["transcript_len"]
        start_site = header["hdr_cds"][0] - 1
        stop_site = header["hdr_cds"][1]
        bounds = [
            ("UTR5_OUTER", 0, start_site - left_span),
            ("START_WINDOW", start_site - left_span, start_site + right_span + 1),
            ("CDS_CORE", start_site + right_span + 1, stop_site - left_span),
            ("STOP_WINDOW", stop_site - left_span, stop_site + right_span + 1),
            ("UTR3_OUTER", stop_site + right_span + 1, length),
        ]
        for label, start, end in bounds:
            start, end = max(0, min(start, length)), max(0, min(end, length))
            if start >= end:
                continue
            rows.append({
                "transcript_id": tid,
                "label": label,
                "ribopy_alias": RIBOPY_ALIASES[label],
                "start": start,
                "end": end,
            })
    return rows

def ribo_bin_provenance(left_span: int, right_span: int, parameter_source: str) -> dict:
    """The attributes stored on /ribo_region_bins. Records HOW the spans were chosen."""
    return {
        "algorithm": "ribopy_get_extended_boundaries",
        "left_span": int(left_span),
        "right_span": int(right_span),
        "start_site_source": "header_cds_start_minus_1",
        "stop_site_source": "header_cds_end_stop_inclusive",
        "parameter_source": parameter_source,
        "note": "Derived junction windows, NOT annotation. Labels carry ribopy_alias for "
                "ribopy's UTR5/UTR5J/CDS/UTR3J/UTR3 names.",
    }

