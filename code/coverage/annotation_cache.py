#!/usr/bin/env python3
"""The sample-independent half of a coverage build, computed once and reused."""
from __future__ import annotations

import hashlib
import pickle
import time
from pathlib import Path

CACHE_VERSION = 1

SOURCE_MODULES = ("transcript_coords.py", "transcript_regions.py", "annotation_cache.py")

def log(message):
    print("[annotation] %s" % message, flush=True)

def _sha256(path):
    digest = hashlib.sha256()
    with open(path, "rb") as handle:
        for chunk in iter(lambda: handle.read(1 << 20), b""):
            digest.update(chunk)
    return digest.hexdigest()

def fingerprint(gtf, appris, regions, left_span, right_span):
    """A content digest of everything that determines the bundle."""
    here = Path(__file__).resolve().parent
    parts = ["version=%d" % CACHE_VERSION,
             "left_span=%d" % left_span, "right_span=%d" % right_span]
    inputs = {}
    for label, path in (("gtf", gtf), ("appris", appris), ("regions", regions)):
        if path is None:
            parts.append("%s=absent" % label)
            continue
        digest = _sha256(path)
        inputs[label] = {"name": Path(path).name, "sha256": digest,
                         "bytes": Path(path).stat().st_size}
        parts.append("%s=%s" % (label, digest))
    for name in SOURCE_MODULES:
        module = here / name
        parts.append("%s=%s" % (name, _sha256(module) if module.exists() else "absent"))
    return hashlib.sha256("\n".join(parts).encode()).hexdigest(), inputs

def build(gtf, appris, regions, left_span, right_span):
    """Parse the annotation into the reusable coordinate bundle."""
    import pandas as pd
    import transcript_coords
    import transcript_regions

    started = time.time()
    log("parsing headers and GTF features")
    headers = transcript_regions.parse_reference_headers(appris)
    features = transcript_coords.parse_gtf_features(gtf, set(headers))
    coords = transcript_coords.build_transcript_coords(features, headers)
    cds_table = transcript_coords.build_cds_exon_table(features)
    transcripts, exons = coords["transcripts"], coords["exons"]
    log("  %d transcripts, %d exons, %d positions"
        % (len(transcripts), len(exons), coords["n_positions"]))

    stop_ids = transcript_regions.parse_stop_codon_transcripts(gtf, set(headers))
    bed = transcript_regions.parse_actual_regions_bed(regions) if regions else None
    region_rows, region_summary = transcript_regions.build_regions(headers, stop_ids, bed)
    bin_rows = transcript_regions.build_ribo_region_bins(headers, left_span, right_span)
    log("  %d region rows (%d relocated, %d without a stop), %d ribo bins"
        % (len(region_rows), region_summary["n_stop_relocated"],
           region_summary["n_no_annotated_stop"], len(bin_rows)))

    index_of_id = {tid: i for i, tid in enumerate(transcripts["transcript_id"])}
    regions_df = pd.DataFrame(region_rows)
    regions_df["transcript_index"] = regions_df["transcript_id"].map(index_of_id)
    regions_df = regions_df.sort_values(["transcript_index", "start"]).reset_index(drop=True)

    ribo_bins = pd.DataFrame(bin_rows) if bin_rows else pd.DataFrame(
        columns=["transcript_id", "label", "ribopy_alias", "start", "end"])
    if len(ribo_bins):
        ribo_bins["transcript_index"] = ribo_bins["transcript_id"].map(index_of_id)
        ribo_bins = ribo_bins.sort_values(
            ["transcript_index", "start"]).reset_index(drop=True)

    log("  built in %.1f s" % (time.time() - started))
    return {
        "headers": headers, "coords": coords, "cds_table": cds_table,
        "transcripts": transcripts, "exons": exons,
        "n_positions": coords["n_positions"],
        "regions": regions_df, "ribo_bins": ribo_bins,
        "region_summary": region_summary, "stop_ids": stop_ids,
        "index_of_id": index_of_id,
        "index_of_base": {tid.split(".", 1)[0]: i for tid, i in index_of_id.items()},
    }

def load_or_build(cache_path, gtf, appris, regions, left_span, right_span):
    """Reuse the cache when its fingerprint matches; otherwise rebuild and write it.

    Returns (bundle, reused). A mismatch is not an error: it means an input changed, and
    rebuilding is the correct response.
    """
    digest, inputs = fingerprint(gtf, appris, regions, left_span, right_span)
    cache_path = Path(cache_path) if cache_path else None

    if cache_path and cache_path.exists():
        try:
            with open(cache_path, "rb") as handle:
                stored = pickle.load(handle)
        except Exception as exc:
            log("cache unreadable (%s); rebuilding" % exc)
        else:
            if stored.get("fingerprint") == digest:
                log("reusing %s" % cache_path)
                return stored["bundle"], True
            log("cache fingerprint differs from the current inputs; rebuilding")

    bundle = build(gtf, appris, regions, left_span, right_span)
    if cache_path:
        cache_path.parent.mkdir(parents=True, exist_ok=True)
        temporary = cache_path.with_suffix(cache_path.suffix + ".tmp-%d" % time.time_ns())
        with open(temporary, "wb") as handle:
            pickle.dump({"fingerprint": digest, "inputs": inputs,
                         "version": CACHE_VERSION, "bundle": bundle}, handle,
                        protocol=pickle.HIGHEST_PROTOCOL)
        temporary.replace(cache_path)
        log("wrote %s (%.1f MB)" % (cache_path, cache_path.stat().st_size / 1e6))
    return bundle, False

