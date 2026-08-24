#!/usr/bin/env python3
"""The shared-coverage HDF5 container: constants, writer, reader, validator.

Schema 3 stores only the four coverage arrays plus CDS bounds (layout: docs/hdf5_schema.md);
regions, counts and keys are derived on read, never stored.
"""
from __future__ import annotations

import argparse
import datetime as _datetime
import json
import os
import sys
from pathlib import Path

import numpy as np

SCHEMA = "riboflow_paper/shared-coverage/3"
SCHEMA_VERSION = 3

SIGNALS = ("genome_psite", "txome_psite", "genome_footprint", "txome_footprint")
PSITE_SIGNALS = ("genome_psite", "txome_psite")
FOOTPRINT_SIGNALS = ("genome_footprint", "txome_footprint")

ROUTES = ("genome", "transcriptome")
SIGNAL_ROUTE = {"genome_psite": "genome", "genome_footprint": "genome",
                "txome_psite": "transcriptome", "txome_footprint": "transcriptome"}
SIGNAL_MEASURE = {"genome_psite": "psite", "txome_psite": "psite",
                  "genome_footprint": "footprint", "txome_footprint": "footprint"}
ASSAYS = ("ribo", "rna")

COORDINATE_SYSTEM = "transcript_5p_to_3p"
REGION_LABELS = ("UTR5", "CDS", "UTR3")
#: cds_start == cds_end == NO_CDS marks a transcript without a CDS.
NO_CDS = -1

COVERAGE_DTYPE = np.int32
COVERAGE_MAX = int(np.iinfo(COVERAGE_DTYPE).max)

DEFAULT_CHUNK = 1 << 16
DEFAULT_GZIP_LEVEL = 9
DEFAULT_SHUFFLE = True

TRANSCRIPT_STR_COLUMNS = ("transcript_id", "gene_id", "gene_name")
TRANSCRIPT_INT_COLUMNS = ("transcript_len", "cds_start", "cds_end")
TRANSCRIPT_COLUMNS = TRANSCRIPT_STR_COLUMNS + TRANSCRIPT_INT_COLUMNS + ("coverage_offset",)

ROOT_ATTRS = ("schema", "schema_version", "sample", "assay", "routes", "coordinate_system",
              "psite_placement", "paper_cds_trim", "n_transcripts", "n_positions",
              "created_utc", "provenance")


class SchemaError(RuntimeError):
    """Raised when a coverage file violates the schema or an invariant."""


def _utc_now() -> str:
    return _datetime.datetime.utcnow().strftime("%Y-%m-%dT%H:%M:%SZ")


def _text(value) -> str:
    """HDF5 hands back `bytes` or `str` depending on how a value was written."""
    if isinstance(value, bytes):
        return value.decode()
    return value if isinstance(value, str) else str(value)


def sha256_file(path) -> str:
    import hashlib
    digest = hashlib.sha256()
    with open(path, "rb") as handle:
        for chunk in iter(lambda: handle.read(1 << 20), b""):
            digest.update(chunk)
    return digest.hexdigest()


CODE_VERSION_MODULES = ("build_shared_coverage.py", "coverage_schema.py",
                        "psite_placement.py", "transcript_coords.py",
                        "transcript_regions.py")


def code_version(here=None) -> dict:
    """{module: sha256} for the pipeline's own sources, plus a combined digest."""
    import hashlib

    here = Path(here) if here else Path(__file__).resolve().parent
    per_module = {}
    for name in CODE_VERSION_MODULES:
        path = here / name
        per_module[name] = sha256_file(path) if path.exists() else None
    combined = hashlib.sha256()
    for name in CODE_VERSION_MODULES:
        combined.update(("%s=%s\n" % (name, per_module[name])).encode())
    return {"modules": per_module, "combined_sha256": combined.hexdigest()}


def invocation(argv=None, record_paths=False) -> dict:
    """The command as invoked, `argv[0]` shortened to a repository-relative name and, unless
    `record_paths`, every path-like argument reduced to its basename."""
    argv = list(sys.argv if argv is None else argv)
    if not argv:
        return {"command": "", "paths_redacted": not record_paths}
    parts = ["code/coverage/" + Path(argv[0]).name]
    for token in argv[1:]:
        if record_paths or not token:
            parts.append(token)
            continue
        if token.startswith("-") and "=" in token:
            flag, _, value = token.partition("=")
            if os.sep in value:
                token = "%s=%s" % (flag, Path(value).name or value)
        elif not token.startswith("-") and os.sep in token:
            token = Path(token).name or token
        parts.append(token)
    return {"command": " ".join(parts), "paths_redacted": not record_paths}


# ── derived quantities ───────────────────────────────────────────────────────

STATE_NO_READS = "no_reads_assigned"
STATE_OUTSIDE_SLICE = "reads_outside_requested_slice"
STATE_COVERED = "covered"


def describe_coverage_state(n_events: int, slice_sum: int) -> str:
    """Which of the three coverage states holds for one transcript and one slice."""
    if n_events == 0:
        return STATE_NO_READS
    if slice_sum == 0:
        return STATE_OUTSIDE_SLICE
    return STATE_COVERED


def regions_from_cds(transcript_len: int, cds_start: int, cds_end: int) -> dict:
    """{label: (start, end)} tiling [0, transcript_len): UTR5 | CDS | UTR3, empty parts
    omitted. A transcript without a CDS has no regions."""
    if cds_start == NO_CDS:
        return {}
    regions = {}
    if cds_start > 0:
        regions["UTR5"] = (0, int(cds_start))
    regions["CDS"] = (int(cds_start), int(cds_end))
    if cds_end < transcript_len:
        regions["UTR3"] = (int(cds_end), int(transcript_len))
    return regions


def window_sums(values, coverage_offset, starts, ends) -> np.ndarray:
    """Sum `values` over one transcript-relative [start, end) window per transcript.

    end <= start sums to 0; out-of-range bounds clamp like a Python slice; exact int prefix sums.
    """
    values = np.asarray(values)
    prefix = np.concatenate([[0], np.cumsum(values, dtype=np.int64)])
    offset = np.asarray(coverage_offset, dtype=np.int64)
    starts, ends = np.asarray(starts, dtype=np.int64), np.asarray(ends, dtype=np.int64)
    lo = np.clip(offset + np.maximum(starts, 0), 0, values.size)
    hi = np.clip(offset + np.maximum(ends, starts), 0, values.size)
    hi = np.maximum(hi, lo)
    return prefix[hi] - prefix[lo]


# ── writer ───────────────────────────────────────────────────────────────────

class CoverageWriter:
    """Incremental, bounded, atomic writer.

        with CoverageWriter(path, sample=..., transcripts=..., provenance=...) as writer:
            writer.write_signal("genome_psite", array)      # one signal at a time
            ...

    `transcripts` is a DataFrame with TRANSCRIPT_COLUMNS, sorted by transcript_id.
    """

    def __init__(self, path, sample, transcripts, provenance, paper_cds_trim=15,
                 chunk=DEFAULT_CHUNK, gzip_level=DEFAULT_GZIP_LEVEL, shuffle=DEFAULT_SHUFFLE,
                 assay="ribo"):
        import h5py

        if assay not in ASSAYS:
            raise SchemaError("unknown assay %r; expected one of %s"
                              % (assay, ", ".join(ASSAYS)))
        missing = [c for c in TRANSCRIPT_COLUMNS if c not in transcripts.columns]
        if missing:
            raise SchemaError("transcript table lacks column(s): %s" % ", ".join(missing))
        self.assay = assay
        self.final_path = Path(path)
        self.tmp_path = self.final_path.with_name(
            "%s.tmp-%d-%s" % (self.final_path.name, os.getpid(),
                              _datetime.datetime.utcnow().strftime("%Y%m%dT%H%M%S")))
        self.sample = sample
        self.transcripts = transcripts
        self.n_positions = int(transcripts["transcript_len"].sum())
        self._written = set()

        self._refuse_stale_temporaries()
        self.final_path.parent.mkdir(parents=True, exist_ok=True)
        self.handle = h5py.File(self.tmp_path, "w")
        self._write_attrs(paper_cds_trim, provenance)
        self._write_transcripts()
        self._create_coverage(chunk, gzip_level, shuffle)

    def _refuse_stale_temporaries(self):
        stale = sorted(self.final_path.parent.glob(self.final_path.name + ".tmp-*")) \
            if self.final_path.parent.exists() else []
        if stale:
            raise SchemaError(
                "%d stale temporary file(s) beside %s, from an interrupted run:\n%s\n"
                "They are never reused, because a partial file cannot be distinguished "
                "from a complete one by inspection. Delete them and re-run."
                % (len(stale), self.final_path,
                   "\n".join("    %s" % p for p in stale)))

    def _write_attrs(self, paper_cds_trim, provenance):
        import psite_placement

        provenance = dict(provenance)
        provenance.setdefault("generation", invocation())
        provenance.setdefault("code_version", code_version())
        attrs = self.handle.attrs
        attrs["schema"] = SCHEMA
        attrs["schema_version"] = SCHEMA_VERSION
        attrs["sample"] = self.sample
        attrs["assay"] = self.assay
        attrs["routes"] = list(ROUTES)
        attrs["coordinate_system"] = COORDINATE_SYSTEM
        attrs["psite_placement"] = psite_placement.PSITE_PLACEMENT
        attrs["paper_cds_trim"] = int(paper_cds_trim)
        attrs["n_transcripts"] = int(len(self.transcripts))
        attrs["n_positions"] = self.n_positions
        attrs["created_utc"] = _utc_now()
        attrs["provenance"] = json.dumps(provenance, indent=2, sort_keys=True)

    def _write_transcripts(self):
        group = self.handle.create_group("transcripts")
        frame = self.transcripts
        n = len(frame)
        for column in TRANSCRIPT_STR_COLUMNS:
            values = frame[column].astype(str).to_numpy()
            width = max((len(v.encode()) for v in values), default=1) or 1
            group.create_dataset(column, data=values.astype("S%d" % width))
        for column in TRANSCRIPT_INT_COLUMNS:
            group.create_dataset(column, data=frame[column].to_numpy(dtype=np.int32))
        group.create_dataset("coverage_offset",
                             data=frame["coverage_offset"].to_numpy(dtype=np.int64))
        del n

    def _create_coverage(self, chunk, gzip_level, shuffle):
        group = self.handle.create_group("coverage")
        chunk = int(min(chunk, max(self.n_positions, 1)))
        for name in SIGNALS:
            group.create_dataset(
                name, shape=(self.n_positions,), dtype=COVERAGE_DTYPE,
                chunks=(chunk,), compression="gzip", compression_opts=gzip_level,
                shuffle=shuffle, fillvalue=0)

    def write_signal(self, name, values, chunk_positions=1 << 22):
        """Write one full-coordinate array, in chunks, with a range check (any integer
        dtype in; int32 on disk)."""
        if name not in SIGNALS:
            raise SchemaError("unknown signal %r; expected one of %s"
                              % (name, ", ".join(SIGNALS)))
        values = np.asarray(values)
        if values.shape != (self.n_positions,):
            raise SchemaError(
                "signal %s has shape %s, expected (%d,) -- the sum of transcript_len"
                % (name, values.shape, self.n_positions))
        if values.size:
            low, high = int(values.min()), int(values.max())
            if low < 0:
                raise SchemaError("signal %s has a negative value (%d)" % (name, low))
            if high > COVERAGE_MAX:
                raise SchemaError(
                    "signal %s reaches %d, which overflows int32 (max %d). Coverage is "
                    "stored as int32; this needs a schema change, not a silent cast."
                    % (name, high, COVERAGE_MAX))
        dataset = self.handle["coverage"][name]
        for start in range(0, self.n_positions, chunk_positions):
            stop = min(start + chunk_positions, self.n_positions)
            dataset[start:stop] = values[start:stop].astype(COVERAGE_DTYPE, copy=False)
        self._written.add(name)

    def finalize(self):
        """Validate the temporary file, then atomically move it into place."""
        missing = [s for s in SIGNALS if s not in self._written]
        if missing:
            raise SchemaError("signal(s) never written: %s" % ", ".join(missing))
        self.handle.close()
        self.handle = None
        problems = validate_file(self.tmp_path)
        if problems:
            raise SchemaError(
                "the file failed validation and was NOT moved into place. It is still at "
                "%s for inspection.\n%s"
                % (self.tmp_path, "\n".join("  - %s" % p for p in problems)))
        os.replace(self.tmp_path, self.final_path)
        return self.final_path

    def abort(self):
        """Close and remove the temporary file. Nothing reaches the final path."""
        if self.handle is not None:
            self.handle.close()
            self.handle = None
        if self.tmp_path.exists():
            self.tmp_path.unlink()

    def __enter__(self):
        return self

    def __exit__(self, exc_type, exc, tb):
        if exc_type is not None:
            self.abort()
            return False
        self.finalize()
        return False


# ── validator ────────────────────────────────────────────────────────────────

def validate_file(path) -> list:
    """Structural validation. Returns a list of problems; empty means valid."""
    import h5py

    problems = []
    path = Path(path)
    if not path.exists():
        return ["file does not exist: %s" % path]
    try:
        handle = h5py.File(path, "r")
    except Exception as exc:
        return ["cannot open as HDF5: %s" % exc]

    with handle:
        schema = _text(handle.attrs.get("schema", ""))
        if schema != SCHEMA:
            problems.append(
                "schema is %r, expected %r. Earlier schemas stored a different P-site rule "
                "(1) or a different layout (2); rebuild with "
                "code/coverage/build_shared_coverage.py." % (schema, SCHEMA))
            return problems
        for attr in ROOT_ATTRS:
            if attr not in handle.attrs:
                problems.append("missing root attribute %r" % attr)
        for group in ("transcripts", "coverage"):
            if group not in handle:
                problems.append("missing group /%s" % group)
        if problems:
            return problems
        if _text(handle.attrs["coordinate_system"]) != COORDINATE_SYSTEM:
            problems.append("coordinate_system is %r, expected %r"
                            % (handle.attrs["coordinate_system"], COORDINATE_SYSTEM))
        if _text(handle.attrs["assay"]) not in ASSAYS:
            problems.append("assay is %r, expected one of %s"
                            % (handle.attrs["assay"], ", ".join(ASSAYS)))
        try:
            json.loads(_text(handle.attrs["provenance"]))
        except ValueError as exc:
            problems.append("provenance is not valid JSON: %s" % exc)

        n_transcripts = int(handle.attrs["n_transcripts"])
        n_positions = int(handle.attrs["n_positions"])
        transcripts = handle["transcripts"]
        for column in TRANSCRIPT_COLUMNS:
            if column not in transcripts:
                problems.append("missing /transcripts/%s" % column)
            elif transcripts[column].shape != (n_transcripts,):
                problems.append("/transcripts/%s has shape %s, expected (%d,)"
                                % (column, transcripts[column].shape, n_transcripts))
        if problems:
            return problems

        transcript_len = transcripts["transcript_len"][:].astype(np.int64)
        coverage_offset = transcripts["coverage_offset"][:]
        cds_start = transcripts["cds_start"][:].astype(np.int64)
        cds_end = transcripts["cds_end"][:].astype(np.int64)
        if int(transcript_len.sum()) != n_positions:
            problems.append("sum(transcript_len) = %d but n_positions = %d"
                            % (transcript_len.sum(), n_positions))
        expected = np.concatenate([[0], np.cumsum(transcript_len)[:-1]]) \
            if n_transcripts else np.zeros(0, dtype=np.int64)
        if not np.array_equal(coverage_offset, expected):
            problems.append("coverage_offset is not the running sum of transcript_len")
        ids = [_text(s) for s in transcripts["transcript_id"][:]]
        if ids != sorted(ids):
            problems.append("transcripts are not in sorted transcript_id order "
                            "(storage order is load-bearing for the pooled-Pearson "
                            "reconstruction)")
        if len(set(ids)) != len(ids):
            problems.append("duplicate transcript_id")
        has_cds = cds_start != NO_CDS
        bad = has_cds & ((cds_start < 0) | (cds_end < cds_start) | (cds_end > transcript_len))
        if bad.any():
            first = int(np.argmax(bad))
            problems.append("CDS bounds outside [0, transcript_len), first at transcript "
                            "%d (%s): [%d, %d) of %d"
                            % (first, ids[first], cds_start[first], cds_end[first],
                               transcript_len[first]))
        if (~has_cds & (cds_end != NO_CDS)).any():
            problems.append("a transcript without a CDS must have cds_start == cds_end == -1")

        for name in SIGNALS:
            if name not in handle["coverage"]:
                problems.append("missing /coverage/%s" % name)
                continue
            dataset = handle["coverage"][name]
            if dataset.shape != (n_positions,):
                problems.append("/coverage/%s has shape %s, expected (%d,)"
                                % (name, dataset.shape, n_positions))
            if dataset.dtype != COVERAGE_DTYPE:
                problems.append("/coverage/%s has dtype %s, expected %s"
                                % (name, dataset.dtype, COVERAGE_DTYPE))
    return problems


# ── reader ───────────────────────────────────────────────────────────────────

class CoverageFile:
    """Read access with the schema checked on open."""

    def __init__(self, path):
        import h5py

        self.path = Path(path)
        problems = validate_file(self.path)
        if problems:
            raise SchemaError("%s failed validation:\n%s"
                              % (self.path, "\n".join("  - %s" % p for p in problems)))
        self.handle = h5py.File(self.path, "r")
        attrs = self.handle.attrs
        self.sample = _text(attrs["sample"])
        self.assay = _text(attrs["assay"])
        self.routes = tuple(_text(r) for r in attrs["routes"])
        self.coordinate_system = _text(attrs["coordinate_system"])
        self.psite_placement = _text(attrs["psite_placement"])
        self.schema_version = int(attrs["schema_version"])
        self.trim = int(attrs["paper_cds_trim"])
        self.n_transcripts = int(attrs["n_transcripts"])
        self.n_positions = int(attrs["n_positions"])
        group = self.handle["transcripts"]
        self._ids = [_text(s) for s in group["transcript_id"][:]]
        self._gene_ids = [_text(s) for s in group["gene_id"][:]]
        self._gene_names = [_text(s) for s in group["gene_name"][:]]
        self._len = group["transcript_len"][:].astype(np.int64)
        self._offset = group["coverage_offset"][:]
        self._cds_start = group["cds_start"][:].astype(np.int64)
        self._cds_end = group["cds_end"][:].astype(np.int64)

    @property
    def provenance(self) -> dict:
        return json.loads(_text(self.handle.attrs["provenance"]))

    def identity(self) -> dict:
        """Everything a consumer needs to decide whether this file is the right one."""
        return {
            "path": str(self.path),
            "schema": _text(self.handle.attrs["schema"]),
            "schema_version": self.schema_version,
            "sample": self.sample,
            "assay": self.assay,
            "routes": list(self.routes),
            "coordinate_system": self.coordinate_system,
            "psite_placement": self.psite_placement,
            "paper_cds_trim": self.trim,
            "n_transcripts": self.n_transcripts,
            "created_utc": _text(self.handle.attrs["created_utc"]),
        }

    # ── per-transcript tables, as arrays ──────────────────────────────────────
    @property
    def transcript_ids(self):
        return list(self._ids)

    @property
    def gene_names(self):
        return list(self._gene_names)

    @property
    def transcript_len(self):
        return self._len

    @property
    def coverage_offset(self):
        return self._offset

    @property
    def cds_start(self):
        return self._cds_start

    @property
    def cds_end(self):
        return self._cds_end

    # ── identifier resolution ─────────────────────────────────────────────────
    def index_of_transcript(self, transcript_id: str) -> int:
        """Exact or unversioned transcript id -> row index. Raises with candidates."""
        if transcript_id in self._ids:
            return self._ids.index(transcript_id)
        base = transcript_id.split(".", 1)[0]
        hits = [i for i, t in enumerate(self._ids) if t.split(".", 1)[0] == base]
        if len(hits) == 1:
            return hits[0]
        if not hits:
            raise SchemaError("transcript %r is not in %s" % (transcript_id, self.path))
        raise SchemaError(
            "transcript %r is ambiguous -- %d versions present: %s"
            % (transcript_id, len(hits), ", ".join(self._ids[i] for i in hits)))

    def resolve_gene(self, gene_id: str, transcript_id: str = None) -> int:
        """Gene id (versioned or not) -> row index. REFUSES to guess when ambiguous."""
        if transcript_id is not None:
            index = self.index_of_transcript(transcript_id)
            resolved = self._gene_ids[index]
            if gene_id and resolved.split(".", 1)[0] != gene_id.split(".", 1)[0]:
                raise SchemaError(
                    "transcript %s belongs to gene %s, not %s"
                    % (self._ids[index], resolved, gene_id))
            return index
        base = gene_id.split(".", 1)[0]
        hits = [i for i, g in enumerate(self._gene_ids) if g.split(".", 1)[0] == base]
        if not hits:
            raise SchemaError("gene %r is not in %s" % (gene_id, self.path))
        if len(hits) == 1:
            return hits[0]
        listing = "\n".join(
            "    %s  transcript_len=%d  cds_len=%d"
            % (self._ids[i], int(self._len[i]),
               max(int(self._cds_end[i] - self._cds_start[i]), 0)) for i in hits)
        raise SchemaError(
            "gene %r maps to %d eligible transcripts and no unique one resolves it. "
            "Pass --transcript-id to choose; this is not guessed.\n%s"
            % (gene_id, len(hits), listing))

    # ── one transcript ────────────────────────────────────────────────────────
    def transcript_info(self, index: int) -> dict:
        """Identity and geometry of one transcript; counts come from `event_counts`."""
        return {"transcript_id": self._ids[index], "gene_id": self._gene_ids[index],
                "gene_name": self._gene_names[index],
                "transcript_len": int(self._len[index]),
                "cds_start": int(self._cds_start[index]),
                "cds_end": int(self._cds_end[index]),
                "cds_len": max(int(self._cds_end[index] - self._cds_start[index]), 0)}

    def get_track(self, index: int, signal: str) -> np.ndarray:
        if signal not in SIGNALS:
            raise SchemaError("unknown signal %r; expected one of %s"
                              % (signal, ", ".join(SIGNALS)))
        start = int(self._offset[index])
        return self.handle["coverage"][signal][start:start + int(self._len[index])]

    def get_tracks(self, index: int, signals=SIGNALS) -> dict:
        return {name: self.get_track(index, name) for name in signals}

    def event_counts(self, index: int) -> dict:
        """{signal: sum over the whole transcript} -- P-site events, footprint bases."""
        return {name: int(self.get_track(index, name).sum(dtype=np.int64))
                for name in SIGNALS}

    def regions_of(self, index: int) -> dict:
        """{label: (start, end)}: UTR5 | CDS | UTR3 tiling the transcript."""
        return regions_from_cds(int(self._len[index]), int(self._cds_start[index]),
                                int(self._cds_end[index]))

    def slice_region(self, index: int, label: str = "CDS", trim: int = None) -> tuple:
        """(start, end) of one region, trimmed by `trim` nt at each end (`None` = the file's
        own `paper_cds_trim`). Collapses to (start, start) rather than inverting."""
        regions = self.regions_of(index)
        if label not in regions:
            raise SchemaError(
                "transcript %s has no %s region (present: %s)"
                % (self._ids[index], label, ", ".join(sorted(regions)) or "none"))
        start, end = regions[label]
        trim = self.trim if trim is None else trim
        if trim:
            start, end = start + trim, end - trim
        return (start, end) if end > start else (start, start)

    # ── whole-file arrays ─────────────────────────────────────────────────────
    def signal(self, name: str) -> np.ndarray:
        """One full-coordinate array."""
        if name not in SIGNALS:
            raise SchemaError("unknown signal %r" % name)
        return self.handle["coverage"][name][:]

    def cds_window_sums(self, values, trim: int = 0) -> np.ndarray:
        """Per-transcript sum of a full-coordinate array over the CDS trimmed by `trim`."""
        return window_sums(values, self._offset, self._cds_start + trim,
                           self._cds_end - trim)

    def close(self):
        self.handle.close()

    def __enter__(self):
        return self

    def __exit__(self, exc_type, exc, tb):
        self.close()
        return False


def open_coverage(path) -> CoverageFile:
    return CoverageFile(path)


def main(argv=None) -> int:
    parser = argparse.ArgumentParser(
        description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument("--validate", type=Path, required=True,
                        help="a .shared_coverage.h5 to validate")
    args = parser.parse_args(argv)
    problems = validate_file(args.validate)
    if problems:
        print("INVALID: %s" % args.validate)
        for problem in problems:
            print("  - %s" % problem)
        return 1
    with open_coverage(args.validate) as coverage:
        print("VALID: %s" % args.validate)
        print("  schema      %s" % SCHEMA)
        print("  sample      %s" % coverage.sample)
        print("  transcripts %d" % coverage.n_transcripts)
        print("  positions   %d" % coverage.n_positions)
        print("  trim        %d" % coverage.trim)
    return 0


if __name__ == "__main__":
    sys.path.insert(0, str(Path(__file__).resolve().parent))
    sys.exit(main())
