#!/usr/bin/env python3
"""The shared-coverage HDF5 container: constants, writer, reader, validator."""
from __future__ import annotations

import argparse
import datetime as _datetime
import json
import os
import sys
from pathlib import Path

import numpy as np

#: 1 -> 2 (2026-08-01): CIGAR-aware P-site placement became the only rule, and the file
SCHEMA = "riboflow_paper/shared-coverage/2"
SCHEMA_VERSION = 2
SUPPORTED_SCHEMA_VERSIONS = (2,)
GENERATOR = "riboflow_paper/code/coverage/build_shared_coverage.py"

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
EXON_SOURCE = "gencode_exon_features"
STOP_CODON_ASSIGNMENT = "utr3"

COVERAGE_DTYPE = np.int32
COVERAGE_MAX = int(np.iinfo(COVERAGE_DTYPE).max)

DEFAULT_CHUNK = 1 << 16
DEFAULT_GZIP_LEVEL = 9
DEFAULT_SHUFFLE = True

TRANSCRIPT_STR_COLUMNS = (
    "transcript_id", "gene_id", "transcript_name", "gene_name", "chrom", "strand")
TRANSCRIPT_INT_COLUMNS = (
    "transcript_len", "n_exons", "cds_len_gtf", "n_cds_exons",
    "coverage_offset", "exon_offset", "region_offset", "bin_offset")
TRANSCRIPT_COUNT_COLUMNS = (
    "n_genome_psite_events", "n_txome_psite_events",
    "n_genome_footprint_bases", "n_txome_footprint_bases")
TRANSCRIPT_BOOL_COLUMNS = (
    "in_transcriptome_reference", "length_filtered", "has_annotated_stop",
    "cds_divisible_by_3",
    "hist_cds_genome_psite_key", "hist_cds_txome_psite_key",
    "hist_cds_genome_footprint_key", "hist_cds_txome_footprint_key")

EXON_STR_COLUMNS = ("chrom",)
EXON_INT_COLUMNS = ("transcript_index", "exon_index", "g_start", "g_end", "tx_start", "tx_end")

REGION_STR_COLUMNS = ("label", "source")
REGION_INT_COLUMNS = (
    "transcript_index", "raw_header_start_1based", "raw_header_end_1based",
    "raw_bed_start", "raw_bed_end", "start", "end")

BIN_STR_COLUMNS = ("label", "ribopy_alias")
BIN_INT_COLUMNS = ("transcript_index", "start", "end")

class SchemaError(RuntimeError):
    """Raised when a coverage file violates the schema or an invariant."""

def _utc_now() -> str:
    return _datetime.datetime.utcnow().strftime("%Y-%m-%dT%H:%M:%SZ")

def _text(value) -> str:
    """HDF5 hands back `bytes` or `str` depending on how an attribute was written."""
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
    """{module: sha256} for the pipeline's own sources, plus a combined digest.

    Distinguishes "built by different code" from "built from different inputs", which
    input digests alone cannot do.
    """
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
    """The command as invoked, with `argv[0]` shortened to a repository-relative name.

    When `record_paths` is false, any argument that LOOKS like a path is reduced to its
    basename. An absolute path names the machine that built the file, the layout of a
    private project and often the operator; the input digests are what identify the data.

    "Looks like a path" is deliberately syntactic -- it contains a separator -- rather
    than "names an existing file". An output path given on the command line does not exist
    yet when the command runs, and an existence test would let exactly those through.
    """
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

def _vlen_str():
    import h5py
    return h5py.special_dtype(vlen=str)

STATE_NO_READS = "no_reads_assigned"
STATE_OUTSIDE_SLICE = "reads_outside_requested_slice"
STATE_ALL_ZERO_SLICE = "slice_all_zero"
STATE_COVERED = "covered"

def describe_coverage_state(n_events: int, slice_sum: int) -> str:
    """Which of the four coverage states holds for one transcript and one requested slice."""
    if n_events == 0:
        return STATE_NO_READS
    if slice_sum == 0:
        return STATE_OUTSIDE_SLICE
    return STATE_COVERED

class CoverageWriter:
    """Incremental, bounded, atomic writer.

    Usage:

        with CoverageWriter(path, sample=..., tables=..., provenance=...) as writer:
            writer.write_signal("genome_psite", genome_psite_array)
            writer.write_signal("txome_psite", txome_psite_array)
            ...            # free those arrays, then build the footprint pair
            writer.set_transcript_counts(...)

    `write_signal` may be called in any order, but calling it one PAIR at a time is what
    bounds peak memory to two full-coordinate arrays rather than four.
    """

    def __init__(self, path, sample, transcripts, exons, regions, ribo_bins,
                 offsets, provenance, paper_cds_trim=15, reference_name="",
                 chunk=DEFAULT_CHUNK, gzip_level=DEFAULT_GZIP_LEVEL,
                 shuffle=DEFAULT_SHUFFLE, ribo_bin_attrs=None, assay="ribo"):
        import h5py

        if assay not in ASSAYS:
            raise SchemaError("unknown assay %r; expected one of %s"
                              % (assay, ", ".join(ASSAYS)))
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

        self._write_attrs(paper_cds_trim, reference_name)
        self._write_transcripts()
        self._write_table("exons", exons, EXON_STR_COLUMNS, EXON_INT_COLUMNS)
        self._write_table("regions", regions, REGION_STR_COLUMNS, REGION_INT_COLUMNS)
        self._write_table("ribo_region_bins", ribo_bins, BIN_STR_COLUMNS, BIN_INT_COLUMNS,
                          attrs=ribo_bin_attrs)
        self._write_offsets(offsets)
        self._write_provenance(provenance)
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

    def _write_attrs(self, paper_cds_trim, reference_name):
        import psite_placement

        attrs = self.handle.attrs
        attrs["schema"] = SCHEMA
        attrs["schema_version"] = SCHEMA_VERSION
        attrs["coordinate_system"] = COORDINATE_SYSTEM
        attrs["exon_source"] = EXON_SOURCE
        attrs["stop_codon_assignment"] = STOP_CODON_ASSIGNMENT
        attrs["sample"] = self.sample
        attrs["assay"] = self.assay
        attrs["routes"] = list(ROUTES)
        attrs["psite_placement"] = psite_placement.PSITE_PLACEMENT
        attrs["reference_name"] = reference_name
        attrs["n_transcripts"] = int(len(self.transcripts))
        attrs["n_positions"] = self.n_positions
        attrs["paper_cds_trim"] = int(paper_cds_trim)
        attrs["created_utc"] = _utc_now()
        attrs["generator"] = GENERATOR

    _TABLE_COMPRESSION = dict(compression="gzip", compression_opts=1)

    def _dataset(self, group, name, data, dtype=None):
        """Create a chunked, lightly compressed dataset. Scalar-length tables are exempt:
        HDF5 requires chunking for compression, and a zero-length dataset cannot be chunked."""
        kwargs = dict(self._TABLE_COMPRESSION) if len(data) else {}
        if kwargs:
            kwargs["chunks"] = (min(len(data), 1 << 16),)
        if dtype is not None:
            kwargs["dtype"] = dtype
        return group.create_dataset(name, data=data, **kwargs)

    def _write_transcripts(self):
        group = self.handle.create_group("transcripts")
        frame = self.transcripts
        for column in TRANSCRIPT_STR_COLUMNS:
            self._dataset(group, column, frame[column].astype(str).to_numpy(),
                          dtype=_vlen_str())
        for column in TRANSCRIPT_INT_COLUMNS:
            self._dataset(group, column, frame[column].to_numpy(dtype=np.int64))
        for column in TRANSCRIPT_COUNT_COLUMNS:
            group.create_dataset(column, shape=(len(frame),), dtype=np.int64,
                                 fillvalue=-1)
        for column in TRANSCRIPT_BOOL_COLUMNS:
            if column in frame.columns:
                self._dataset(group, column, frame[column].to_numpy(dtype=bool))
            else:
                group.create_dataset(column, shape=(len(frame),), dtype=bool,
                                     fillvalue=False)

    def _write_table(self, name, frame, str_columns, int_columns, attrs=None):
        group = self.handle.create_group(name)
        for column in str_columns:
            self._dataset(group, column, frame[column].astype(str).to_numpy(),
                          dtype=_vlen_str())
        for column in int_columns:
            self._dataset(group, column, frame[column].to_numpy(dtype=np.int64))
        for key, value in (attrs or {}).items():
            group.attrs[key] = value

    def _write_offsets(self, offsets):
        group = self.handle.create_group("offsets")
        for route in ("genome", "transcriptome"):
            mapping = offsets.get(route, {})
            lengths = np.array(sorted(mapping), dtype=np.int32)
            values = np.array([mapping[int(k)] for k in lengths], dtype=np.int32)
            sub = group.create_group(route)
            sub.create_dataset("read_length", data=lengths)
            sub.create_dataset("psite_offset", data=values)

    def _write_provenance(self, provenance):
        provenance = dict(provenance)
        provenance.setdefault("generation", invocation())
        provenance.setdefault("code_version", code_version())

        group = self.handle.create_group("provenance")
        group.attrs["json"] = json.dumps(provenance, indent=2, sort_keys=True)
        for key in ("paper_cds_trim", "genome_uniqueness", "txome_uniqueness",
                    "left_span", "right_span", "psite_placement",
                    "appris_principal_ranks_consumed"):
            if key in provenance.get("parameters", {}):
                group.attrs[key] = provenance["parameters"][key]
        group.attrs["command"] = provenance["generation"].get("command", "")
        group.attrs["paths_redacted"] = bool(
            provenance["generation"].get("paths_redacted", True))
        group.attrs["code_version"] = provenance["code_version"].get("combined_sha256", "")

    def _create_coverage(self, chunk, gzip_level, shuffle):
        group = self.handle.create_group("coverage")
        chunk = int(min(chunk, max(self.n_positions, 1)))
        for name in SIGNALS:
            dataset = group.create_dataset(
                name, shape=(self.n_positions,), dtype=COVERAGE_DTYPE,
                chunks=(chunk,), compression="gzip", compression_opts=gzip_level,
                shuffle=shuffle, fillvalue=0)
            dataset.attrs["route"] = SIGNAL_ROUTE[name]
            dataset.attrs["measure"] = SIGNAL_MEASURE[name]
            dataset.attrs["assay"] = self.assay
        group.attrs["chunk"] = chunk
        group.attrs["gzip_level"] = gzip_level
        group.attrs["shuffle"] = bool(shuffle)
        group.attrs["routes"] = list(ROUTES)

    def write_signal(self, name, values, chunk_positions=1 << 22):
        """Write one full-coordinate array, in chunks, with an overflow check.

        `values` may be any integer dtype; it is range-checked and cast to int32. Passing
        int64 accumulators is normal and is exactly why the check exists.
        """
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

    def set_transcript_counts(self, column, values):
        """Fill one of the per-transcript event-count or coverage-state datasets."""
        if column in TRANSCRIPT_COUNT_COLUMNS:
            data = np.asarray(values, dtype=np.int64)
        elif column in TRANSCRIPT_BOOL_COLUMNS:
            data = np.asarray(values, dtype=bool)
        else:
            raise SchemaError("unknown per-transcript column %r" % column)
        if data.shape != (len(self.transcripts),):
            raise SchemaError("%s has %d values, expected %d"
                              % (column, data.size, len(self.transcripts)))
        self.handle["transcripts"][column][...] = data

    def finalize(self):
        """Validate the temporary file, then atomically move it into place."""
        missing = [s for s in SIGNALS if s not in self._written]
        if missing:
            raise SchemaError("signal(s) never written: %s" % ", ".join(missing))
        unfilled = [c for c in TRANSCRIPT_COUNT_COLUMNS
                    if int(self.handle["transcripts"][c][0]) == -1
                    and len(self.transcripts)]
        if unfilled:
            raise SchemaError(
                "per-transcript event counts never set: %s. They distinguish 'no reads "
                "assigned' from 'reads assigned outside the slice you asked for', so an "
                "unset value is not a safe default." % ", ".join(unfilled))
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

def validate_file(path) -> list:
    """Structural validation. Returns a list of problems; empty means valid.

    Deliberately checks invariants rather than just presence, so a truncated or
    partially-written file is caught structurally and not only by digest.
    """
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
        schema = handle.attrs.get("schema")
        if schema != SCHEMA:
            problems.append(
                "schema is %r, expected %r. A schema-1 file's genome P-sites were placed "
                "with a reference-offset rule, so it is not readable as a "
                "schema-2 file; rebuild it with code/coverage/build_shared_coverage.py."
                % (schema, SCHEMA))
            return problems

        for group in ("transcripts", "coverage", "exons", "regions", "offsets", "provenance"):
            if group not in handle:
                problems.append("missing group /%s" % group)
        if problems:
            return problems

        for attr in ("sample", "assay", "coordinate_system", "psite_placement",
                     "schema_version", "reference_name"):
            if attr not in handle.attrs:
                problems.append("missing root attribute %r" % attr)
        if handle.attrs.get("coordinate_system") not in (None, COORDINATE_SYSTEM):
            problems.append("coordinate_system is %r, expected %r"
                            % (handle.attrs.get("coordinate_system"), COORDINATE_SYSTEM))
        if handle.attrs.get("assay") is not None and handle.attrs["assay"] not in ASSAYS:
            problems.append("assay is %r, expected one of %s"
                            % (handle.attrs["assay"], ", ".join(ASSAYS)))
        provenance_attrs = handle["provenance"].attrs
        for attr in ("json", "command", "code_version"):
            if attr not in provenance_attrs:
                problems.append("missing /provenance attribute %r" % attr)

        n_transcripts = int(handle.attrs["n_transcripts"])
        n_positions = int(handle.attrs["n_positions"])
        transcripts = handle["transcripts"]

        for column in (TRANSCRIPT_STR_COLUMNS + TRANSCRIPT_INT_COLUMNS
                       + TRANSCRIPT_COUNT_COLUMNS + TRANSCRIPT_BOOL_COLUMNS):
            if column not in transcripts:
                problems.append("missing /transcripts/%s" % column)
            elif transcripts[column].shape != (n_transcripts,):
                problems.append("/transcripts/%s has shape %s, expected (%d,)"
                                % (column, transcripts[column].shape, n_transcripts))
        if problems:
            return problems

        transcript_len = transcripts["transcript_len"][:]
        coverage_offset = transcripts["coverage_offset"][:]

        if int(transcript_len.sum()) != n_positions:
            problems.append("sum(transcript_len) = %d but n_positions = %d"
                            % (transcript_len.sum(), n_positions))
        expected = np.concatenate([[0], np.cumsum(transcript_len)[:-1]]) \
            if n_transcripts else np.zeros(0, dtype=np.int64)
        if not np.array_equal(coverage_offset, expected):
            problems.append("coverage_offset is not the running sum of transcript_len")

        ids = [s.decode() if isinstance(s, bytes) else s
               for s in transcripts["transcript_id"][:]]
        if ids != sorted(ids):
            problems.append("transcripts are not in sorted transcript_id order "
                            "(storage order is load-bearing for the pooled-Pearson "
                            "reconstruction)")
        if len(set(ids)) != len(ids):
            problems.append("duplicate transcript_id")

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
            if dataset.attrs.get("route") != SIGNAL_ROUTE[name]:
                problems.append("/coverage/%s declares route %r, expected %r"
                                % (name, dataset.attrs.get("route"), SIGNAL_ROUTE[name]))
            if dataset.attrs.get("measure") != SIGNAL_MEASURE[name]:
                problems.append("/coverage/%s declares measure %r, expected %r"
                                % (name, dataset.attrs.get("measure"), SIGNAL_MEASURE[name]))

        exons = handle["exons"]
        exon_tx = exons["transcript_index"][:]
        spliced = np.bincount(exon_tx, weights=(exons["g_end"][:] - exons["g_start"][:]),
                              minlength=n_transcripts).astype(np.int64)
        if not np.array_equal(spliced, transcript_len):
            bad = int(np.argmax(spliced != transcript_len))
            problems.append(
                "spliced exon length != transcript_len, first at transcript %d (%s): %d vs %d"
                % (bad, ids[bad], spliced[bad], transcript_len[bad]))

        regions = handle["regions"]
        region_tx = regions["transcript_index"][:]
        covered = np.bincount(region_tx,
                              weights=(regions["end"][:] - regions["start"][:]),
                              minlength=n_transcripts).astype(np.int64)
        if not np.array_equal(covered, transcript_len):
            bad = int(np.argmax(covered != transcript_len))
            problems.append(
                "regions do not tile [0, transcript_len), first at transcript %d (%s): "
                "%d vs %d" % (bad, ids[bad], covered[bad], transcript_len[bad]))

        if "json" not in handle["provenance"].attrs:
            problems.append("/provenance has no `json` attribute")
        else:
            try:
                json.loads(handle["provenance"].attrs["json"])
            except ValueError as exc:
                problems.append("/provenance/json is not valid JSON: %s" % exc)

    return problems

class CoverageFile:
    """Read access with the schema checked on open, not on first surprising failure."""

    def __init__(self, path):
        import h5py

        self.path = Path(path)
        problems = validate_file(self.path)
        if problems:
            raise SchemaError("%s failed validation:\n%s"
                              % (self.path, "\n".join("  - %s" % p for p in problems)))
        self.handle = h5py.File(self.path, "r")
        attrs = self.handle.attrs
        self.sample = attrs["sample"]
        self.assay = _text(attrs.get("assay", ""))
        self.routes = tuple(_text(r) for r in attrs.get("routes", ROUTES))
        self.coordinate_system = _text(attrs.get("coordinate_system", ""))
        self.psite_placement = _text(attrs.get("psite_placement", ""))
        self.schema_version = int(attrs.get("schema_version", 0))
        self.reference_name = _text(attrs.get("reference_name", ""))
        self.trim = int(attrs["paper_cds_trim"])
        self.n_transcripts = int(attrs["n_transcripts"])
        self._ids = [s.decode() if isinstance(s, bytes) else s
                     for s in self.handle["transcripts"]["transcript_id"][:]]
        self._gene_ids = [s.decode() if isinstance(s, bytes) else s
                          for s in self.handle["transcripts"]["gene_id"][:]]

    @property
    def provenance(self) -> dict:
            return json.loads(self.handle["provenance"].attrs["json"])

    def identity(self) -> dict:
        """Everything a consumer needs to decide whether this file is the right one."""
        return {
            "path": str(self.path),
            "schema": _text(self.handle.attrs.get("schema", "")),
            "schema_version": self.schema_version,
            "sample": _text(self.sample),
            "assay": self.assay,
            "routes": list(self.routes),
            "coordinate_system": self.coordinate_system,
            "psite_placement": self.psite_placement,
            "reference_name": self.reference_name,
            "paper_cds_trim": self.trim,
            "n_transcripts": self.n_transcripts,
            "has_ribo_region_bins": self.has_ribo_region_bins,
            "created_utc": _text(self.handle.attrs.get("created_utc", "")),
        }

    @property
    def has_ribo_region_bins(self) -> bool:
        group = self.handle.get("ribo_region_bins")
        return bool(group is not None and len(group.get("start", ())))

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

        lengths = self.handle["transcripts"]["transcript_len"]
        cds = self.handle["transcripts"]["cds_len_gtf"]
        listing = "\n".join(
            "    %s  transcript_len=%d  cds_len_gtf=%d"
            % (self._ids[i], int(lengths[i]), int(cds[i])) for i in hits)
        raise SchemaError(
            "gene %r maps to %d eligible transcripts and no unique one resolves it. "
            "Pass --transcript-id to choose; this is not guessed.\n%s"
            % (gene_id, len(hits), listing))

    def transcript_info(self, index: int) -> dict:
        group = self.handle["transcripts"]
        info = {}
        for column in TRANSCRIPT_STR_COLUMNS:
            value = group[column][index]
            info[column] = value.decode() if isinstance(value, bytes) else value
        for column in TRANSCRIPT_INT_COLUMNS + TRANSCRIPT_COUNT_COLUMNS:
            info[column] = int(group[column][index])
        for column in TRANSCRIPT_BOOL_COLUMNS:
            info[column] = bool(group[column][index])
        return info

    def get_track(self, index: int, signal: str) -> np.ndarray:
        if signal not in SIGNALS:
            raise SchemaError("unknown signal %r; expected one of %s"
                              % (signal, ", ".join(SIGNALS)))
        group = self.handle["transcripts"]
        start = int(group["coverage_offset"][index])
        length = int(group["transcript_len"][index])
        return self.handle["coverage"][signal][start:start + length]

    def get_tracks(self, index: int, signals=SIGNALS) -> dict:
        return {name: self.get_track(index, name) for name in signals}

    def regions_of(self, index: int) -> dict:
        """{label: (start, end)} normalized, for one transcript."""
        group = self.handle["transcripts"]
        start = int(group["region_offset"][index])
        region_index = self.handle["regions"]["transcript_index"]
        stop = start
        while stop < region_index.shape[0] and int(region_index[stop]) == index:
            stop += 1
        labels = self.handle["regions"]["label"][start:stop]
        starts = self.handle["regions"]["start"][start:stop]
        ends = self.handle["regions"]["end"][start:stop]
        return {(l.decode() if isinstance(l, bytes) else l): (int(s), int(e))
                for l, s, e in zip(labels, starts, ends)}

    def ribo_bins_of(self, index: int) -> list:
        """The five-way ribopy overlay for one transcript, 5'->3'.

        Returns `[(ribopy_alias, label, start, end), …]` -- for example
        `[("UTR5", "UTR5_OUTER", 0, 116), ("UTR5J", "START_WINDOW", 116, 162), …]`.

        The alias is ribopy's `UTR5 / UTR5J / CDS / UTR3J / UTR3` naming and the
        label is what the bin actually is. They are kept apart on purpose: a bare `CDS`
        always means the canonical region from `/regions`, never ribopy's shrunken one.

        Empty for a transcript whose bins were all clipped away, and for any file written
        without them.
        """
        group = self.handle.get("ribo_region_bins")
        if group is None or "start" not in group or not len(group["start"]):
            return []
        start = int(self.handle["transcripts"]["bin_offset"][index])
        bin_index = group["transcript_index"]
        stop = start
        while stop < bin_index.shape[0] and int(bin_index[stop]) == index:
            stop += 1
        rows = zip(group["ribopy_alias"][start:stop], group["label"][start:stop],
                   group["start"][start:stop], group["end"][start:stop])
        return [(_text(alias), _text(label), int(lo), int(hi))
                for alias, label, lo, hi in rows]

    def slice_region(self, index: int, label: str = "CDS", trim: int = None) -> tuple:
        """(start, end) of one region, optionally trimmed by `trim` nt at each end.

        `trim=None` uses the file's own `paper_cds_trim`, so the published Figure-3 framing
        is reproduced without the caller hardcoding 15.
        """
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
        print("  positions   %d" % int(coverage.handle.attrs["n_positions"]))
        print("  trim        %d" % coverage.trim)
    return 0

if __name__ == "__main__":
    sys.path.insert(0, str(Path(__file__).resolve().parent))
    sys.exit(main())
