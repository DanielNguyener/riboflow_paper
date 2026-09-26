#!/usr/bin/env python3
"""A sparse per-read alignment-state store: the two BAM passes, done once, kept on disk.

`gene_read_partition_lib.compute_partition` reads BOTH BAMs end to end on every call --
~136 million pysam records -- before it looks at a single gene, because it was written for
a handful of hand-picked genes. Asking it about every gene means paying that per batch.

This stores exactly the state its classifier consumes, so the passes happen once per
sample and every later question is a slice:

    python code/clustering/read_state.py --sample HeLa         # ~20 min, once
    python code/clustering/build_gene_counts.py --sample HeLa   # minutes, any gene set

What is stored is what `classify_union` asks for and nothing else: for every reported
genome alignment its locus, score, NH and CIGAR blocks; for every read its primary, its
uniqueness, and the transcript its transcriptome primary landed on when it has one
(presence only: the post-dedup BAM is already RiboFlow_v2's MAPQ >= 10 set). Read IDs are interned to int32
and the strings are then dropped -- the classifier treats a read id as an opaque key
(`tie_biotype_lib.LOCUS_FRAME_COLUMNS` says so explicitly), so integers serve.

Sparse throughout, following `genome_coverage/`: nothing absent is stored. A read with no
transcriptome primary has no row in `txome/` rather than a row of sentinels; positions are
delta-encoded within a chromosome and every dataset is gzip-9 + shuffle.
"""
from __future__ import annotations

import argparse
import array
import datetime
import json
import os
import sys
import time
from pathlib import Path

import numpy as np

HERE = Path(__file__).resolve().parent
REPO = HERE.parents[1]

SCHEMA = "riboflow_paper/read_state/3"

#: Same sentinel `gene_read_partition_lib` uses: far below any real alignment score, so an
#: unscored alignment can never tie with a scored one in the pseudogene-tie test.
MISSING_AS = -(10 ** 9)

FLAG_SECONDARY = 1
FLAG_REVERSE = 2
FLAG_PRESENT = 1
FLAG_UNIQUE = 2

#: gzip-4, not -9: these are hundreds of millions of delta-encoded integers, and -9
#: buys a few per cent for several times the write time.
_COMPRESSION = dict(compression="gzip", compression_opts=4, shuffle=True)


def _delta(values):
    """First value, then gaps. Sorted positions compress far better this way."""
    out = np.asarray(values, dtype=np.int64).copy()
    if len(out) > 1:
        out[1:] -= out[:-1]
    return out


def _undelta(values):
    return np.cumsum(np.asarray(values, dtype=np.int64))


# ── building ─────────────────────────────────────────────────────────────────

class ReadKeys:
    """Read name -> int64, without keeping 99 million Python strings alive.

    Illumina/SRA names in this cohort are `<run>.<serial>` (`SRR3306589.144250553`), so a
    read is identified exactly by a small run code and an integer serial, packed into one
    int64. That turns interning into `np.unique` + `searchsorted` over an array instead of
    a dict of strings -- on the pre-dedup BAM (237 million alignments, ~99 million reads)
    the dict alone would want well over ten gigabytes.

    A name that does not parse falls the whole run back to a plain dictionary, which is
    correct for any naming scheme and merely heavier. There is no third behaviour: the
    encoding is either exact for every read or not used at all.
    """

    SERIAL_BITS = 40
    MAX_SERIAL = (1 << SERIAL_BITS) - 1

    def __init__(self):
        self.runs = {}
        self.run_names = []
        self.fallback = None          # name -> int, only if parsing ever fails

    def encode(self, name):
        if self.fallback is not None:
            return self._encode_fallback(name)
        run, _, serial = name.rpartition(".")
        if run and serial.isdigit():
            value = int(serial)
            if value <= self.MAX_SERIAL:
                code = self.runs.get(run)
                if code is None:
                    code = len(self.run_names)
                    self.runs[run] = code
                    self.run_names.append(run)
                return (code << self.SERIAL_BITS) | value
        self._start_fallback()
        return self._encode_fallback(name)

    def _start_fallback(self):
        self.fallback = {}

    def _encode_fallback(self, name):
        value = self.fallback.get(name)
        if value is None:
            value = len(self.fallback)
            self.fallback[name] = value
        return value

    @property
    def parsed(self):
        return self.fallback is None


def scan_genome(genome_bam, keys, log):
    """One pass. Per-alignment columns, plus the CIGAR blocks of every primary.

    Nothing is keyed by read id here: read indices do not exist until both BAMs have been
    seen, so alignments carry their packed name key and are resolved afterwards.
    """
    import pysam

    al_key = array.array("q")
    al_chrom = array.array("h")
    al_start = array.array("i")
    al_end = array.array("i")
    al_pos5 = array.array("i")
    al_score = array.array("i")
    al_nh = array.array("h")
    al_flags = array.array("B")
    block_start = array.array("i")
    block_end = array.array("i")
    block_count = array.array("h")          # per primary, in scan order

    chrom_names, chrom_codes = [], {}
    started, n_records = time.time(), 0
    bam = pysam.AlignmentFile(str(genome_bam), "rb")
    try:
        for read in bam.fetch(until_eof=True):
            if read.is_unmapped or read.is_supplementary:
                continue
            n_records += 1
            if n_records % 25_000_000 == 0:
                log("  %d million alignments, %.0f s" % (n_records // 1_000_000,
                                                         time.time() - started))
            try:
                nh = int(read.get_tag("NH"))
            except KeyError:
                raise SystemExit(
                    "%s has an alignment with no NH tag. This analysis is about "
                    "multimapping, so a BAM that cannot report NH cannot answer it."
                    % genome_bam)
            blocks = read.get_blocks()
            if not blocks:
                continue
            name = read.reference_name
            code = chrom_codes.get(name)
            if code is None:
                code = len(chrom_names)
                chrom_codes[name] = code
                chrom_names.append(name)

            reverse = bool(read.is_reverse)
            start = blocks[0][0]
            end = blocks[-1][1]
            secondary = bool(read.is_secondary)

            al_key.append(keys.encode(read.query_name))
            al_chrom.append(code)
            al_start.append(start)
            al_end.append(end)
            al_pos5.append((end - 1) if reverse else start)
            al_score.append(int(read.get_tag("AS")) if read.has_tag("AS") else MISSING_AS)
            al_nh.append(nh)
            al_flags.append((FLAG_SECONDARY if secondary else 0)
                            | (FLAG_REVERSE if reverse else 0))
            if not secondary:
                for begin, finish in blocks:
                    block_start.append(begin)
                    block_end.append(finish)
                block_count.append(len(blocks))
    finally:
        bam.close()

    log("genome: %d alignments, %.0f s" % (n_records, time.time() - started))
    columns = {"key": np.frombuffer(al_key, dtype=np.int64),
               "chrom": np.frombuffer(al_chrom, dtype=np.int16),
               "start": np.frombuffer(al_start, dtype=np.int32),
               "end": np.frombuffer(al_end, dtype=np.int32),
               "pos5": np.frombuffer(al_pos5, dtype=np.int32),
               "score": np.frombuffer(al_score, dtype=np.int32),
               "nh": np.frombuffer(al_nh, dtype=np.int16),
               "flags": np.frombuffer(al_flags, dtype=np.uint8)}
    blocks = {"start": np.frombuffer(block_start, dtype=np.int32),
              "end": np.frombuffer(block_end, dtype=np.int32),
              "count": np.frombuffer(block_count, dtype=np.int16).astype(np.int32)}
    return columns, blocks, chrom_names


def scan_txome(txome_bam, base2ver, keys, log):
    """One pass over the transcriptome primaries. Sparse: only APPRIS-resolving reads.

    Mirrors `reference_lib.read_txome_primary` exactly -- every primary, the reference
    name resolved through `base2ver`.
    """
    import pysam

    tx_key = array.array("q")
    tx_tid = array.array("i")

    tid_codes, tid_names = {}, []
    started, n_primary = time.time(), 0
    bam = pysam.AlignmentFile(str(txome_bam), "rb")
    try:
        for read in bam.fetch(until_eof=True):
            if read.is_unmapped or read.is_secondary or read.is_supplementary:
                continue
            n_primary += 1
            base = read.reference_name.split(".", 1)[0].split("|", 1)[0]
            tid = base2ver.get(base)
            if tid is None:
                continue
            code = tid_codes.get(tid)
            if code is None:
                code = len(tid_names)
                tid_codes[tid] = code
                tid_names.append(tid)
            tx_key.append(keys.encode(read.query_name))
            tx_tid.append(code)
    finally:
        bam.close()

    log("transcriptome: %d primaries, %d resolving to an APPRIS transcript, %.0f s"
        % (n_primary, len(tx_key), time.time() - started))
    return {"key": np.frombuffer(tx_key, dtype=np.int64),
            "tid": np.frombuffer(tx_tid, dtype=np.int32)}, tid_names


def _regroup_blocks(counts, order):
    """Row indices that reorder a ragged block list by `order`, vectorised."""
    offsets = np.zeros(len(counts) + 1, dtype=np.int64)
    np.cumsum(counts, out=offsets[1:])
    sizes = counts[order]
    starts = offsets[order]
    total = int(sizes.sum())
    within = np.arange(total, dtype=np.int64)
    group_start = np.zeros(len(sizes) + 1, dtype=np.int64)
    np.cumsum(sizes, out=group_start[1:])
    within -= np.repeat(group_start[:-1], sizes)
    return np.repeat(starts, sizes) + within, sizes


def write_state(path, sample, genome_bam, txome_bam, columns, blocks,
                chrom_names, txome, tid_names, keys, log):
    """Resolve read indices, sort, and write. Everything sparse and delta-encoded."""
    import h5py

    started = time.time()
    # Read indices exist only once both BAMs have been seen: a read may appear on the
    # transcriptome route and nowhere in the genome BAM, and it still needs an index.
    unique_keys = np.unique(np.concatenate([columns["key"], txome["key"]]))
    n_reads = int(len(unique_keys))
    al_read = np.searchsorted(unique_keys, columns["key"]).astype(np.int32)
    tx_read = np.searchsorted(unique_keys, txome["key"]).astype(np.int32)
    log("resolved %d distinct read ids (%s)"
        % (n_reads, "packed name keys" if keys.parsed else "dictionary fallback"))

    primary_rows = np.flatnonzero((columns["flags"] & FLAG_SECONDARY) == 0)
    primary_reads = al_read[primary_rows]
    primary_order = np.argsort(primary_reads, kind="stable")
    primary_reads = primary_reads[primary_order]
    if len(primary_reads) and (np.diff(primary_reads) == 0).any():
        raise SystemExit("a read has two primary genome alignments; the BAM is malformed "
                         "or the read-name key is colliding")
    primary_rows = primary_rows[primary_order]
    primary_flags = (np.uint8(FLAG_PRESENT)
                     | np.where(columns["nh"][primary_rows] == 1,
                                np.uint8(FLAG_UNIQUE), np.uint8(0)).astype(np.uint8))
    gather, sizes = _regroup_blocks(blocks["count"], primary_order)
    flat_start = blocks["start"][gather]
    flat_end = blocks["end"][gather]
    block_offsets = np.zeros(len(sizes) + 1, dtype=np.int64)
    np.cumsum(sizes, out=block_offsets[1:])

    # Alignments sorted by (chromosome, start): a gene's reads are then a binary-search
    # slice instead of a scan.
    order = np.lexsort((columns["start"], columns["chrom"]))
    al_read = al_read[order]
    for key in ("chrom", "start", "end", "pos5", "score", "nh", "flags"):
        columns[key] = columns[key][order]
    new_row = np.empty(len(order), dtype=np.int64)
    new_row[order] = np.arange(len(order), dtype=np.int64)
    primary_align = new_row[primary_rows]
    del new_row, order

    n_chrom = len(chrom_names)
    chrom_offsets = np.searchsorted(columns["chrom"], np.arange(n_chrom + 1), side="left")
    chrom_offsets[-1] = len(columns["chrom"])
    spans = columns["end"] - columns["start"]
    max_span = np.zeros(n_chrom, dtype=np.int32)
    for code in range(n_chrom):
        lo, hi = int(chrom_offsets[code]), int(chrom_offsets[code + 1])
        if hi > lo:
            max_span[code] = int(spans[lo:hi].max())
    del spans

    # Alignments grouped by read, so a multimapper's other loci are one slice.
    by_read_order = np.argsort(al_read, kind="stable").astype(np.int64)
    read_offsets = np.searchsorted(al_read[by_read_order], np.arange(n_reads + 1),
                                   side="left")
    read_offsets[-1] = len(by_read_order)

    txome_order = np.argsort(tx_read, kind="stable")
    tx_read = tx_read[txome_order]
    txome["tid"] = txome["tid"][txome_order]
    tid_order = np.argsort(txome["tid"], kind="stable").astype(np.int64)
    tid_offsets = np.searchsorted(txome["tid"][tid_order],
                                  np.arange(len(tid_names) + 1), side="left")
    tid_offsets[-1] = len(tid_order)

    path.parent.mkdir(parents=True, exist_ok=True)
    temporary = path.with_suffix(".tmp")
    with h5py.File(temporary, "w") as handle:
        handle.attrs["schema"] = SCHEMA
        handle.attrs["schema_version"] = 3
        handle.attrs["sample"] = sample
        handle.attrs["genome_bam"] = str(genome_bam)
        handle.attrs["transcriptome_bam"] = str(txome_bam)
        handle.attrs["genome_uniqueness"] = "NH==1 (primary)"
        handle.attrs["txome_presence"] = ("primary alignment in the post-dedup BAM "
                                          "(RiboFlow_v2 MAPQ>=10)")
        handle.attrs["n_reads"] = n_reads
        handle.attrs["n_alignments"] = int(len(al_read))
        handle.attrs["n_primary"] = int(len(primary_reads))
        handle.attrs["created_utc"] = datetime.datetime.utcnow().isoformat() + "Z"
        handle.attrs["provenance"] = json.dumps(
            {"generator": "code/clustering/read_state.py",
             "missing_as": MISSING_AS,
             "read_key": "packed" if keys.parsed else "dictionary",
             "runs": list(keys.run_names),
             "note": "read ids are interned to int32; the names are not stored"})

        def store(name, values, delta=False):
            handle.create_dataset(
                name, data=(_delta(values) if delta else np.asarray(values)),
                **_COMPRESSION)

        reference = handle.create_group("reference")
        reference.create_dataset("chrom_names", data=np.array(chrom_names, dtype="S"),
                                 **_COMPRESSION)
        reference.create_dataset("transcript_ids", data=np.array(tid_names, dtype="S"),
                                 **_COMPRESSION)
        store("reference/chrom_offsets", chrom_offsets)
        store("reference/chrom_max_span", max_span)

        store("alignments/read", al_read)
        store("alignments/chrom", columns["chrom"])
        store("alignments/start", columns["start"], delta=True)
        store("alignments/end", columns["end"])
        store("alignments/pos5", columns["pos5"])
        store("alignments/score", columns["score"])
        store("alignments/nh", columns["nh"])
        store("alignments/flags", columns["flags"])

        # The packed name key of every read, in read-index order. A read index is a rank in
        # THIS store's key array, so it means nothing in another store; the key does. Two
        # stores built from different BAMs of the same sample are joined on this.
        store("reads/key", unique_keys, delta=True)

        store("by_read/align_idx", by_read_order)
        store("by_read/offsets", read_offsets)

        store("primary/read", primary_reads, delta=True)
        store("primary/align_idx", primary_align)
        store("primary/flags", primary_flags)
        store("primary/block_offsets", block_offsets)
        store("blocks/start", flat_start)
        store("blocks/end", flat_end)

        store("txome/read", tx_read, delta=True)
        store("txome/tid", txome["tid"])
        store("txome/by_tid_idx", tid_order)
        store("txome/by_tid_offsets", tid_offsets)

    os.replace(temporary, path)
    log("wrote %s (%.0f MB, %.0f s)"
        % (path, path.stat().st_size / 1e6, time.time() - started))


def build(sample, genome_bam, txome_bam, base2ver, path, log):
    keys = ReadKeys()
    log("pass 1/2: genome BAM %s" % genome_bam)
    columns, blocks, chrom_names = scan_genome(genome_bam, keys, log)
    log("pass 2/2: transcriptome BAM %s" % txome_bam)
    txome, tid_names = scan_txome(txome_bam, base2ver, keys, log)
    write_state(Path(path), sample, genome_bam, txome_bam, columns, blocks,
                chrom_names, txome, tid_names, keys, log)
    return Path(path)


class _Membership:
    """`qname in genome_present` over a sorted int array, without a 13-million-entry set."""

    def __init__(self, sorted_values):
        self._values = sorted_values

    def __contains__(self, key):
        index = np.searchsorted(self._values, key)
        return index < len(self._values) and self._values[index] == key

    def __len__(self):
        return len(self._values)


class _TxomeMap:
    """`qname in txome_present` and `txome_present[qname]` -> transcript id."""

    def __init__(self, reads, tids, tid_names):
        self._reads = reads
        self._tids = tids
        self._names = tid_names

    def _row(self, key):
        index = np.searchsorted(self._reads, key)
        if index < len(self._reads) and self._reads[index] == key:
            return int(index)
        return None

    def __contains__(self, key):
        return self._row(key) is not None

    def __getitem__(self, key):
        row = self._row(key)
        if row is None:
            raise KeyError(key)
        return self._names[self._tids[row]]

    def get(self, key, default=None):
        try:
            return self[key]
        except KeyError:
            return default

    def __len__(self):
        return len(self._reads)


class ReadState:
    """Everything `classify_union` needs, sliced out of the HDF5 on demand."""

    def __init__(self, path):
        import h5py

        self.path = Path(path)
        self._handle = h5py.File(self.path, "r")
        if self._handle.attrs.get("schema") != SCHEMA:
            raise SystemExit("%s is not a %s file" % (path, SCHEMA))
        self.sample = str(self._handle.attrs["sample"])
        self.n_reads = int(self._handle.attrs["n_reads"])

        self.chrom_names = [n.decode() for n in self._handle["reference/chrom_names"][:]]
        self.chrom_code = {name: code for code, name in enumerate(self.chrom_names)}
        self.transcript_ids = [t.decode() for t in
                               self._handle["reference/transcript_ids"][:]]
        self.tid_code = {tid: code for code, tid in enumerate(self.transcript_ids)}
        self._chrom_offsets = self._handle["reference/chrom_offsets"][:]
        self._chrom_max_span = self._handle["reference/chrom_max_span"][:]

        self._al_read = self._handle["alignments/read"][:]
        self._al_start = _undelta(self._handle["alignments/start"][:]).astype(np.int64)
        self._al_end = self._handle["alignments/end"][:]
        self._al_pos5 = self._handle["alignments/pos5"][:]
        self._al_score = self._handle["alignments/score"][:]
        self._al_nh = self._handle["alignments/nh"][:]
        self._al_flags = self._handle["alignments/flags"][:]
        self._al_chrom = self._handle["alignments/chrom"][:]

        self.read_key = _undelta(self._handle["reads/key"][:])
        self._by_read = self._handle["by_read/align_idx"][:]
        self._read_offsets = self._handle["by_read/offsets"][:]

        self._p_read = _undelta(self._handle["primary/read"][:])
        self._p_align = self._handle["primary/align_idx"][:]
        self._p_flags = self._handle["primary/flags"][:]
        self._p_block_offsets = self._handle["primary/block_offsets"][:]
        self._b_start = self._handle["blocks/start"][:]
        self._b_end = self._handle["blocks/end"][:]

        self._tx_read = _undelta(self._handle["txome/read"][:])
        self._tx_tid = self._handle["txome/tid"][:]
        self._tx_by_tid = self._handle["txome/by_tid_idx"][:]
        self._tx_tid_offsets = self._handle["txome/by_tid_offsets"][:]

        # Which alignments can place their read at a gene: the primary, or a secondary
        # whose AS equals the primary's (both present). `gene_read_partition_lib.qualifies`.
        primary_score = np.full(self.n_reads, MISSING_AS, dtype=np.int64)
        has_primary = np.zeros(self.n_reads, dtype=bool)
        primary_score[self._p_read] = self._al_score[self._p_align]
        has_primary[self._p_read] = True
        read_score = primary_score[self._al_read]
        self._qualifies = (((self._al_flags & FLAG_SECONDARY) == 0)
                           | (has_primary[self._al_read] & (read_score != MISSING_AS)
                              & (self._al_score != MISSING_AS)
                              & (self._al_score == read_score)))
        del primary_score, has_primary, read_score

        self.genome_present = _Membership(self._p_read)
        self.genome_unique = _Membership(
            self._p_read[(self._p_flags & FLAG_UNIQUE) != 0])
        self.txome_present = _TxomeMap(self._tx_read, self._tx_tid, self.transcript_ids)

    def close(self):
        self._handle.close()

    def __enter__(self):
        return self

    def __exit__(self, *_exc):
        self.close()

    # -- the gene-assignment rule -------------------------------------------
    def gene_side(self, chrom, start, end):
        """Read ids with a QUALIFYING alignment overlapping [start, end) on `chrom`.

        Qualifying = the primary, or a secondary score-tied with it: the population
        `gene_read_partition_lib.resolve_genome_side` keeps from an indexed BAM fetch.
        Strand is ignored. Found by binary search, over alignments sorted by (chromosome,
        start); the lower bound backs off by the chromosome's longest alignment span (over
        all alignments, so a valid bound for any subset) so a spliced read starting before
        the gene and reaching into it is not missed.
        """
        code = self.chrom_code.get(str(chrom))
        if code is None:
            return np.empty(0, dtype=np.int32)
        lo, hi = int(self._chrom_offsets[code]), int(self._chrom_offsets[code + 1])
        if hi <= lo:
            return np.empty(0, dtype=np.int32)
        starts = self._al_start[lo:hi]
        first = lo + int(np.searchsorted(
            starts, int(start) - int(self._chrom_max_span[code]), side="left"))
        last = lo + int(np.searchsorted(starts, int(end), side="left"))
        if last <= first:
            return np.empty(0, dtype=np.int32)
        window = slice(first, last)
        keep = (self._al_end[window] > int(start)) & self._qualifies[window]
        return np.unique(self._al_read[window][keep])

    def transcript_side(self, tid):
        """Read ids whose transcriptome primary is `tid`."""
        code = self.tid_code.get(tid)
        if code is None:
            return np.empty(0, dtype=np.int32)
        lo = int(self._tx_tid_offsets[code])
        hi = int(self._tx_tid_offsets[code + 1])
        return self._tx_read[self._tx_by_tid[lo:hi]]

    # -- the two mappings the chain consumes, served lazily -----------------
    # `classify_union` only asks for the reads it actually needs: `primary.get(q)` for the
    # reads at the gene, `records[q]` for the multimappers among them. Views keep it that
    # way -- building either eagerly for a whole gene would rebuild block lists for every
    # uniquely-mapped read the tie test never looks at.

    @property
    def primary(self):
        return _PrimaryView(self)

    @property
    def records(self):
        return _RecordsView(self)

    # -- the two dictionaries the chain consumes ----------------------------
    def _primary_row(self, read_index):
        index = np.searchsorted(self._p_read, read_index)
        if index < len(self._p_read) and self._p_read[index] == read_index:
            return int(index)
        return None

    def primary_dict(self, read_indices):
        """read id -> (chromosome, strand, blocks, NH, AS), as `classify_union` expects."""
        out = {}
        for read_index in read_indices:
            row = self._primary_row(read_index)
            if row is None:
                continue
            align = int(self._p_align[row])
            begin = int(self._p_block_offsets[row])
            finish = int(self._p_block_offsets[row + 1])
            blocks = list(zip(self._b_start[begin:finish].tolist(),
                              self._b_end[begin:finish].tolist()))
            out[int(read_index)] = (
                self.chrom_names[int(self._al_chrom[align])],
                "-" if self._al_flags[align] & FLAG_REVERSE else "+",
                blocks, int(self._al_nh[align]), int(self._al_score[align]))
        return out

    def records_dict(self, read_indices):
        """read id -> [(chromosome, pos5, AS, is_secondary), ...] over EVERY locus."""
        out = {}
        for read_index in read_indices:
            lo = int(self._read_offsets[read_index])
            hi = int(self._read_offsets[read_index + 1])
            if hi <= lo:
                continue
            rows = self._by_read[lo:hi]
            out[int(read_index)] = [
                (self.chrom_names[int(self._al_chrom[row])], int(self._al_pos5[row]),
                 int(self._al_score[row]), bool(self._al_flags[row] & FLAG_SECONDARY))
                for row in rows]
        return out


class _PrimaryView:
    """`primary.get(read)` -> (chromosome, strand, blocks, NH, AS)."""

    def __init__(self, state):
        self._state = state

    def get(self, key, default=None):
        found = self._state.primary_dict([key])
        return found.get(int(key), default)

    def __getitem__(self, key):
        found = self.get(key)
        if found is None:
            raise KeyError(key)
        return found

    def __contains__(self, key):
        return self.get(key) is not None


class _RecordsView:
    """`records[read]` -> [(chromosome, pos5, AS, is_secondary), ...] over every locus."""

    def __init__(self, state):
        self._state = state

    def get(self, key, default=None):
        found = self._state.records_dict([key])
        return found.get(int(key), default)

    def __getitem__(self, key):
        found = self.get(key)
        if found is None:
            raise KeyError(key)
        return found

    def __contains__(self, key):
        return self.get(key) is not None


def default_path(output_root, sample, label=""):
    """One store per (sample, BAM set): post-dedup and pre-dedup coexist."""
    stem = "%s.%s" % (sample, label) if label else sample
    return Path(output_root) / ("%s.read_state.h5" % stem)


def main(argv=None):
    parser = argparse.ArgumentParser(
        description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument("--sample", required=True)
    parser.add_argument("--genome-bam")
    parser.add_argument("--transcriptome-bam")
    parser.add_argument("--bams")
    parser.add_argument("--gtf")
    parser.add_argument("--appris")
    parser.add_argument("--output", help="default: results/clustering")
    parser.add_argument("--label", default="post_dedup",
                        help="BAM-set tag; the store is <sample>.<label>.read_state.h5")
    parser.add_argument("--force", action="store_true", help="rebuild an existing file")
    args = parser.parse_args(argv)

    sys.path.insert(0, str(REPO / "code" / "common"))
    import inputs as input_resolver
    log = input_resolver.make_log("read_state")

    paths = input_resolver.resolve_external_inputs(
        args.bams, args.gtf, args.appris, sample=args.sample)
    os.environ["RIBOFLOW_PAPER_GTF"] = str(Path(paths["gtf"]).resolve())
    os.environ["RIBOFLOW_PAPER_APPRIS"] = str(Path(paths["appris"]).resolve())
    genome_bam = Path(args.genome_bam or paths["ribo_genome"])
    txome_bam = Path(args.transcriptome_bam or paths["ribo_txome"])

    output = Path(args.output) if args.output else REPO / "results" / "clustering"
    destination = default_path(output, args.sample, args.label)
    if destination.exists() and not args.force:
        raise SystemExit("%s exists; pass --force to rebuild" % destination)

    saved = list(sys.path)
    for entry in (REPO / "code" / "read_taxonomy", REPO / "code" / "common",
                  REPO / "code" / "common" / "ribo_seq_qc"):
        sys.path.insert(0, str(entry))
    import bam_inputs
    import reference_lib
    sys.path[:] = saved + [p for p in sys.path if p not in saved]

    log("resolving the APPRIS transcript table")
    base2ver = reference_lib.build_transcript_table()["base2ver"]
    build(args.sample, genome_bam, txome_bam, base2ver, destination, log)
    return 0


if __name__ == "__main__":
    sys.exit(main())
