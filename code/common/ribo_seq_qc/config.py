#!/usr/bin/env python3
"""Shared configuration for the multi-sample Ribo-Seq QC pipeline."""

import os

_CODE_DIR = os.path.dirname(os.path.abspath(__file__))
_REPO = os.path.dirname(os.path.dirname(_CODE_DIR))
_REPO = os.path.dirname(_REPO)

class AnnotationError(RuntimeError):
    pass

def out_dir():
    """Genome QC output root."""
    return os.environ.get("RIBOFLOW_PAPER_QC_OUT",
                          os.path.join(_REPO, "results", "ribo_seq_qc", "genome"))

def tx_out_dir():
    """Transcriptome QC output root, kept distinct so the genome masters are never clobbered."""
    return os.environ.get("RIBOFLOW_PAPER_QC_TX_OUT",
                          os.path.join(_REPO, "results", "ribo_seq_qc", "transcriptome"))

def gtf_path(required=True):
    """GENCODE annotation GTF. No default; not redistributed with this repository."""
    value = os.environ.get("RIBOFLOW_PAPER_GTF")
    if not value and required:
        raise AnnotationError(
            "no GTF configured. Set RIBOFLOW_PAPER_GTF to a GENCODE annotation GTF "
            "(v34 was used for the published cohort), or pass --gtf to "
            "code/make_tables.py. GENCODE is not redistributed with this repository; "
            "download it from https://www.gencodegenes.org/.")
    return value

def appris_path(required=True):
    """APPRIS transcript-lengths TSV. No default; not redistributed here."""
    value = os.environ.get("RIBOFLOW_PAPER_APPRIS")
    if not value and required:
        raise AnnotationError(
            "no APPRIS transcript-lengths table configured. Set RIBOFLOW_PAPER_APPRIS, "
            "or pass --appris to code/make_tables.py. It is not redistributed with this "
            "repository; it comes with the RiboFlow transcriptome reference.")
    return value

def cache_root():
    """`<output root>/.cache` -- internal, safe to delete, rebuilt on the next run."""
    return os.path.join(
        os.environ.get("RIBOFLOW_PAPER_OUT", os.path.join(_REPO, "results")), ".cache")

def cache_dir():
    return os.path.join(cache_root(), "annotation")

def tables_dir():
    return os.path.join(out_dir(), "tables")

_tables_dir = tables_dir

def staging_dir():
    return os.path.join(tables_dir(), "_staging")

def plots_dir():
    return os.path.join(out_dir(), "plots")

# One fingerprinted bundle file, written atomically; never five files reused on existence.

BUNDLE_PAYLOADS = ("appris_cds", "appris_meta", "appris_utr",
                   "appris_gene_body", "all_gene_bodies")

CACHE_SCHEMA_VERSION = 1

def bundle_path():
    """The single cache file. `.pkl` so no parquet engine is required."""
    return os.path.join(cache_dir(), "qc_annotation_bundle.pkl")

def _digest_file(path):
    import hashlib
    digest = hashlib.sha256()
    with open(path, "rb") as handle:
        for chunk in iter(lambda: handle.read(1 << 20), b""):
            digest.update(chunk)
    return digest.hexdigest()

def bundle_fingerprint(gtf=None, appris=None):
    """A content digest of everything that determines the bundle."""
    import hashlib

    gtf = gtf or gtf_path()
    appris = appris or appris_path()
    parts = [
        "schema=%d" % CACHE_SCHEMA_VERSION,
        "gtf=%s" % _digest_file(gtf),
        "appris=%s" % _digest_file(appris),
        "min_five_utr=%d" % MIN_FIVE_UTR,
        "min_cds=%d" % MIN_CDS,
        "min_three_utr=%d" % MIN_THREE_UTR,
        "builder=%s" % _digest_file(os.path.abspath(__file__)),
    ]
    return hashlib.sha256("\n".join(parts).encode()).hexdigest()

MIN_FIVE_UTR = 30
MIN_CDS = 150
MIN_THREE_UTR = 30

MIN_LEN, MAX_LEN = 20, 45

# Uniqueness policy lives in `code/common/bam_inputs.py` (NH == 1 / MAPQ >= 42).

FRAME_COLORS = {0: "#F8766D", 1: "#00BA38", 2: "#619CFF"}

BAM_SUFFIXES = (".post_dedup.bam", ".bam")

def sample_from_bam(path):
    base = os.path.basename(path)
    for suf in BAM_SUFFIXES:
        if base.endswith(suf):
            return base[: -len(suf)]
    return os.path.splitext(base)[0]

def build_annotation_cache(gtf=None, appris=None):
    """Parse the GTF + APPRIS once into the five BUNDLE_PAYLOADS tables and cache them.

    Returns the payload dict and writes it as one fingerprinted bundle.
    """
    import gzip
    import re
    import pandas as pd

    gtf = gtf or gtf_path()
    appris = appris or appris_path()
    os.makedirs(cache_dir(), exist_ok=True)

    meta_rows = []
    with open(appris) as fh:
        for line in fh:
            line = line.strip()
            if not line:
                continue
            fields = line.split("|")
            tx_id = fields[0]
            utr5_len = utr3_len = cds_len = 0
            has_utr5 = has_utr3 = False
            for f in fields[1:]:
                if f.startswith("UTR5:"):
                    c = f[5:].split("-")
                    if len(c) == 2:
                        utr5_len, has_utr5 = int(c[1]) - int(c[0]) + 1, True
                elif f.startswith("CDS:"):
                    c = f[4:].split("-")
                    if len(c) == 2:
                        cds_len = int(c[1]) - int(c[0]) + 1
                elif f.startswith("UTR3:"):
                    c = f[5:].split("-")
                    if len(c) == 2:
                        utr3_len, has_utr3 = int(c[1]) - int(c[0]) + 1, True
            length_filtered = (
                (not has_utr5 or utr5_len < MIN_FIVE_UTR)
                or cds_len < MIN_CDS
                or (not has_utr3 or utr3_len < MIN_THREE_UTR)
            )
            meta_rows.append({
                "transcript_id":   tx_id,
                "appris_cds_len":  cds_len,
                "length_filtered": length_filtered,
            })
    meta_df = pd.DataFrame(meta_rows).drop_duplicates("transcript_id")
    appris_ids = set(meta_df["transcript_id"])

    # gene: all gene types (no filter) — separates intronic from intergenic.
    _FEAT = {"CDS", "UTR", "gene"}
    cds_rows  = []
    utr_rows  = []
    gene_rows = []
    _open = gzip.open if gtf.endswith(".gz") else open
    with _open(gtf, "rt") as fh:
        for line in fh:
            if line[0] == "#":
                continue
            f = line.rstrip("\n").split("\t")
            if len(f) < 9 or f[2] not in _FEAT:
                continue
            start  = int(f[3]) - 1
            end    = int(f[4])
            strand = f[6]
            chrom  = f[0]

            if f[2] == "gene":
                gene_rows.append({"Chromosome": chrom, "Start": start,
                                   "End": end, "Strand": strand})
                continue

            if 'gene_type "protein_coding"' not in f[8]:
                continue
            m_tx = re.search(r'transcript_id "([^"]+)"', f[8])
            if not m_tx or m_tx.group(1) not in appris_ids:
                continue
            tx_id = m_tx.group(1)
            if f[2] == "CDS":
                m_gene = re.search(r'gene_id "([^"]+)"', f[8])
                m_name = re.search(r'gene_name "([^"]+)"', f[8])
                cds_rows.append({
                    "Chromosome":    chrom,
                    "Start":         start,
                    "End":           end,
                    "Strand":        strand,
                    "Phase":         int(f[7]) if f[7] != "." else 0,
                    "transcript_id": tx_id,
                    "gene_id":       m_gene.group(1) if m_gene else "",
                    "gene_name":     m_name.group(1) if m_name else "",
                })
            else:
                utr_rows.append({
                    "Chromosome":    chrom,
                    "Start":         start,
                    "End":           end,
                    "Strand":        strand,
                    "transcript_id": tx_id,
                })

    cds_df = pd.DataFrame(cds_rows)

    if utr_rows:
        utr_df = pd.DataFrame(utr_rows)
        tx_extent = cds_df.groupby("transcript_id").agg(
            min_cds=("Start", "min"),
            max_cds=("End",   "max"),
        )
        utr_df = utr_df.join(tx_extent, on="transcript_id")
        plus = utr_df["Strand"] == "+"
        utr_df["utr_type"] = None
        utr_df.loc[ plus & (utr_df["End"]   <= utr_df["min_cds"]), "utr_type"] = "five"
        utr_df.loc[ plus & (utr_df["Start"] >= utr_df["max_cds"]), "utr_type"] = "three"
        utr_df.loc[~plus & (utr_df["Start"] >= utr_df["max_cds"]), "utr_type"] = "five"
        utr_df.loc[~plus & (utr_df["End"]   <= utr_df["min_cds"]), "utr_type"] = "three"
        utr_df = utr_df.dropna(subset=["utr_type"])
        utr_df = utr_df[["Chromosome", "Start", "End", "Strand",
                          "utr_type", "transcript_id"]].reset_index(drop=True)
    else:
        utr_df = pd.DataFrame(
            columns=["Chromosome", "Start", "End", "Strand", "utr_type", "transcript_id"]
        )

    plus = cds_df["Strand"] == "+"
    plus_start  = cds_df[plus].groupby("transcript_id")["Start"].min()
    minus_start = cds_df[~plus].groupby("transcript_id")["End"].max() - 1
    cds_df["cds_genomic_start"] = cds_df["transcript_id"].map(
        pd.concat([plus_start, minus_start])
    )

    exon_coords = pd.concat([
        cds_df[["Chromosome", "Start", "End", "Strand", "transcript_id"]],
        utr_df[["Chromosome", "Start", "End", "Strand", "transcript_id"]],
    ], ignore_index=True)
    gene_body_df = (
        exon_coords
        .groupby(["Chromosome", "Strand", "transcript_id"], as_index=False)
        .agg(Start=("Start", "min"), End=("End", "max"))
    )[["Chromosome", "Start", "End", "Strand", "transcript_id"]]

    _cols = ["Chromosome", "Start", "End", "Strand"]
    all_gene_bodies_df = pd.DataFrame(gene_rows) if gene_rows else pd.DataFrame(columns=_cols)

    payloads = {"appris_cds": cds_df, "appris_meta": meta_df, "appris_utr": utr_df,
                "appris_gene_body": gene_body_df, "all_gene_bodies": all_gene_bodies_df}
    _write_bundle(payloads, gtf, appris)
    return payloads

def _write_bundle(payloads, gtf, appris):
    """Write the bundle atomically (temp file + os.replace) so readers never see a partial pickle."""
    import pickle
    import tempfile

    target = bundle_path()
    os.makedirs(os.path.dirname(target), exist_ok=True)
    document = {"schema_version": CACHE_SCHEMA_VERSION,
                "fingerprint": bundle_fingerprint(gtf, appris),
                "payloads": payloads}
    handle, temporary = tempfile.mkstemp(dir=os.path.dirname(target), suffix=".tmp")
    try:
        with os.fdopen(handle, "wb") as stream:
            pickle.dump(document, stream, protocol=4)
        os.replace(temporary, target)
    except BaseException:
        if os.path.exists(temporary):
            os.unlink(temporary)
        raise
    return target

def annotation_fingerprint(extra=()):
    """Content digest of the GTF + APPRIS + `extra` terms, for annotation-derived caches
    held outside the bundle (pass the deriving function's source digest in `extra`)."""
    import hashlib

    parts = [bundle_fingerprint()] + [str(e) for e in extra]
    return hashlib.sha256("\n".join(parts).encode()).hexdigest()

def cached_frame(path, fingerprint, build):
    """Read `path` if it carries `fingerprint`, else `build()` and write atomically.

    Absent, unreadable, corrupt or mismatched all mean rebuild.
    """
    import pickle
    import tempfile

    path = str(path)
    if os.path.exists(path):
        try:
            with open(path, "rb") as stream:
                document = pickle.load(stream)
            if document.get("fingerprint") == fingerprint:
                return document["payload"]
        except Exception:
            pass
    payload = build()
    os.makedirs(os.path.dirname(path), exist_ok=True)
    handle, temporary = tempfile.mkstemp(dir=os.path.dirname(path), suffix=".tmp")
    try:
        with os.fdopen(handle, "wb") as stream:
            pickle.dump({"fingerprint": fingerprint, "payload": payload}, stream,
                        protocol=4)
        os.replace(temporary, path)
    except BaseException:
        if os.path.exists(temporary):
            os.unlink(temporary)
        raise
    return payload

def load_bundle(gtf=None, appris=None):
    """The five tables, rebuilding whenever the cached bundle is not the current one
    (absent, unreadable, different schema version, inputs or builder)."""
    import pickle

    path = bundle_path()
    if os.path.exists(path):
        try:
            with open(path, "rb") as stream:
                document = pickle.load(stream)
            if (document.get("schema_version") == CACHE_SCHEMA_VERSION
                    and document.get("fingerprint") == bundle_fingerprint(gtf, appris)
                    and set(document.get("payloads", {})) == set(BUNDLE_PAYLOADS)):
                return document["payloads"]
        except Exception:
            pass
    return build_annotation_cache(gtf=gtf, appris=appris)

def load_annotation(gtf=None, appris=None):
    """The per-CDS-exon annotation table for APPRIS principal isoforms."""
    return load_bundle(gtf, appris)["appris_cds"]

def load_appris_meta(gtf=None, appris=None):
    """The per-transcript APPRIS universe, with the `length_filtered` flag."""
    return load_bundle(gtf, appris)["appris_meta"]

def load_appris_utr(gtf=None, appris=None):
    """The per-UTR-exon table for APPRIS principal isoforms."""
    return load_bundle(gtf, appris)["appris_utr"]

def load_gene_bodies(gtf=None, appris=None):
    """The per-transcript genomic locus for APPRIS principal isoforms."""
    return load_bundle(gtf, appris)["appris_gene_body"]

def load_all_gene_bodies(gtf=None, appris=None):
    """Every GTF gene locus, all gene types, one row per gene."""
    return load_bundle(gtf, appris)["all_gene_bodies"]
