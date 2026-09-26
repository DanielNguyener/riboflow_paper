"""Raw-input resolution for programs that open BAMs: flag > RIBOFLOW_PAPER_* > config/local.yaml.

Nothing has a built-in default; a missing input names all three ways to supply it.
Functions only -- the environment is read when called, never at import.
"""
from __future__ import annotations

import hashlib
import os
from pathlib import Path

REPO = Path(__file__).resolve().parents[2]
LOCAL_CONFIG = REPO / "config" / "local.yaml"


def die(message):
    raise SystemExit("error: %s" % message)


def make_log(prefix):
    """`log("...")` -> `[prefix] ...` on stdout, flushed."""
    def log(message):
        print("[%s] %s" % (prefix, message), flush=True)
    return log


def import_from(directory, name, extra=()):
    """Import module `name` from `directory` (+ `extra` dirs) without leaving sys.path changed."""
    import sys
    saved = list(sys.path)
    try:
        for entry in [str(directory)] + [str(e) for e in extra]:
            if entry not in sys.path:
                sys.path.insert(0, entry)
        return __import__(name)
    finally:
        sys.path[:] = saved + [p for p in sys.path if p not in saved]


def repo_path(relative):
    """A repository-relative path made absolute; absolute paths pass through."""
    path = Path(relative)
    return path if path.is_absolute() else REPO / path


def _local_config():
    if not LOCAL_CONFIG.exists():
        return {}
    import yaml
    with open(LOCAL_CONFIG) as handle:
        return yaml.safe_load(handle) or {}


def resolve_external_inputs(bams=None, gtf=None, appris=None, sample="HeLa"):
    """{ribo_genome, ribo_txome, gtf, appris} as existing paths, or a hard stop."""
    local = _local_config()
    bams = bams or os.environ.get("RIBOFLOW_PAPER_BAMS") or local.get("bams")
    gtf = gtf or os.environ.get("RIBOFLOW_PAPER_GTF") or local.get("gtf")
    appris = appris or os.environ.get("RIBOFLOW_PAPER_APPRIS") or local.get("appris")
    missing = [name for name, value in
               (("--bams / RIBOFLOW_PAPER_BAMS / local.yaml:bams", bams),
                ("--gtf / RIBOFLOW_PAPER_GTF / local.yaml:gtf", gtf),
                ("--appris / RIBOFLOW_PAPER_APPRIS / local.yaml:appris", appris))
               if not value]
    if missing:
        die("raw inputs not configured: %s (see config/inputs.example.yaml)"
            % "; ".join(missing))
    bam_inputs = import_from(Path(__file__).resolve().parent, "bam_inputs")
    paths = {"ribo_genome": os.path.join(
                 str(bams), bam_inputs.GENOME_BAM_TEMPLATE.format(s=sample)),
             "ribo_txome": os.path.join(
                 str(bams), bam_inputs.TXOME_BAM_TEMPLATE.format(s=sample)),
             "gtf": str(gtf), "appris": str(appris)}
    for key in ("ribo_genome", "ribo_txome", "gtf", "appris"):
        if not os.path.exists(paths[key]):
            die("%s does not exist: %s" % (key, paths[key]))
    return paths


def sha256_of(path, chunk=1 << 22):
    digest = hashlib.sha256()
    with open(path, "rb") as handle:
        for block in iter(lambda: handle.read(chunk), b""):
            digest.update(block)
    return digest.hexdigest()
