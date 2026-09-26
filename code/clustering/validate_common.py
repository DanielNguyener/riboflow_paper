"""What the three `validate_cluster_*.py` scripts share: constants, the GTF line reader,
interval merging and the per-cluster statistics. They write tables only; the Figure 6
validation row is drawn from those tables by `panels/plot_read_fate_clusters.py`."""
from __future__ import annotations

import argparse
import gzip
import itertools
import os
import re
import sys
from pathlib import Path

import numpy as np
import pandas as pd

HERE = Path(__file__).resolve().parent
REPO = HERE.parents[1]
sys.path.insert(0, str(REPO / "code" / "common"))
import inputs  # noqa: E402

KEY = ["gene", "gene_id", "transcript_id"]
ATTR = re.compile(r'(\S+) "([^"]*)"')


def describe(path, **extra):
    path = Path(path)
    return {"path": str(path), **extra, "bytes": path.stat().st_size,
            "sha256": inputs.sha256_of(path)}


# ── command line ─────────────────────────────────────────────────────────────────────────
def parser(doc):
    p = argparse.ArgumentParser(description=doc, formatter_class=argparse.RawDescriptionHelpFormatter)
    p.add_argument("--clusters", required=True, type=Path, help="<stem>.clusters_k<K>.tsv")
    p.add_argument("--output", required=True, type=Path)
    p.add_argument("--gtf", type=Path, help="GENCODE v34 GTF (default: env / config/local.yaml)")
    return p


def resolve_gtf(a):
    gtf = a.gtf or os.environ.get("RIBOFLOW_PAPER_GTF") or inputs._local_config().get("gtf")
    if not gtf or not Path(gtf).exists():
        inputs.die("GTF not configured or missing: --gtf / RIBOFLOW_PAPER_GTF / local.yaml:gtf")
    return gtf


def stem_of(clusters_path):
    return clusters_path.name.split(".clusters_k")[0]


def load_clusters(clusters_path):
    """(stem, cluster table) from a `<stem>.clusters_k<K>.tsv` path."""
    return stem_of(clusters_path), pd.read_csv(clusters_path, sep="\t")


# ── annotation ───────────────────────────────────────────────────────────────────────────
def gtf_records(gtf, features):
    """The tab-split fields of every `features` line of a (gzipped) GTF."""
    opener = gzip.open if str(gtf).endswith(".gz") else open
    with opener(gtf, "rt") as handle:
        for line in handle:
            if line.startswith("#"):
                continue
            f = line.rstrip("\n").split("\t")
            if len(f) >= 9 and f[2] in features:
                yield f


def merge(intervals):
    """Union of half-open [start, end) intervals as a sorted, non-overlapping tuple list."""
    out = []
    for s, e in sorted(intervals):
        if out and s <= out[-1][1]:
            out[-1][1] = max(out[-1][1], e)
        else:
            out.append([s, e])
    return [tuple(x) for x in out]


def merged_length(intervals):
    return sum(end - start for start, end in merge(intervals))


# ── statistics ───────────────────────────────────────────────────────────────────────────
def quartiles(values):
    n = len(values)
    return {"q25": np.percentile(values, 25) if n else np.nan,
            "median": np.median(values) if n else np.nan,
            "q75": np.percentile(values, 75) if n else np.nan,
            "mean": values.mean() if n else np.nan}


def odds_ratio(inside, outside):
    """inside/outside: boolean Series. Woolf CI with Haldane 0.5 when a cell is 0."""
    from scipy.stats import fisher_exact
    a, b = int(inside.sum()), int((~inside).sum())
    c, d = int(outside.sum()), int((~outside).sum())
    table = np.array([[a, b], [c, d]])
    _stat, p = fisher_exact(table)
    cells = table.astype(float)
    if (cells == 0).any():
        cells = cells + 0.5
    orr = (cells[0, 0] * cells[1, 1]) / (cells[0, 1] * cells[1, 0])
    se = np.sqrt((1.0 / cells).sum())
    return orr, np.exp(np.log(orr) - 1.96 * se), np.exp(np.log(orr) + 1.96 * se), p, (a, b, c, d)


def holm(pvalues):
    order = np.argsort(pvalues)
    adjusted = np.empty(len(pvalues))
    running = 0.0
    for rank, idx in enumerate(order):
        running = max(running, min(1.0, pvalues[idx] * (len(pvalues) - rank)))
        adjusted[idx] = running
    return adjusted


def cluster_tests(groups, clusters):
    """Kruskal-Wallis over `groups` ({cluster: values}) and Holm-adjusted pairwise
    Mann-Whitney U, as rows for a tests table."""
    from scipy.stats import kruskal, mannwhitneyu
    stat, p = kruskal(*[groups[c] for c in clusters])
    tests = [{"test": "kruskal_wallis", "a": "all", "b": "all", "statistic": stat, "p": p, "p_holm": p}]
    raw = []
    for a, b in itertools.combinations(clusters, 2):
        u, pu = mannwhitneyu(groups[a], groups[b], alternative="two-sided")
        raw.append((a, b, u, pu))
    for (a, b, u, pu), ph in zip(raw, holm([r[3] for r in raw])):
        tests.append({"test": "mann_whitney_u", "a": str(a), "b": str(b), "statistic": u, "p": pu, "p_holm": ph})
    return tests
