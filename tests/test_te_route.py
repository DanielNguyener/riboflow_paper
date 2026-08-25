"""Figure 4's numbers: the shipped tables carry the published values, and the two R
programs reproduce those tables from the shipped count matrices.

The R layer runs only when `Rscript` is on PATH (base R, no packages). Everything is
written to `tmp_path`; the repository's results/ is never touched.

Run with `python` (3.9).
"""
from __future__ import annotations

import shutil
import subprocess
from pathlib import Path

import pandas as pd
import pytest

REPO = Path(__file__).resolve().parents[1]
TE_ROUTE = REPO / "code" / "te_route"
COUNTS = REPO / "data" / "ribo_rna" / "counts"
TABLES = REPO / "data" / "te_route" / "tables"
ORF_CATALOG = REPO / "data" / "annotation" / "orf_catalog.tsv"

#: The published numbers (manuscript Figure 4; docs/numeric_claims.tsv C13-C18).
N_TRANSCRIPTS = 19736
N_GATED = 11589
N_ESTIMATION = 7864
N_LINES = 24
MEDIAN_DTE = -0.0634
IQR_DTE = (-0.0992, -0.0408)
N_LOWER, N_HIGHER = 10727, 862
N_SIGNIFICANT, N_SIG_NEGATIVE, N_SIG_POSITIVE = 80, 79, 1
MEDIAN_SPEARMAN, MEDIAN_PEARSON = 0.9686, 0.9660


@pytest.fixture(scope="module")
def per_gene():
    return pd.read_csv(TABLES / "per_gene_delta.tsv", sep="\t")


@pytest.fixture(scope="module")
def correlation():
    return pd.read_csv(TABLES / "route_correlation.tsv", sep="\t")


# ── the shipped count matrices ────────────────────────────────────────────────

@pytest.mark.parametrize("name", ["ribo_counts_genome.csv", "ribo_counts_txome.csv",
                                  "rna_counts_genome.csv", "rna_counts_txome.csv"])
def test_each_count_matrix_is_the_full_reference_set_over_the_cohort(name):
    frame = pd.read_csv(COUNTS / name)
    assert frame.columns[0] == "transcript_id"
    assert frame.shape == (N_TRANSCRIPTS, 1 + N_LINES), frame.shape
    assert not frame["transcript_id"].duplicated().any()
    values = frame.drop(columns="transcript_id")
    assert (values.dtypes == "int64").all(), "counts are integers"
    assert (values >= 0).all().all()


def test_the_four_matrices_share_one_index_and_one_column_order():
    frames = [pd.read_csv(COUNTS / n) for n in sorted(COUNTS.glob("*.csv"))]
    assert len(frames) == 4
    for frame in frames[1:]:
        assert frame["transcript_id"].equals(frames[0]["transcript_id"])
        assert list(frame.columns) == list(frames[0].columns)


# ── the shipped statistics tables ─────────────────────────────────────────────

def test_the_gate_and_the_published_delta_te_summary(per_gene):
    assert len(per_gene) == N_GATED
    assert not per_gene["transcript_id"].duplicated().any()
    assert per_gene["n_lines"].between(1, N_LINES).all()
    dte = per_gene["dte_mean"]
    assert round(float(dte.median()), 4) == MEDIAN_DTE
    assert tuple(round(float(v), 4) for v in dte.quantile([0.25, 0.75])) == IQR_DTE
    assert int((dte < 0).sum()) == N_LOWER
    assert int((dte > 0).sum()) == N_HIGHER


def test_the_eighty_significant_transcripts(per_gene):
    big = per_gene["dte_padj"].notna() & (per_gene["dte_padj"] < 0.05) \
        & (per_gene["dte_mean"].abs() > 1)
    assert int(big.sum()) == N_SIGNIFICANT
    assert int((big & (per_gene["dte_mean"] < 0)).sum()) == N_SIG_NEGATIVE
    assert int((big & (per_gene["dte_mean"] > 0)).sum()) == N_SIG_POSITIVE


def test_delta_te_is_delta_ribo_minus_delta_rna(per_gene):
    """The identity the assay plane (panel C) draws as y = x."""
    delta = per_gene["dribo_mean"] - per_gene["drna_mean"] - per_gene["dte_mean"]
    assert float(delta.abs().max()) < 1e-9


def test_the_confidence_interval_brackets_the_mean(per_gene):
    assert (per_gene["dte_ci_low"] <= per_gene["dte_mean"]).all()
    assert (per_gene["dte_mean"] <= per_gene["dte_ci_high"]).all()


def test_gene_names_come_from_the_orf_catalog(per_gene):
    catalog = pd.read_csv(ORF_CATALOG, sep="\t").drop_duplicates("transcript_id")
    names = catalog.set_index("transcript_id")["gene_name"]
    merged = per_gene.set_index("transcript_id")["gene_name"]
    assert merged.index.isin(names.index).all()
    assert (merged == names.reindex(merged.index)).all()
    for gene in ("GAPDH", "COMT"):
        assert (merged == gene).sum() == 1, gene


def test_the_route_correlation_table(correlation):
    assert len(correlation) == N_LINES
    assert not correlation["sample"].duplicated().any()
    assert "HeLa" in set(correlation["sample"])
    assert round(float(correlation["spearman_rho"].median()), 4) == MEDIAN_SPEARMAN
    assert round(float(correlation["pearson_r"].median()), 4) == MEDIAN_PEARSON
    assert correlation["n_transcripts"].between(1, N_GATED).all()


# ── the R programs reproduce the shipped tables ───────────────────────────────

RSCRIPT = shutil.which("Rscript")


def run_r(script, *args):
    return subprocess.run([RSCRIPT, str(TE_ROUTE / script), *map(str, args)],
                          capture_output=True, text=True, cwd=str(REPO))


@pytest.mark.skipif(not RSCRIPT, reason="Rscript is not on PATH")
def test_the_r_programs_print_help_and_refuse_odd_arguments(tmp_path):
    for script in ("normalization.R", "te_statistics.R"):
        result = run_r(script, "--help")
        assert result.returncode == 0 and "usage:" in result.stdout, result.stderr
        result = run_r(script, "--output")
        assert result.returncode != 0
        result = run_r(script, "--nope", str(tmp_path))
        assert result.returncode != 0 and "unknown argument" in result.stderr


@pytest.fixture(scope="module")
def rerun(tmp_path_factory):
    if not RSCRIPT:
        pytest.skip("Rscript is not on PATH")
    out = tmp_path_factory.mktemp("te_route")
    normalized, tables = out / "normalized", out / "tables"
    result = run_r("normalization.R", "--counts", COUNTS, "--output", normalized)
    assert result.returncode == 0, result.stdout + result.stderr
    result = run_r("te_statistics.R", "--normalized", normalized,
                   "--orf-catalog", ORF_CATALOG, "--output", tables)
    assert result.returncode == 0, result.stdout + result.stderr
    return normalized, tables, result.stdout


def test_normalization_asserts_its_two_counts_and_writes_one_factor_per_library(rerun):
    normalized, _tables, _log = rerun
    factors = pd.read_csv(normalized / "size_factors.csv")
    assert list(factors.columns) == ["sample", "ribo", "rna"]
    assert len(factors) == N_LINES
    assert (factors[["ribo", "rna"]] > 0).all().all()
    for name in ("ribo_scaled_genome.csv", "ribo_scaled_txome.csv",
                 "rna_scaled_genome.csv", "rna_scaled_txome.csv"):
        scaled = pd.read_csv(normalized / name)
        assert scaled.shape == (N_GATED, 1 + N_LINES), name


def test_the_statistics_reproduce_the_shipped_tables_exactly(rerun):
    _normalized, tables, log = rerun
    for name in ("per_gene_delta.tsv", "route_correlation.tsv"):
        assert (tables / name).read_bytes() == (TABLES / name).read_bytes(), \
            "%s differs from the shipped table" % name
    assert "padj < 0.05 and |delta TE| > 1: 80  (79 negative, 1 positive)" in log
    assert "median Spearman 0.9686  median Pearson 0.9660" in log
