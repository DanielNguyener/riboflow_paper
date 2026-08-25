"""`build_s1_table.py`: the selection rule, the quality override, and byte-identity.

Two layers, as everywhere else here:

  * **synthetic** -- a hand-written QC table and a hand-written xlsx, small enough that the
    expected panel can be reasoned about directly. Every filter, the score, the ordering,
    the override, the declared-curation guard and each failure message. No R, no RiboBase.
  * **real** -- the actual `.rda` and the actual supplement, asserting the produced CSV is
    byte-identical to the published `samples.csv`. **Skipped** unless both
    are supplied:

        RIBOBASER_RDA=/path/to/ribobaser/data/Ribobase_QC_dedup_data.rda \\
        RIBOFLOW_PAPER_S1_XLSX=/path/to/41587_2025_2718_MOESM3_ESM.xlsx \\
            python -m pytest tests/test_s1_table.py -q

The synthetic layer monkeypatches `CURATED_ANNOTATION`, because that table is keyed by the
real panel's GSMs. That is the point of it being declared data: it is not derivable, so a
test cannot conjure it either.

Run with `python` (3.9).
"""
from __future__ import annotations

import hashlib
import importlib.util
import os
import subprocess
import sys
from pathlib import Path

import pandas as pd
import pytest

REPO = Path(__file__).resolve().parents[1]
GENERATOR = REPO / "supporting_information" / "S1_Table" / "build_s1_table.py"
PUBLISHED = REPO / "supporting_information" / "S1_Table" / "samples.csv"

RDA = os.environ.get("RIBOBASER_RDA")
XLSX = os.environ.get("RIBOFLOW_PAPER_S1_XLSX")


def _load():
    name = "build_s1_table"
    if name in sys.modules:
        return sys.modules[name]
    spec = importlib.util.spec_from_file_location(name, str(GENERATOR))
    module = importlib.util.module_from_spec(spec)
    sys.modules[name] = module
    spec.loader.exec_module(module)
    return module


@pytest.fixture(scope="module")
def s1():
    return _load()


# ── the synthetic universe ────────────────────────────────────────────────────
# Six samples. What each one is for:
#
#   GSM_A  cell line A, deepest, run of 4 -> the top pick
#   GSM_B  cell line B, shallower, run of 5 -> second
#   GSM_A2 cell line A again, weaker -> must lose to GSM_A (one row per cell line)
#   GSM_C  cell line C, qualifies, lowest score -> the override's drop target
#   GSM_D  cell line D, BELOW the depth floor -> filtered out
#   GSM_E  cell line E, qualifies on depth but has no matched RNA -> filtered out
#   GSM_F  cell line F, deep, but only 2 consecutive lengths above 0.60 -> filtered out,
#          and it is what the synthetic override adds back (the MCF10A analogue)
SYNTHETIC_QC = [
    # Experiment, Study,  Cell line, Species, Start, End, PeriScore, distr, reads
    ("GSM_A", "GSE_A", "LineA", "human", 26, 31, 0.700000111, "[0.7, 0.72, 0.65, 0.61, 0.4]", 9_000_000),
    ("GSM_B", "GSE_B", "LineB", "human", 25, 30, 0.680000222, "[0.61, 0.62, 0.63, 0.64, 0.65]", 3_000_000),
    ("GSM_A2", "GSE_A", "LineA", "human", 26, 31, 0.640000333, "[0.61, 0.62, 0.63, 0.4]", 2_000_000),
    ("GSM_C", "GSE_C", "LineC", "human", 24, 29, 0.620000444, "[0.61, 0.62, 0.63, 0.4]", 1_600_000),
    ("GSM_D", "GSE_D", "LineD", "human", 26, 31, 0.900000555, "[0.9, 0.91, 0.92, 0.93]", 900_000),
    ("GSM_E", "GSE_E", "LineE", "human", 26, 31, 0.700000666, "[0.7, 0.71, 0.72, 0.73]", 5_000_000),
    ("GSM_F", "GSE_F", "LineF", "human", 26, 28, 0.950000777, "[0.5, 0.96, 0.97]", 1_450_000),
    ("GSM_M", "GSE_M", "LineM", "mouse", 26, 31, 0.990000888, "[0.9, 0.91, 0.92, 0.93]", 8_000_000),
]
QC_COLUMNS = ["Experiment", "Study", "Cell line", "Species", "Start length", "End length",
              "Periodicity score", "Periodicity distr", "Read counts (dynamic)"]
NO_RNA = {"GSM_E"}


def _qc_frame():
    frame = pd.DataFrame(SYNTHETIC_QC, columns=QC_COLUMNS)
    # the three coverage columns and the two other depth windows, derived so they are
    # distinguishable per row but need no separate table
    frame["CDS coverage (15-40)"] = 0.90 + frame.index * 0.001111
    frame["CDS coverage (27-30)"] = 0.91 + frame.index * 0.001111
    frame["CDS coverage (dynamic)"] = 0.92 + frame.index * 0.001111
    frame["Read counts (15-40)"] = frame["Read counts (dynamic)"] + 100_000
    frame["Read counts (27-30)"] = frame["Read counts (dynamic)"] - 100_000
    return frame


@pytest.fixture
def qc_csv(tmp_path):
    path = tmp_path / "qc.csv"
    _qc_frame().to_csv(path, index=False)
    return path


@pytest.fixture
def xlsx(tmp_path):
    """A supplement with the sheet name, the second-row header, and the three columns."""
    rows = []
    for record in SYNTHETIC_QC:
        gsm = record[0]
        rows.append({
            "experiment_alias": gsm,
            "matched_RNA-seq_experiment_alias": None if gsm in NO_RNA else gsm + "_rna",
            "cell_info_GEO": "cell line: %s | note: synthetic" % record[2],
        })
    frame = pd.DataFrame(rows)
    path = tmp_path / "supplement.xlsx"
    with pd.ExcelWriter(path, engine="openpyxl") as writer:
        # the real sheet carries a title row above the header; startrow=1 reproduces that
        frame.to_excel(writer, sheet_name="S1_RiboBase_metadata", index=False, startrow=1)
    return path


@pytest.fixture
def curated(s1, monkeypatch):
    annotation = {gsm: ("tissue_%s" % gsm[-1], "type_%s" % gsm[-1], "disease, with comma")
                  for gsm in ("GSM_A", "GSM_B", "GSM_C", "GSM_F")}
    monkeypatch.setattr(s1, "CURATED_ANNOTATION", annotation)
    return annotation


@pytest.fixture
def synthetic_override(s1, monkeypatch):
    """Drop LineC (weakest qualifier), add LineF (excluded by the run-length rule)."""
    monkeypatch.setattr(s1, "OVERRIDE_DROP", "GSM_C")
    monkeypatch.setattr(s1, "OVERRIDE_ADD", "GSM_F")


def _build(s1, qc_csv, xlsx, **kwargs):
    return s1.build_table(s1.load_qc(qc_csv=qc_csv), s1.load_matched_rna(xlsx), **kwargs)


# ── the pieces of the rule ────────────────────────────────────────────────────
def test_parse_distr_round_trips_the_stored_format(s1):
    assert s1.parse_distr("[0.68812, 0.6241, 0.60067]") == [0.68812, 0.6241, 0.60067]
    assert s1.parse_distr("[0.5]") == [0.5]


@pytest.mark.parametrize("distr,threshold,expected", [
    ("[0.7, 0.72, 0.65, 0.61, 0.4]", 0.60, 4),
    ("[0.61, 0.5, 0.62, 0.63, 0.64]", 0.60, 3),          # the run RESTARTS after a dip
    ("[0.50697, 0.90148, 0.87483]", 0.60, 2),            # MCF10A's actual vector
    ("[0.50697, 0.90148, 0.87483]", 0.50, 3),            # ... and why the floor matters
    ("[0.4, 0.4]", 0.60, 0),
    ("[0.6, 0.6]", 0.60, 2),                             # >= is inclusive
])
def test_max_consecutive_above_counts_runs_not_totals(s1, distr, threshold, expected):
    assert s1.max_consecutive_above(distr, threshold) == expected


def test_pool_applies_every_hard_filter(s1, qc_csv, xlsx):
    pool = s1.build_pool(s1.load_qc(qc_csv=qc_csv), s1.load_matched_rna(xlsx))
    assert set(pool["Experiment"]) == {"GSM_A", "GSM_B", "GSM_A2", "GSM_C"}
    # GSM_M non-human, GSM_D below MIN_READS, GSM_E no matched RNA, GSM_F run of 2


def test_norm_depth_is_min_max_over_the_pool(s1, qc_csv, xlsx):
    pool = s1.build_pool(s1.load_qc(qc_csv=qc_csv), s1.load_matched_rna(xlsx))
    assert pool["norm_depth"].min() == pytest.approx(0.0)
    assert pool["norm_depth"].max() == pytest.approx(1.0)
    deepest = pool.loc[pool["Read counts (dynamic)"].idxmax(), "Experiment"]
    assert deepest == "GSM_A"


def test_one_row_per_cell_line_keeps_the_best(s1, qc_csv, xlsx):
    pool = s1.build_pool(s1.load_qc(qc_csv=qc_csv), s1.load_matched_rna(xlsx))
    panel = s1.select_panel(pool)
    assert list(panel["Cell line"]) == sorted(set(panel["Cell line"]),
                                              key=list(panel["Cell line"]).index)
    assert panel["Cell line"].duplicated().sum() == 0
    assert set(panel["Experiment"]) == {"GSM_A", "GSM_B", "GSM_C"}    # not GSM_A2
    assert panel["score"].is_monotonic_decreasing


def test_tied_scores_are_refused_rather_than_ordered_arbitrarily(s1, qc_csv, xlsx,
                                                                 monkeypatch):
    """`groupby(...).first()` on tied scores is not deterministic; the pool asserts it."""
    frame = _qc_frame()
    frame.loc[frame.Experiment == "GSM_B", "Read counts (dynamic)"] = \
        frame.loc[frame.Experiment == "GSM_A", "Read counts (dynamic)"].iloc[0]
    frame.loc[frame.Experiment == "GSM_B", "Periodicity distr"] = \
        frame.loc[frame.Experiment == "GSM_A", "Periodicity distr"].iloc[0]
    path = qc_csv.parent / "tied.csv"
    frame.to_csv(path, index=False)
    with pytest.raises(SystemExit) as error:
        s1.build_pool(s1.load_qc(qc_csv=path), s1.load_matched_rna(xlsx))
    assert "tied composite score" in str(error.value)


# ── the assembled table ───────────────────────────────────────────────────────
def test_table_has_the_published_shape_and_column_order(s1, qc_csv, xlsx, curated,
                                                        synthetic_override):
    table = _build(s1, qc_csv, xlsx)
    assert list(table.columns) == s1.COLUMNS
    assert len(table.columns) == 23
    assert list(table["ribo_GSM"]) == ["GSM_A", "GSM_B", "GSM_F"]     # override applied


def test_qc_passthrough_is_verbatim_full_precision(s1, qc_csv, xlsx, curated,
                                                   synthetic_override):
    """The 10 passthrough columns are the .rda's own values, not rounded copies."""
    table = _build(s1, qc_csv, xlsx).set_index("ribo_GSM")
    source = _qc_frame().set_index("Experiment")
    for gsm in ("GSM_A", "GSM_B"):
        for column in s1.QC_PASSTHROUGH:
            assert table.loc[gsm, column] == source.loc[gsm, column]
    assert table.loc["GSM_A", "Periodicity score"] == 0.700000111
    assert table.loc["GSM_A", "periodicity"] == 0.7                   # the rounded twin


def test_the_override_row_is_derived_not_declared(s1, qc_csv, xlsx, curated,
                                                  synthetic_override):
    """Only WHICH sample is declared. Its numbers are computed, including norm_depth."""
    import numpy as np
    table = _build(s1, qc_csv, xlsx).set_index("ribo_GSM")
    added = table.loc["GSM_F"]
    source = _qc_frame().set_index("Experiment").loc["GSM_F"]

    assert added["max_consec"] == 2                                   # computed, not typed
    assert added["Periodicity score"] == round(source["Periodicity score"], 3)
    assert added["Periodicity score"] != source["Periodicity score"]  # rounded, unlike the rest

    kept = table.drop("GSM_F")
    depth = np.log10(kept["Read counts (dynamic)"].astype(float))
    expected = ((np.log10(float(source["Read counts (dynamic)"])) - depth.min())
                / (depth.max() - depth.min()))
    assert added["norm_depth"] == pytest.approx(round(expected, 3))


def test_the_override_row_norm_depth_is_negative_when_it_is_the_shallowest(
        s1, qc_csv, xlsx, curated, synthetic_override):
    """The published -0.15: a min-max value below 0 can only mean an out-of-pool row."""
    table = _build(s1, qc_csv, xlsx).set_index("ribo_GSM")
    assert table.loc["GSM_F", "norm_depth"] < 0
    assert (table.drop("GSM_F")["norm_depth"] >= 0).all()


def test_no_override_yields_the_automatic_panel(s1, qc_csv, xlsx, curated,
                                                synthetic_override):
    table = _build(s1, qc_csv, xlsx, apply_override=False)
    assert list(table["ribo_GSM"]) == ["GSM_A", "GSM_B", "GSM_C"]


def test_override_refuses_when_its_drop_target_was_not_selected(s1, qc_csv, xlsx, curated,
                                                                monkeypatch):
    """If the inputs change under it, the override must fail loudly, not silently skip."""
    monkeypatch.setattr(s1, "OVERRIDE_DROP", "GSM_NOT_SELECTED")
    monkeypatch.setattr(s1, "OVERRIDE_ADD", "GSM_F")
    with pytest.raises(SystemExit) as error:
        _build(s1, qc_csv, xlsx)
    assert "did not choose it" in str(error.value)


def test_missing_curation_is_named_not_filled_in(s1, qc_csv, xlsx, synthetic_override,
                                                 monkeypatch):
    monkeypatch.setattr(s1, "CURATED_ANNOTATION", {"GSM_A": ("t", "c", "d")})
    with pytest.raises((SystemExit, KeyError)) as error:
        _build(s1, qc_csv, xlsx)
    assert "GSM_B" in str(error.value) or "GSM_B" in repr(error.value)


def test_integer_columns_stay_integers(s1, qc_csv, xlsx, curated, synthetic_override,
                                       tmp_path):
    """Concatenating the override row must not float-ify the read counts."""
    table = _build(s1, qc_csv, xlsx)
    for column in s1.INTEGER_COLUMNS:
        assert str(table[column].dtype) == "int64", column
    path = tmp_path / "out.csv"
    table.to_csv(path, index=False)
    written = pd.read_csv(path, dtype=str).set_index("ribo_GSM")
    for gsm in ("GSM_A", "GSM_F"):
        for column in s1.INTEGER_COLUMNS:
            assert "." not in written.loc[gsm, column], (gsm, column)
    assert written.loc["GSM_A", "Read counts (dynamic)"] == "9000000"


# ── inputs and their failure messages ─────────────────────────────────────────
def test_qc_missing_columns_are_all_named_at_once(s1, tmp_path):
    path = tmp_path / "bad.csv"
    pd.DataFrame({"Experiment": ["GSM_A"], "Species": ["human"]}).to_csv(path, index=False)
    with pytest.raises(SystemExit) as error:
        s1.load_qc(qc_csv=path)
    message = str(error.value)
    assert "Cell line" in message and "Periodicity distr" in message


def test_duplicate_experiment_ids_are_refused(s1, tmp_path):
    frame = _qc_frame()
    frame = pd.concat([frame, frame.iloc[[0]]], ignore_index=True)
    path = tmp_path / "dup.csv"
    frame.to_csv(path, index=False)
    with pytest.raises(SystemExit) as error:
        s1.load_qc(qc_csv=path)
    assert "duplicate Experiment" in str(error.value)


def test_no_rda_and_no_csv_says_where_the_rda_comes_from(s1):
    with pytest.raises(SystemExit) as error:
        s1.load_qc(rda=None, qc_csv=None)
    message = str(error.value)
    assert "RIBOBASER_RDA" in message and "ribobaser" in message


def test_a_missing_rda_is_reported_by_path(s1, tmp_path):
    with pytest.raises(SystemExit) as error:
        s1.load_qc(rda=tmp_path / "nope.rda")
    assert "no such .rda" in str(error.value)


def test_xlsx_sheet_columns_are_validated(s1, tmp_path):
    path = tmp_path / "wrong.xlsx"
    with pd.ExcelWriter(path, engine="openpyxl") as writer:
        pd.DataFrame({"experiment_alias": ["GSM_A"]}).to_excel(
            writer, sheet_name="S1_RiboBase_metadata", index=False, startrow=1)
    with pytest.raises(SystemExit) as error:
        s1.load_matched_rna(path)
    assert "matched_RNA-seq_experiment_alias" in str(error.value)


def test_cli_refuses_to_overwrite_without_force(s1, qc_csv, xlsx, curated,
                                                synthetic_override, tmp_path):
    output = tmp_path / "s1.csv"
    output.write_text("existing\n")
    with pytest.raises(SystemExit) as error:
        s1.main(["--qc-csv", str(qc_csv), "--xlsx", str(xlsx), "--output", str(output)])
    assert "refusing to overwrite" in str(error.value)


def test_cli_verify_reports_a_difference_instead_of_just_failing(s1, qc_csv, xlsx, curated,
                                                                 synthetic_override,
                                                                 tmp_path, capsys):
    output = tmp_path / "s1.csv"
    reference = tmp_path / "reference.csv"
    assert s1.main(["--qc-csv", str(qc_csv), "--xlsx", str(xlsx),
                    "--output", str(output)]) == 0
    reference.write_text(output.read_text().replace("tissue_A", "tissue_Z"))
    code = s1.main(["--qc-csv", str(qc_csv), "--xlsx", str(xlsx), "--output", str(output),
                    "--force", "--verify", str(reference)])
    assert code == 1
    captured = capsys.readouterr()
    assert "GSM_A" in captured.err and "tissue" in captured.err


def test_cli_is_deterministic(s1, qc_csv, xlsx, curated, synthetic_override, tmp_path):
    digests = []
    for index in range(2):
        output = tmp_path / ("run%d.csv" % index)
        s1.main(["--qc-csv", str(qc_csv), "--xlsx", str(xlsx), "--output", str(output)])
        digests.append(hashlib.sha256(output.read_bytes()).hexdigest())
    assert digests[0] == digests[1]


# ── the published table ───────────────────────────────────────────────────────
def test_the_published_csv_still_has_the_recorded_checksum(s1):
    """If this fails, the reference moved and every claim about it needs re-checking."""
    assert hashlib.sha256(PUBLISHED.read_bytes()).hexdigest() == s1.PUBLISHED_SHA256


def test_the_curated_annotation_covers_exactly_the_published_panel(s1):
    published = pd.read_csv(PUBLISHED)
    assert set(s1.CURATED_ANNOTATION) == set(published["ribo_GSM"])
    for _, row in published.iterrows():
        assert s1.CURATED_ANNOTATION[row["ribo_GSM"]] == (row["tissue"], row["cell_type"],
                                                          row["disease"])


def test_the_curated_annotation_agrees_with_docs_accessions(s1):
    """Two published copies of the same curation; they must not drift apart."""
    accessions = pd.read_csv(REPO / "docs" / "accessions.tsv", sep="\t")
    for _, row in accessions.iterrows():
        assert s1.CURATED_ANNOTATION[row["ribo_GSM"]] == (row["tissue"], row["cell_type"],
                                                          row["disease"])


@pytest.mark.skipif(not (RDA and XLSX),
                    reason="set RIBOBASER_RDA and RIBOFLOW_PAPER_S1_XLSX to run the real "
                           "end-to-end S1 Table reproduction")
def test_real_inputs_reproduce_the_published_table_byte_for_byte(tmp_path):
    """The claim the whole script exists to support."""
    output = tmp_path / "samples.csv"
    result = subprocess.run(
        [sys.executable, str(GENERATOR), "--rda", RDA, "--xlsx", XLSX,
         "--output", str(output), "--verify"],
        capture_output=True, text=True)
    assert result.returncode == 0, result.stdout + result.stderr
    assert output.read_bytes() == PUBLISHED.read_bytes()
    assert hashlib.sha256(output.read_bytes()).hexdigest() == \
        "19a3811f73f6217c922a4054f6830e8f0f66dabc5df55026f353e4a3cc8a2e9b"


@pytest.mark.skipif(not (RDA and XLSX), reason="needs the real .rda and supplement")
def test_real_inputs_place_the_override_sample_where_the_docstring_says(tmp_path):
    """MCF10A is 2nd of 799 on periodicity -- the stated reason for the override."""
    s1 = _load()
    qc_dump = tmp_path / "qc.csv"
    qc = s1.load_qc(rda=RDA, dump_to=qc_dump)
    metadata = s1.load_matched_rna(XLSX)

    universe = qc[qc["Species"] == "human"].copy()
    universe["rna"] = universe["Experiment"].map(
        metadata["matched_RNA-seq_experiment_alias"])
    universe = universe[universe["rna"].notna() & (universe["rna"].astype(str) != "nan")]
    universe = universe[universe["Read counts (dynamic)"] >= 1_000_000]
    ranked = universe.sort_values("Periodicity score", ascending=False).reset_index(drop=True)
    rank = int(ranked.index[ranked["Experiment"] == s1.OVERRIDE_ADD][0]) + 1

    assert rank == 2, "MCF10A's periodicity rank moved: %d of %d" % (rank, len(ranked))
    assert len(ranked) == 799
    dropped = qc.set_index("Experiment").loc[s1.OVERRIDE_DROP, "Periodicity score"]
    added = qc.set_index("Experiment").loc[s1.OVERRIDE_ADD, "Periodicity score"]
    assert added > dropped
