#!/usr/bin/env python3
"""Produce the S1 Table (`samples.csv`) from its two sources."""
from __future__ import annotations

import argparse
import hashlib
import os
import subprocess
import sys
import tempfile
from pathlib import Path

import numpy as np
import pandas as pd

REPO = Path(__file__).resolve().parents[2]
PUBLISHED = REPO / "supporting_information" / "S1_Table" / "samples.csv"
PUBLISHED_SHA256 = "19a3811f73f6217c922a4054f6830e8f0f66dabc5df55026f353e4a3cc8a2e9b"

SHEET = "S1_RiboBase_metadata"
SHEET_HEADER_ROW = 1

MIN_READS = 1_500_000
DISTR_THR = 0.60
DISTR_N = 3
W_DEP = 0.5
W_DISTR = 0.5

OVERRIDE_DROP = "GSM4483531"
OVERRIDE_ADD = "GSM1608276"
OVERRIDE_ROUND = 3

CURATED_ANNOTATION = {
    "GSM1248729": ("bone", "cancer cell line", "osteosarcoma"),
    "GSM4793173": ("kidney", "cancer cell line", "normal (transformed)"),
    "GSM1187138": ("bone", "cybrid cell line (osteosarcoma \u00d7 patient enucleated cells)",
                   "osteosarcoma (143B-derived cybrid, patient mtDNA)"),
    "GSM4192110": ("embryo", "embryonic stem cell", "normal"),
    "GSM2100602": ("cervix", "cancer cell line", "cervical adenocarcinoma"),
    "GSM2760255": ("skeletal muscle", "cancer cell line", "rhabdomyosarcoma"),
    "GSM4192441": ("brain", "cancer cell line", "glioblastoma (Grade IV)"),
    "GSM2838831": ("liver", "cancer cell line", "hepatocellular carcinoma"),
    "GSM4832654": ("blood", "primary cell", "normal"),
    "GSM4293622": ("embryo", "embryonic stem cell", "normal"),
    "GSM2082526": ("embryo", "embryonic stem cell", "normal"),
    "GSM4504140": ("prostate", "patient-derived xenograft", "prostate adenocarcinoma"),
    "GSM4192431": ("brain", "cancer cell line", "glioblastoma (Grade IV)"),
    "GSM3490787": ("heart", "primary cell", "normal"),
    "GSM4263928": ("breast", "circulating tumor cell line", "breast cancer (metastatic)"),
    "GSM1585244": ("skin", "primary cell", "normal"),
    "GSM4736516": ("brain", "iPSC-derived", "normal"),
    "GSM3720427": ("ovary", "cancer cell line", "ovarian carcinoma (cisplatin-resistant)"),
    "GSM3739081": ("prostate", "primary tissue", "normal"),
    "GSM3753507": ("kidney", "cancer cell line", "normal (transformed)"),
    "GSM1047586": ("skin", "immortalized primary cell", "normal"),
    "GSM4736518": ("brain", "iPSC-derived", "normal"),
    "GSM2714737": ("lung", "cancer cell line", "lung adenocarcinoma"),
    "GSM1608276": ("breast", "normal epithelial cell line", "normal"),
}

QC_PASSTHROUGH = [
    "Start length", "End length", "Periodicity score", "Periodicity distr",
    "CDS coverage (15-40)", "CDS coverage (27-30)", "CDS coverage (dynamic)",
    "Read counts (15-40)", "Read counts (27-30)", "Read counts (dynamic)",
]
COLUMNS = (["cell_line", "ribo_GSM", "ribo_GSE", "rna_GSM", "periodicity", "cds_cov",
            "reads", "max_consec", "norm_depth"] + QC_PASSTHROUGH
           + ["tissue", "cell_type", "disease", "cell_info_GEO"])
INTEGER_COLUMNS = ["reads", "max_consec", "Start length", "End length",
                   "Read counts (15-40)", "Read counts (27-30)", "Read counts (dynamic)"]
ROUNDED_COLUMNS = ["periodicity", "cds_cov", "norm_depth"]

def load_qc(rda=None, qc_csv=None, dump_to=None):
    """The ribobaser QC table: either a previous CSV dump, or `Rscript` over the .rda.

    The .rda is an R binary; `Rscript`'s `write.csv` is how the original selection read it,
    and reproducing its 15-significant-digit output is part of reproducing the published
    floats. `--qc-csv` exists so a machine without R can still run this from a dump.
    """
    if qc_csv:
        frame = pd.read_csv(qc_csv)
    else:
        if not rda:
            raise SystemExit(
                "no QC table: pass --rda (or set RIBOBASER_RDA) or --qc-csv.\n"
                "The .rda is third-party MIT data from CenikLab/ribobaser and is not "
                "redistributed here; set RIBOBASER_RDA to your own copy.")
        rda = Path(rda).resolve()
        if not rda.exists():
            raise SystemExit("no such .rda: %s" % rda)
        if not _which("Rscript"):
            raise SystemExit(
                "Rscript is not on PATH, and it is what reads the .rda.\nEither install R, "
                "or dump the table elsewhere and pass --qc-csv:\n"
                "  Rscript -e 'load(\"%s\"); write.csv(Ribobase_QC_dedup_data, \"qc.csv\", "
                "row.names=FALSE)'" % rda)
        handle, path = tempfile.mkstemp(suffix=".csv")
        os.close(handle)
        try:
            script = ('load("%s"); write.csv(Ribobase_QC_dedup_data, "%s", row.names=FALSE)'
                      % (rda, path))
            result = subprocess.run(["Rscript", "-e", script], capture_output=True)
            if result.returncode != 0:
                raise SystemExit("Rscript failed reading %s:\n%s"
                                 % (rda, result.stderr.decode("utf-8", "replace")))
            frame = pd.read_csv(path)
            if dump_to:
                Path(dump_to).parent.mkdir(parents=True, exist_ok=True)
                Path(dump_to).write_bytes(Path(path).read_bytes())
        finally:
            os.unlink(path)

    required = ["Experiment", "Study", "Cell line", "Species"] + QC_PASSTHROUGH
    missing = [c for c in required if c not in frame.columns]
    if missing:
        raise SystemExit("the QC table is missing %d column(s): %s\nPresent: %s"
                         % (len(missing), ", ".join(missing), ", ".join(frame.columns)))
    if frame["Experiment"].duplicated().any():
        raise SystemExit("the QC table has duplicate Experiment ids; the join would be "
                         "ambiguous")
    return frame

def load_matched_rna(xlsx):
    """GSM -> matched RNA-seq GSM, plus `cell_info_GEO`, from the supplement."""
    sheet = pd.read_excel(xlsx, sheet_name=SHEET, header=SHEET_HEADER_ROW)
    for column in ("experiment_alias", "matched_RNA-seq_experiment_alias", "cell_info_GEO"):
        if column not in sheet.columns:
            raise SystemExit("%s sheet %r has no %r column" % (xlsx, SHEET, column))
    if sheet["experiment_alias"].duplicated().any():
        raise SystemExit("%s sheet %r has duplicate experiment_alias values" % (xlsx, SHEET))
    return sheet.set_index("experiment_alias")

def _which(program):
    for directory in os.environ.get("PATH", "").split(os.pathsep):
        candidate = Path(directory) / program
        if candidate.exists() and os.access(str(candidate), os.X_OK):
            return str(candidate)
    return None

def parse_distr(value):
    """`"[0.68812, 0.6241, …]"` -> a list of floats, exactly as the original parsed it."""
    return [float(x) for x in str(value).strip("[]").split(", ")]

def max_consecutive_above(value, threshold):
    """Longest run of consecutive read lengths whose periodicity is >= `threshold`."""
    longest = run = 0
    for periodicity in parse_distr(value):
        run = run + 1 if periodicity >= threshold else 0
        longest = max(longest, run)
    return longest

def build_pool(qc, metadata):
    """Human samples with matched RNA-seq that clear both hard filters, scored."""
    pool = qc[(qc["Species"] == "human")
              & (qc["Read counts (dynamic)"] >= MIN_READS)].copy()
    pool["rna_GSM"] = pool["Experiment"].map(
        metadata["matched_RNA-seq_experiment_alias"])
    pool = pool[pool["rna_GSM"].notna() & (pool["rna_GSM"].astype(str) != "nan")]

    pool["max_consec"] = pool["Periodicity distr"].apply(
        lambda s: max_consecutive_above(s, DISTR_THR))
    pool = pool[pool["max_consec"] >= DISTR_N].copy()
    if pool.empty:
        raise SystemExit("the filters left no samples; check the QC table and the xlsx")

    depth = np.log10(pool["Read counts (dynamic)"])
    pool["norm_depth"] = (depth - depth.min()) / (depth.max() - depth.min())
    spread = pool["max_consec"].max() - pool["max_consec"].min()
    pool["norm_distr"] = 0.0 if spread == 0 else \
        (pool["max_consec"] - pool["max_consec"].min()) / spread
    pool["score"] = W_DEP * pool["norm_depth"] + W_DISTR * pool["norm_distr"]

    tied = pool["score"].duplicated().sum()
    if tied:
        raise SystemExit("%d tied composite score(s) in the pool: the selection would not "
                         "be deterministic. Add an explicit tie-break before trusting this."
                         % tied)
    return pool

def select_panel(pool):
    """Best-scoring sample per cell line, ordered by score."""
    ordered = pool.sort_values("score", ascending=False)
    best = (ordered.groupby("Cell line", sort=False).first().reset_index()
            .sort_values("score", ascending=False).reset_index(drop=True))
    return best

def _row(record, rna_gsm, norm_depth, round_qc=False):
    """One output row. `round_qc` reproduces the added row's 3-dp QC values."""
    def qc_value(column):
        value = record[column]
        if round_qc and isinstance(value, float):
            return round(value, OVERRIDE_ROUND)
        return value

    gsm = record["Experiment"]
    tissue, cell_type, disease = CURATED_ANNOTATION[gsm]
    row = {
        "cell_line": record["Cell line"],
        "ribo_GSM": gsm,
        "ribo_GSE": record["Study"],
        "rna_GSM": rna_gsm,
        "periodicity": round(float(record["Periodicity score"]), 3),
        "cds_cov": round(float(record["CDS coverage (dynamic)"]), 3),
        "reads": int(record["Read counts (dynamic)"]),
        "max_consec": int(record["max_consec"]),
        "norm_depth": round(float(norm_depth), 3),
        "tissue": tissue,
        "cell_type": cell_type,
        "disease": disease,
        "cell_info_GEO": record["cell_info_GEO"],
    }
    for column in QC_PASSTHROUGH:
        row[column] = qc_value(column)
    return row

def _override_row(qc, metadata, panel):
    """The MCF10A row: everything derived, including its panel-range `norm_depth`.

    Its depth sits BELOW the panel's minimum, so this normalization is negative. That is
    the published value (-0.15) and it is why the row demonstrably post-dates the panel.
    """
    record = qc[qc["Experiment"] == OVERRIDE_ADD]
    if record.empty:
        raise SystemExit("the override sample %s is not in the QC table" % OVERRIDE_ADD)
    record = record.iloc[0].to_dict()
    record["max_consec"] = max_consecutive_above(record["Periodicity distr"], DISTR_THR)
    record["cell_info_GEO"] = metadata.loc[OVERRIDE_ADD, "cell_info_GEO"]

    depth = np.log10(panel["Read counts (dynamic)"].astype(float))
    norm_depth = ((np.log10(float(record["Read counts (dynamic)"])) - depth.min())
                  / (depth.max() - depth.min()))
    rna_gsm = metadata.loc[OVERRIDE_ADD, "matched_RNA-seq_experiment_alias"]
    return _row(record, rna_gsm, norm_depth, round_qc=True)

def build_table(qc, metadata, apply_override=True):
    """The published 24 x 23 table."""
    pool = build_pool(qc, metadata)
    best = select_panel(pool)
    print("[s1] pool: %d samples over %d cell lines"
          % (len(pool), pool["Cell line"].nunique()))

    selected, dropped_line = best.copy(), None
    if apply_override:
        if OVERRIDE_DROP not in set(selected["Experiment"]):
            raise SystemExit(
                "the override expects to drop %s, but the automatic selection did not "
                "choose it. The inputs differ from the ones the published table was built "
                "from; re-derive the override before trusting the result." % OVERRIDE_DROP)
        dropped_line = selected.loc[selected["Experiment"] == OVERRIDE_DROP,
                                    "Cell line"].iloc[0]
        selected = selected[selected["Experiment"] != OVERRIDE_DROP]

    rows = []
    for record in selected.to_dict("records"):
        record["cell_info_GEO"] = metadata.loc[record["Experiment"], "cell_info_GEO"]
        rows.append(_row(record, record["rna_GSM"], record["norm_depth"]))
    if apply_override:
        added = _override_row(qc, metadata, selected)
        rows.append(added)
        print("[s1] quality override: -%s (%s)  +%s (%s)"
              % (OVERRIDE_DROP, dropped_line, OVERRIDE_ADD, added["cell_line"]))

    table = pd.DataFrame(rows)[COLUMNS]
    for column in INTEGER_COLUMNS:
        table[column] = table[column].astype("int64")
    for column in ROUNDED_COLUMNS:
        table[column] = table[column].astype("float64")

    unknown = set(table["ribo_GSM"]) - set(CURATED_ANNOTATION)
    if unknown:
        raise SystemExit("no curated tissue/cell_type/disease for: %s"
                         % ", ".join(sorted(unknown)))
    if table["ribo_GSM"].duplicated().any() or table["rna_GSM"].duplicated().any():
        raise SystemExit("the panel has duplicate ribo or RNA accessions")
    return table

def sha256(path):
    digest = hashlib.sha256()
    with open(path, "rb") as handle:
        for block in iter(lambda: handle.read(1 << 20), b""):
            digest.update(block)
    return digest.hexdigest()

def main(argv=None):
    parser = argparse.ArgumentParser(
        description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument("--rda", default=os.environ.get("RIBOBASER_RDA"),
                        help="Ribobase_QC_dedup_data.rda (default: $RIBOBASER_RDA)")
    parser.add_argument("--qc-csv", help="a previously dumped QC table, instead of --rda")
    parser.add_argument("--xlsx", required=True, help="41587_2025_2718_MOESM3_ESM.xlsx")
    parser.add_argument("--output", type=Path, required=True,
                        help="where to write the regenerated table. REQUIRED: give a "
                             "temporary path to inspect it, or "
                             "supporting_information/S1_Table/samples.csv to replace the "
                             "shipped copy.")
    parser.add_argument("--dump-qc", type=Path,
                        help="also save the Rscript dump of the QC table here")
    parser.add_argument("--no-override", action="store_true",
                        help="the automatic selection alone, WITHOUT the MCF10A quality "
                             "override -- 24 rows that are not the published panel")
    parser.add_argument("--verify", nargs="?", const=str(PUBLISHED), default=None,
                        help="compare the result byte-for-byte against the published CSV "
                             "(default: %s)" % PUBLISHED.relative_to(REPO))
    parser.add_argument("--force", action="store_true", help="overwrite --output")
    args = parser.parse_args(argv)

    qc = load_qc(args.rda, args.qc_csv, args.dump_qc)
    metadata = load_matched_rna(args.xlsx)
    table = build_table(qc, metadata, apply_override=not args.no_override)

    if args.output.exists() and not args.force:
        raise SystemExit("refusing to overwrite %s -- pass --force" % args.output)
    args.output.parent.mkdir(parents=True, exist_ok=True)
    table.to_csv(args.output, index=False)
    print("[s1] wrote %s  (%d rows x %d columns, sha256 %s)"
          % (args.output, len(table), len(table.columns), sha256(args.output)[:16]))

    if args.verify:
        reference = Path(args.verify)
        if not reference.exists():
            raise SystemExit("no reference to verify against: %s" % reference)
        produced, expected = args.output.read_bytes(), reference.read_bytes()
        if produced == expected:
            print("[s1] VERIFIED byte-identical to %s" % reference)
            if sha256(reference) != PUBLISHED_SHA256 and reference == PUBLISHED:
                print("[s1] note: the published copy's checksum has changed since this "
                      "script was written")
            return 0
        print("[s1] DIFFERS from %s" % reference, file=sys.stderr)
        _report_difference(table, reference)
        return 1
    return 0

def _report_difference(table, reference):
    """Name what differs, rather than leaving a reader to diff 24 rows by hand."""
    expected = pd.read_csv(reference, dtype=str)
    produced = table.astype(str)
    if list(produced.columns) != list(expected.columns):
        print("  columns differ:\n    produced: %s\n    expected: %s"
              % (list(produced.columns), list(expected.columns)), file=sys.stderr)
        return
    if len(produced) != len(expected):
        print("  row count: produced %d, expected %d" % (len(produced), len(expected)),
              file=sys.stderr)
    produced_gsms, expected_gsms = list(produced["ribo_GSM"]), list(expected["ribo_GSM"])
    if produced_gsms != expected_gsms:
        only_produced = [g for g in produced_gsms if g not in expected_gsms]
        only_expected = [g for g in expected_gsms if g not in produced_gsms]
        print("  membership/order: produced-only %s, expected-only %s"
              % (only_produced or "-", only_expected or "-"), file=sys.stderr)
    shared = [g for g in produced_gsms if g in expected_gsms]
    produced_i, expected_i = produced.set_index("ribo_GSM"), expected.set_index("ribo_GSM")
    for gsm in shared:
        for column in produced.columns:
            if column == "ribo_GSM":
                continue
            a, b = produced_i.loc[gsm, column], expected_i.loc[gsm, column]
            if a != b:
                print("  %s %-24s produced=%r expected=%r" % (gsm, column, a, b),
                      file=sys.stderr)

if __name__ == "__main__":
    sys.exit(main())
