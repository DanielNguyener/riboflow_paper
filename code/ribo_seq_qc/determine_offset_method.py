#!/usr/bin/env python3
"""Validate the shipped P-site offsets, and record which detector branch produced each."""
from __future__ import annotations

import argparse
import json
import re
import subprocess
import sys
from collections import defaultdict
from concurrent.futures import ThreadPoolExecutor, as_completed
from pathlib import Path

import numpy as np
import pandas as pd

_HERE = Path(__file__).resolve().parent
_CODE = _HERE.parent
for _d in (_CODE / "common" / "ribo_seq_qc", _CODE / "common"):
    if str(_d) not in sys.path:
        sys.path.insert(0, str(_d))
import config
from psite_offset import ribotish_get_offset, get_offset_periodicity
import bam_inputs as fc
if str(_HERE) not in sys.path:
    sys.path.insert(0, str(_HERE))
import qc_core

PRE_WIN_UP, PRE_WIN_DN = 50, 30
WIN_CODONS, MIN_DOWN, DOM_FRAC = 10, 200, 0.40
WIN = WIN_CODONS * 3
_CDS_RE = re.compile(r"\|CDS:(\d+)-(\d+)\|")
QC_FOR = {"genome": fc.genome_qc_path, "transcriptome": fc.txome_qc_path}

def default_out_tsv():
    return fc.output_root() / "ribo_seq_qc" / "offsets" / "offset_method_per_length.tsv"

_CDS_CORE = None

def _genome_cds_core():
    """RiboPy's CDS core in genomic coordinates, parsed once per process (a GTF parse is costly)."""
    global _CDS_CORE
    if _CDS_CORE is None:
        _CDS_CORE = qc_core.genome_cds_core_intervals()
    return _CDS_CORE

def _phase1(cds_length_counts):
    """The selected window, from `qc_core.select_read_lengths` -- NOT a local copy that could drift."""
    return qc_core.select_read_lengths(cds_length_counts)[0]

def genome_metagene(sample):
    import pysam, pyranges as pr
    rec = {"Chromosome": [], "pos5": [], "Strand": [], "length": []}
    bam = pysam.AlignmentFile(str(fc.genome_bam(sample)), "rb")
    for r in bam.fetch():
        if r.is_unmapped or r.is_secondary or r.is_supplementary:
            continue
        if not fc.is_unique_genome_read(r):       # NH == 1
            continue
        rlen = r.query_length
        if not (config.MIN_LEN <= rlen <= config.MAX_LEN):
            continue
        if r.is_reverse:
            p5, strand = r.reference_end - 1, "-"
        else:
            p5, strand = r.reference_start, "+"
        rec["Chromosome"].append(r.reference_name); rec["pos5"].append(p5)
        rec["Strand"].append(strand); rec["length"].append(rlen)
    bam.close()
    reads = pd.DataFrame(rec)
    phase1 = _phase1(qc_core.cds_length_hist_genome(
        reads, _genome_cds_core(), qc_core.SELECT_MIN_LEN, qc_core.SELECT_MAX_LEN))

    ann = config.load_annotation()
    tx = (ann.groupby("transcript_id", sort=False)
          .agg(Chromosome=("Chromosome", "first"), Strand=("Strand", "first"),
               cds_genomic_start=("cds_genomic_start", "first")).reset_index(drop=True))
    txp, txm = tx[tx["Strand"] == "+"].copy(), tx[tx["Strand"] == "-"].copy()
    txp["Start"] = (txp["cds_genomic_start"] - PRE_WIN_UP).clip(lower=0)
    txp["End"] = txp["cds_genomic_start"] + PRE_WIN_DN
    txm["Start"] = (txm["cds_genomic_start"] - PRE_WIN_DN + 1).clip(lower=0)
    txm["End"] = txm["cds_genomic_start"] + PRE_WIN_UP + 1
    windows = pd.concat([txp, txm], ignore_index=True)

    pos5_pr = pr.PyRanges(pd.DataFrame({
        "Chromosome": reads["Chromosome"].values, "Start": reads["pos5"].values,
        "End": reads["pos5"].values + 1, "Strand": reads["Strand"].values,
        "read_idx": reads.index.values}))
    win_pr = pr.PyRanges(windows[["Chromosome", "Start", "End", "Strand", "cds_genomic_start"]])
    joined = pos5_pr.join(win_pr, strandedness="same")
    reads["rel_pos"] = np.nan
    if not joined.df.empty:
        j = joined.df.copy()
        plus = j["Strand"] == "+"
        j["rel_pos"] = np.where(plus, j["Start"] - j["cds_genomic_start"],
                                j["cds_genomic_start"] - j["Start"]).astype(float)
        j = j.drop_duplicates("read_idx", keep="first")
        reads.loc[j["read_idx"].values, "rel_pos"] = j["rel_pos"].values
    return _pre_counts(reads), phase1

def txome_metagene(sample):
    import pysam
    bam = pysam.AlignmentFile(str(fc.txome_bam(sample)), "rb")
    ref_cds0 = {}
    for ref in bam.references:
        m = _CDS_RE.search(ref)
        ref_cds0[ref] = (int(m.group(1)) - 1) if m else None
    rec = {"ref": [], "pos5": [], "length": []}
    for r in bam.fetch(until_eof=True):
        if not fc.is_unique_txome_read(r):        # MAPQ >= 42
            continue
        rlen = r.query_length
        if not (config.MIN_LEN <= rlen <= config.MAX_LEN):
            continue
        rec["ref"].append(r.reference_name); rec["pos5"].append(r.reference_start)
        rec["length"].append(rlen)
    bam.close()
    reads = pd.DataFrame(rec)
    with pysam.AlignmentFile(str(fc.txome_bam(sample)), "rb") as cds_bam:
        phase1 = _phase1(qc_core.cds_length_hist_transcriptome(
            cds_bam, qc_core.SELECT_MIN_LEN, qc_core.SELECT_MAX_LEN))
    cds0 = reads["ref"].map(ref_cds0)
    reads["rel_pos"] = (reads["pos5"] - cds0).astype(float)
    return _pre_counts(reads), phase1

def _pre_counts(reads):
    pre = reads[reads["rel_pos"].between(-PRE_WIN_UP, PRE_WIN_DN - 1)].copy()
    pre["rel_int"] = pre["rel_pos"].astype(int)
    cnt = pre.groupby(["length", "rel_int"]).size()
    out = defaultdict(lambda: defaultdict(int))
    for (rlen, pos), c in cnt.items():
        out[int(rlen)][int(pos)] = int(c)
    return out

def classify(counts):
    upstream = {p: c for p, c in counts.items() if p < 0}
    if not upstream or max(upstream.values()) == 0:
        return "default_no_upstream", 0, float("nan"), -1
    mass = [0, 0, 0]
    for p, c in counts.items():
        if 0 <= p < WIN:
            mass[p % 3] += c
    total = sum(mass)
    if total == 0:
        return "fallback_argmax", 0, float("nan"), -1
    dom_frame = int(max(range(3), key=lambda r: mass[r]))
    dom_frac = max(mass) / total
    method = "periodicity" if (total >= MIN_DOWN and dom_frac >= DOM_FRAC) else "fallback_argmax"
    return method, int(total), float(dom_frac), dom_frame

def worker(sample, alignment):
    build = genome_metagene if alignment == "genome" else txome_metagene
    pre_counts, phase1 = build(sample)
    tbl = pd.read_csv(QC_FOR[alignment]())
    tbl = tbl[tbl["sample"] == sample].set_index("read_length")
    rows = []
    for rlen in phase1:
        counts = {p: pre_counts[rlen].get(p, 0) for p in range(-PRE_WIN_UP, PRE_WIN_DN)}
        method, total, dom_frac, dom_frame = classify(counts)
        off_per = get_offset_periodicity(counts)
        off_arg = ribotish_get_offset(counts)
        off_tbl = int(tbl.loc[rlen, "psite_offset"]) if rlen in tbl.index else None
        rows.append({
            "sample": sample, "alignment": alignment, "read_length": int(rlen),
            "method": method, "n_down_first10": total, "dom_frac": round(dom_frac, 4),
            "dom_frame": dom_frame, "offset_periodicity": off_per, "offset_argmax": off_arg,
            "offset_in_table": off_tbl,
            "reproduces_table": (off_tbl is not None and off_per == off_tbl),
            "periodicity_changed_value": (off_per != off_arg),
        })
    return rows

def run_worker_subproc(sample, alignment):
    cmd = [sys.executable, str(_HERE / "determine_offset_method.py"), "--worker", sample, alignment]
    p = subprocess.run(cmd, capture_output=True, text=True)
    if p.returncode != 0:
        sys.stderr.write(p.stderr)
        raise RuntimeError(f"worker failed {sample}/{alignment} (exit {p.returncode})")
    line = [ln for ln in p.stdout.splitlines() if ln.strip().startswith("[")][-1]
    return json.loads(line)

def main():
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--worker", nargs=2, metavar=("SAMPLE", "ALIGNMENT"))
    ap.add_argument("--workers", type=int, default=10)
    ap.add_argument("--out", default=None)
    args = ap.parse_args()

    if args.worker:
        sample, alignment = args.worker
        print(json.dumps(worker(sample, alignment)))
        return 0

    samples = fc.discover_samples()
    jobs = [(s, a) for a in ("genome", "transcriptome") for s in samples]
    print(f"{len(samples)} samples × 2 alignments = {len(jobs)} jobs", flush=True)
    all_rows = []
    with ThreadPoolExecutor(max_workers=args.workers) as pool:
        futs = {pool.submit(run_worker_subproc, s, a): (s, a) for s, a in jobs}
        for f in as_completed(futs):
            s, a = futs[f]
            all_rows.extend(f.result())
            print(f"  done {s}/{a}", flush=True)

    df = pd.DataFrame(all_rows).sort_values(["alignment", "sample", "read_length"])
    out = Path(args.out) if args.out else default_out_tsv()
    out.parent.mkdir(parents=True, exist_ok=True)
    df.to_csv(out, sep="\t", index=False)
    print(f"\nwrote {out}  ({len(df)} cells)", flush=True)

    n_bad = int((~df["reproduces_table"]).sum())
    print("\n===== validation =====", flush=True)
    if n_bad:
        print(f"  [FAIL] {n_bad} cells did NOT reproduce the table offset — metagene drift; "
              f"method labels NOT trustworthy:", flush=True)
        print(df[~df["reproduces_table"]].to_string(index=False), flush=True)
        return 1
    print(f"  [PASS] all {len(df)} cells reproduce the stored table offset "
          f"(rebuilt metagene == pipeline's)", flush=True)

    print("\n===== method counts (per alignment) =====", flush=True)
    print(df.groupby(["alignment", "method"]).size().to_string(), flush=True)

    fb = df[df["method"] != "periodicity"]
    print(f"\n===== fallback / non-periodicity cells ({len(fb)}) =====", flush=True)
    if len(fb):
        print(fb[["sample", "alignment", "read_length", "method", "n_down_first10",
                  "dom_frac", "offset_in_table", "offset_argmax",
                  "periodicity_changed_value"]].to_string(index=False), flush=True)
    else:
        print("  none — every phase-1 cell passed the periodicity gate", flush=True)

    changed = df[(df["method"] == "periodicity") & (df["periodicity_changed_value"])]
    print(f"\nperiodicity ran AND changed the offset vs argmax in {len(changed)} cells", flush=True)
    if len(changed):
        print(changed[["sample", "alignment", "read_length", "offset_argmax",
                       "offset_in_table"]].to_string(index=False), flush=True)
    return 0

if __name__ == "__main__":
    raise SystemExit(main())
