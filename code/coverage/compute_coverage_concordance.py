#!/usr/bin/env python3
"""Genome-versus-transcriptome concordance tables, computed from coverage HDF5 files."""
from __future__ import annotations

import argparse
import gzip
import sys
from pathlib import Path

import numpy as np
import pandas as pd

HERE = Path(__file__).resolve().parent

PER_SAMPLE_COLUMNS = [
    "sample", "n_regions_total", "n_regions_covered_either", "n_bases_inter",
    "n_bases_union", "n_bases_all", "n_bases_g_only", "n_bases_t_only",
    "spearman_inter", "pearson_inter", "spearman_union", "pearson_union",
    "spearman_all", "pearson_all", "frame_agree", "g_psite_in_cds", "t_psite_in_cds"]
PER_TRANSCRIPT_COLUMNS = [
    "sample", "transcript_id", "gene_name", "cds_len_interior", "n_cov_g", "n_cov_t",
    "n_inter", "spearman", "pearson"]

SIGNAL_SETS = {
    "psite": {
        "signals": ("genome_psite", "txome_psite"),
        "per_sample": "region_concordance_per_sample.tsv",
        "per_transcript": "region_concordance_per_transcript.tsv",
        "with_frame": True,
        "totals_trimmed": False,      # the P-site count is over the UNtrimmed CDS
    },
    "footprint": {
        "signals": ("genome_footprint", "txome_footprint"),
        "per_sample": "region_coverage_per_sample.tsv",
        "per_transcript": "region_coverage_per_transcript.tsv",
        "with_frame": False,
        "totals_trimmed": True,
    },
}

def log(message):
    print("[concordance] %s" % message, flush=True)

def _spear(genome, txome):
    from scipy.stats import spearmanr
    if len(genome) > 2 and np.std(genome) > 0 and np.std(txome) > 0:
        return float(spearmanr(genome, txome).correlation)
    return np.nan

def _pe_log2(genome, txome, pseudocount=1.0):
    """Pearson r on log2(count + pseudocount), hand-centred in float64.

    Never delegate to `scipy.stats.pearsonr`: it differs at ~1e-15 and clamps pooled arrays to 1.0.
    """
    if len(genome) <= 2:
        return np.nan
    x = np.log2(genome.astype(np.float64) + pseudocount)
    y = np.log2(txome.astype(np.float64) + pseudocount)
    x -= x.mean()
    y -= y.mean()
    sxx = float(np.dot(x, x))
    syy = float(np.dot(y, y))
    if sxx <= 0 or syy <= 0:
        return np.nan
    return max(-1.0, min(1.0, float(np.dot(x, y) / np.sqrt(sxx * syy))))

def _codon_view(vector):
    usable = (len(vector) // 3) * 3
    return vector[:usable].reshape(-1, 3)

def load_sample(coverage_path, kind):
    """Interior slices, per-transcript metadata, and the pooled arrays for one sample."""
    sys.path.insert(0, str(HERE))
    import coverage_schema

    spec = SIGNAL_SETS[kind]
    wrapper = coverage_schema.open_coverage(coverage_path)
    trim = wrapper.trim
    offsets = wrapper.coverage_offset
    ids = np.array(wrapper.transcript_ids)
    gene_names = np.array(wrapper.gene_names)
    cds_start = wrapper.cds_start.copy()
    cds_end = wrapper.cds_end.copy()
    # A transcript without a CDS reads as an empty window at 0, as before schema 3.
    no_cds = cds_start == coverage_schema.NO_CDS
    cds_start[no_cds] = 0
    cds_end[no_cds] = 0

    genome = wrapper.signal(spec["signals"][0])
    txome = wrapper.signal(spec["signals"][1])
    # The CDS coverage key: P-sites key a transcript on any read in the UNtrimmed CDS;
    # footprints only on a non-zero trimmed interior. Derived here, never stored.
    key_trim = trim if spec["totals_trimmed"] else 0
    key_genome = coverage_schema.window_sums(
        genome, offsets, cds_start + key_trim, cds_end - key_trim) > 0
    key_txome = coverage_schema.window_sums(
        txome, offsets, cds_start + key_trim, cds_end - key_trim) > 0
    sample = wrapper.sample
    wrapper.close()

    interior = (cds_end - cds_start) - 2 * trim
    eligible = interior > 0
    covered = eligible & (key_genome | key_txome)

    lo = offsets + cds_start + trim
    hi = offsets + cds_end - trim
    return {
        "sample": sample, "trim": trim, "ids": ids, "gene_names": gene_names,
        "genome": genome, "txome": txome, "lo": lo, "hi": hi,
        "interior": interior, "eligible": eligible, "covered": covered,
        "key_genome": key_genome, "key_txome": key_txome,
        "cds_start": cds_start, "cds_end": cds_end, "offsets": offsets,
    }

def pooled_row(data, kind):
    """The per-sample row: every interior base of every transcript, concatenated."""
    spec = SIGNAL_SETS[kind]
    covered_index = np.where(data["covered"])[0]
    if covered_index.size:
        genome_parts = [data["genome"][data["lo"][i]:data["hi"][i]] for i in covered_index]
        txome_parts = [data["txome"][data["lo"][i]:data["hi"][i]] for i in covered_index]
        pooled_g = np.concatenate(genome_parts)
        pooled_t = np.concatenate(txome_parts)
        del genome_parts, txome_parts
    else:
        pooled_g = pooled_t = np.zeros(0, dtype=np.int32)

    zero_block = int(data["interior"][data["eligible"] & ~data["covered"]].sum())
    inter = (pooled_g > 0) & (pooled_t > 0)
    union = (pooled_g > 0) | (pooled_t > 0)
    all_g = np.concatenate([pooled_g, np.zeros(zero_block, np.int32)]) if zero_block else pooled_g
    all_t = np.concatenate([pooled_t, np.zeros(zero_block, np.int32)]) if zero_block else pooled_t

    row = {
        "sample": data["sample"],
        "n_regions_total": int(data["eligible"].sum()),
        "n_regions_covered_either": int(covered_index.size),
        "n_bases_inter": int(inter.sum()),
        "n_bases_union": int(union.sum()),
        "n_bases_all": int(all_g.size),
        "n_bases_g_only": int(((pooled_g > 0) & (pooled_t == 0)).sum()),
        "n_bases_t_only": int(((pooled_t > 0) & (pooled_g == 0)).sum()),
        "spearman_inter": _spear(pooled_g[inter], pooled_t[inter]),
        "pearson_inter": _pe_log2(pooled_g[inter], pooled_t[inter]),
        "spearman_union": _spear(pooled_g[union], pooled_t[union]),
        "pearson_union": _pe_log2(pooled_g[union], pooled_t[union]),
        "spearman_all": _spear(all_g, all_t),
        "pearson_all": _pe_log2(all_g, all_t),
    }
    if spec["with_frame"]:
        both = np.where(data["key_genome"] & data["key_txome"] & data["eligible"])[0]
        agree = total = 0
        for i in both:
            g = _codon_view(data["genome"][data["lo"][i]:data["hi"][i]])
            t = _codon_view(data["txome"][data["lo"][i]:data["hi"][i]])
            covered_codons = (g.sum(1) > 0) & (t.sum(1) > 0)
            if covered_codons.any():
                agree += int((g[covered_codons].argmax(1) == t[covered_codons].argmax(1)).sum())
                total += int(covered_codons.sum())
        row["frame_agree"] = (agree / total) if total else np.nan

    # the diagnostic in-CDS totals: untrimmed for P-sites, trimmed interior for footprints
    if spec["totals_trimmed"]:
        lo, hi = data["lo"], data["hi"]
    else:
        lo = data["offsets"] + data["cds_start"]
        hi = data["offsets"] + data["cds_end"]
    for signal, column in (("genome", "g_psite_in_cds"), ("txome", "t_psite_in_cds")):
        values = data[signal]
        prefix = np.concatenate([[0], np.cumsum(values, dtype=np.int64)])
        row[column] = int((prefix[np.maximum(hi, lo)] - prefix[lo]).sum())
    return row

def transcript_rows(data):
    """One row per transcript covered on BOTH routes, in sorted transcript_id order."""
    both = np.where(data["key_genome"] & data["key_txome"] & data["eligible"])[0]
    rows = []
    for i in both:
        g = data["genome"][data["lo"][i]:data["hi"][i]]
        t = data["txome"][data["lo"][i]:data["hi"][i]]
        rows.append({
            "sample": data["sample"],
            "transcript_id": data["ids"][i],
            "gene_name": data["gene_names"][i],
            "cds_len_interior": int(len(g)),
            "n_cov_g": int((g > 0).sum()),
            "n_cov_t": int((t > 0).sum()),
            "n_inter": int(((g > 0) & (t > 0)).sum()),
            "spearman": _spear(g, t),
            "pearson": _pe_log2(g, t),
        })
    return rows

def write_tsv(frame, path, gzip_it=False):
    path.parent.mkdir(parents=True, exist_ok=True)
    text = frame.to_csv(sep="\t", index=False, lineterminator="\n")
    if gzip_it:
        with gzip.GzipFile(str(path), "wb", mtime=0) as handle:
            handle.write(text.encode())
    else:
        path.write_text(text)
    return path

def main(argv=None):
    parser = argparse.ArgumentParser(
        description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument("--coverage", required=True, type=Path,
                        help="directory of <sample>.shared_coverage.h5 files")
    parser.add_argument("--samples", help="comma-separated subset")
    parser.add_argument("--output", type=Path, default=Path("results/coverage/concordance"))
    parser.add_argument("--gzip-per-transcript", action="store_true", default=True)
    args = parser.parse_args(argv)

    files = sorted(args.coverage.glob("*.shared_coverage.h5"))
    if args.samples:
        wanted = {s.strip() for s in args.samples.split(",") if s.strip()}
        files = [f for f in files if f.name.split(".")[0] in wanted]
    if not files:
        raise SystemExit("no coverage files found in %s" % args.coverage)
    log("%d coverage file(s)" % len(files))

    results = {}
    for kind, spec in SIGNAL_SETS.items():
        per_sample, per_transcript = [], []
        for path in files:
            data = load_sample(path, kind)
            per_sample.append(pooled_row(data, kind))
            per_transcript.extend(transcript_rows(data))
            log("  %-22s %-10s %d transcript rows"
                % (data["sample"], kind, len(per_transcript)))
            del data
        sample_frame = pd.DataFrame(per_sample).sort_values("sample")
        sample_frame = sample_frame.reindex(
            columns=[c for c in PER_SAMPLE_COLUMNS if c in sample_frame.columns])
        transcript_frame = pd.DataFrame(per_transcript).sort_values(
            ["sample", "transcript_id"])[PER_TRANSCRIPT_COLUMNS]

        sample_path = write_tsv(sample_frame, args.output / spec["per_sample"])
        transcript_path = write_tsv(
            transcript_frame,
            args.output / (spec["per_transcript"] + (".gz" if args.gzip_per_transcript else "")),
            gzip_it=args.gzip_per_transcript)
        log("wrote %s and %s" % (sample_path.name, transcript_path.name))
        results[kind] = {"per_sample": str(sample_path),
                         "per_transcript": str(transcript_path)}

    return 0

if __name__ == "__main__":
    sys.path.insert(0, str(HERE))
    sys.exit(main())
