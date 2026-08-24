# RiboFlow configurations — the 24 published samples

One YAML per cell line, plus [`../run_riboflow_cohort.sh`](../run_riboflow_cohort.sh): the
exact configuration used to produce the alignments behind Figures 2–5.

## These are sanitized copies

Machine-specific absolute paths were replaced with placeholders — **213 substitutions, every
one inside an absolute path**. Nothing scientifically load-bearing was touched.

| Placeholder | Set to | Occurrences |
|---|---|---|
| `${REFERENCES_FOR_RIBOFLOW}` | a clone of `ribosomeprofiling/references_for_riboflow` — the rRNA filter, the APPRIS transcriptome, the regions BED and the transcript-lengths TSV | 96 |
| `${STAR_INDEX_DIR}` | the STAR GRCh38 genome index; some YAMLs append `/GRCh38` | 48 |
| `${FASTQ_DIR}` | input FASTQ root, laid out `$FASTQ_DIR/{ribo,rna}/<GSM>/<SRR>_1.fastq.gz` | 36 ribo + 33 RNA |
| `${RIBOFLOW_REPO}` | the `RiboFlow_v2` checkout containing `main.nf` | run script |
| `${NXF_PROFILE}` | optional Nextflow profile; defaults to `ls6`, which is **site-specific** | run script |

Two further changes to the run script: a private IP address was removed from a header
comment, and its `fastqs_ready()` pre-flight check was retargeted at `${FASTQ_DIR}` — it greps
each YAML for FASTQ paths, so without that change the sanitized script would have found zero
FASTQs and silently skipped every sample.

**Preserved verbatim:** every GSM and SRR accession; per-sample `clip_arguments` including the
3′ adapter; per-sample `ribo.read_length.min`/`max`; `left_span: 35`, `right_span: 10`,
`metagene_radius: 50`; every STAR `ribo_arguments` and `rnaseq_arguments` flag; the pipeline's
write-time filters (`mapping_quality_cutoff: 0` genome, `10` transcriptome,
`ribo_filter_flags: 2052`); `dedup_method: position`; and the repository-relative
`intermediates/<sample>` and `output/<sample>` bases.

Those write-time cutoffs are a **different layer** from the analysis filters this repository
applies when reading the BAMs (genome `MAPQ > 4`, transcriptome `MAPQ >= 42` — see
the QC tables under `data/ribo_seq_qc/`). Both are real; neither is the other.

## `ref_name: appris-v1` is not a version conflict

Every YAML sets `ribo.ref_name: appris-v1` while its reference paths point at
`appris_human_v2_*`. `ref_name` is a **label written inside the `.ribo` file**, not a data
path — the annotation actually used is APPRIS human **v2**. Preserved verbatim and flagged
here so it is not misread as an inconsistency or "fixed" by a future reader.

## Per-sample parameters

The read-length window is per sample, and several samples have multiple runs:

| Sample | Read length | Ribo FASTQ | RNA FASTQ |
|---|---|---|---|
| A2780 | 26–29 | 1 | 1 |
| A549 | 22–29 | 1 | 1 |
| BJ | 26–31 | 1 | 1 |
| BRx-142 | 26–31 | 4 | 4 |
| cardiac_fibroblasts | 24–29 | 1 | 1 |
| Cybrid_Cells | 23–30 | 1 | 1 |
| early_neurons | 23–29 | 2 | 2 |
| fibroblast | 23–30 | 1 | 1 |
| H1-hESC | 25–29 | 1 | 1 |
| H9-hESC | 25–29 | 2 | 2 |
| HEK293 | 26–31 | 1 | 1 |
| HEK293T | 25–29 | 1 | 1 |
| HeLa | 26–29 | 1 | 1 |
| hESC | 23–29 | 1 | 1 |
| Huh7 | 24–29 | 1 | 1 |
| LuCaP-PDX | 26–30 | 1 | 1 |
| MCF10A | 26–28 | 1 | 1 |
| megakaryocytes | 25–30 | **6** | 3 |
| neurons | 23–29 | 1 | 1 |
| normal_prostate | 25–30 | 1 | 1 |
| RD-CCL-136 | 24–29 | 3 | 3 |
| U-251 | 25–30 | 1 | 1 |
| U-343 | 24–30 | 1 | 1 |
| U2OS | 25–31 | 1 | 1 |

`megakaryocytes` is the sample that behaves as a data-quality outlier throughout the analysis;
its 6 ribo runs against 3 RNA runs is one visible reason.

## Running it

```bash
export RIBOFLOW_REPO=/path/to/RiboFlow_v2
export FASTQ_DIR=/path/to/fastqs
export REFERENCES_FOR_RIBOFLOW=/path/to/references_for_riboflow
export STAR_INDEX_DIR=/path/to/star_index
export NXF_PROFILE=              # unset for a local run; the published run used the site profile "ls6"

cp config/published_cohort/riboflow_configs/*.yaml "$RIBOFLOW_REPO/RiboFlow_YAMLs_cohort/"
bash config/published_cohort/run_riboflow_cohort.sh
```

The script is idempotent: finished samples are skipped and `-resume` lets partial ones
continue. `FORCE=1` re-submits every sample, relying on the Nextflow work cache so alignment
is not repeated.

## RiboFlow is not vendored, and its version for these runs is unknown

RiboFlow is referenced, not included.
The deployment that produced the Figure 2–5 alignments is **not identified**: the in-tree copy
reported `version = '0.0.0'`, inconsistent with the DSL2 v26.04.2 pipeline the Methods
describe. This remains unresolved.

A commit *is* pinned for the Table 1 benchmark — a different, later run. Do not read it as the
version behind these configurations; the benchmark's own run records are in
[`benchmark/runs/`](../../../benchmark/runs/).

## Licensing

The run script is author-written code under this repository's [`LICENSE`](../../../LICENSE).
The YAMLs are configuration; the reference *paths* they name point at third-party data not
redistributed here.
