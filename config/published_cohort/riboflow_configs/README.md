# RiboFlow_v2 configurations for the 24 published samples

One YAML per cell line, plus [`../run_riboflow_cohort.sh`](../run_riboflow_cohort.sh): the
configuration used to produce the alignments behind Figures 2–6. The alignments were
produced with RiboFlow_v2 at commit
[`e5e041c6`](https://github.com/ribosomeprofiling/riboflow/commit/e5e041c6fa842c27fabe46d2ca87d8aff3696874).
RiboFlow_v2 is not included in this repository. The Table 1 benchmark is a separate, later
run whose records are in [`benchmark/runs/`](../../../benchmark/runs/).

## Sanitization

Machine-specific absolute paths were replaced with placeholders (213 substitutions, all
inside absolute paths). A private IP address was removed from a header comment of the run
script, and its `fastqs_ready()` check was retargeted at `${FASTQ_DIR}`.

| Placeholder | Meaning | Occurrences |
|---|---|---|
| `${REFERENCES_FOR_RIBOFLOW}` | a clone of `ribosomeprofiling/references_for_riboflow` (rRNA filter, APPRIS transcriptome, regions BED, transcript-lengths TSV) | 96 |
| `${STAR_INDEX_DIR}` | the STAR GRCh38 genome index; some YAMLs append `/GRCh38` | 48 |
| `${FASTQ_DIR}` | input FASTQ root, laid out `$FASTQ_DIR/{ribo,rna}/<GSM>/<SRR>_1.fastq.gz` | 36 ribo + 33 RNA |
| `${RIBOFLOW_REPO}` | the RiboFlow_v2 checkout containing `main.nf` | run script |
| `${NXF_PROFILE}` | optional Nextflow profile; the published run used the site profile `ls6` | run script |

Preserved verbatim: every GSM and SRR accession; per-sample `clip_arguments` including the
3′ adapter; per-sample `ribo.read_length.min`/`max`; `left_span: 35`, `right_span: 10`,
`metagene_radius: 50`; every STAR `ribo_arguments` and `rnaseq_arguments` flag; the
pipeline's write-time filters (`mapping_quality_cutoff: 0` genome, `10` transcriptome,
`ribo_filter_flags: 2052`); `dedup_method: position`; and the repository-relative
`intermediates/<sample>` and `output/<sample>` bases.

The write-time cutoffs are distinct from the analysis filters this repository applies when
reading the BAMs (genome MAPQ > 4, transcriptome MAPQ ≥ 42).

Every YAML sets `ribo.ref_name: appris-v1` while its reference paths point at
`appris_human_v2_*`. `ref_name` is a label written inside the `.ribo` file, not a data path;
the annotation used is APPRIS human v2.

## Per-sample parameters

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
| megakaryocytes | 25–30 | 6 | 3 |
| neurons | 23–29 | 1 | 1 |
| normal_prostate | 25–30 | 1 | 1 |
| RD-CCL-136 | 24–29 | 3 | 3 |
| U-251 | 25–30 | 1 | 1 |
| U-343 | 24–30 | 1 | 1 |
| U2OS | 25–31 | 1 | 1 |

## Running

```bash
export RIBOFLOW_REPO=/path/to/RiboFlow_v2
export FASTQ_DIR=/path/to/fastqs
export REFERENCES_FOR_RIBOFLOW=/path/to/references_for_riboflow
export STAR_INDEX_DIR=/path/to/star_index
export NXF_PROFILE=              # unset for a local run

cp config/published_cohort/riboflow_configs/*.yaml "$RIBOFLOW_REPO/RiboFlow_YAMLs_cohort/"
bash config/published_cohort/run_riboflow_cohort.sh
```

Finished samples are skipped and `-resume` continues partial ones; `FORCE=1` re-submits
every sample using the Nextflow work cache.

## Licensing

The run script is covered by this repository's [`LICENSE`](../../../LICENSE). The YAMLs are
configuration; the reference paths they name point at third-party data not redistributed here.
