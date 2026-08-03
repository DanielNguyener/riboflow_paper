# riboflow_paper

## Layout

| | |
|---|---|
| [`code/`](code/) | the programs, one directory per scientific function |
| [`config/`](config/) | sample manifest, panel manifest, and the published cohort's RiboFlow configuration |
| [`data/`](data/) | the shipped analysis tables the panels read, mirroring `results/` |
| [`results/`](results/) | regenerated output — **gitignored** |
| [`docs/`](docs/) | [HDF5 schema](docs/hdf5_schema.md) · [accessions](docs/accessions.tsv) · [numeric claims](docs/numeric_claims.tsv) |
| [`figures/`](figures/) | panel references and the author-assembled figures |
| [`supporting_information/`](supporting_information/S1_Table/) | the sample table and its generator |
| [`benchmark/`](benchmark/) | the runtime traces behind Table 1 |

Inside `code/`, `results/` and `data/`, the same names mean the same thing:

| | |
|---|---|
| [`ribo_seq_qc/`](code/ribo_seq_qc/) | read-length window and P-site offsets |
| [`coverage/`](code/coverage/) | BAMs → shared-coordinate HDF5 → concordance |
| [`ribo_rna/`](code/ribo_rna/) | Ribo-seq vs RNA-seq counts and route comparison |
| [`read_taxonomy/`](code/read_taxonomy/) | read taxonomy, concordance, reach, multimapper biotypes |
| [`alignment_fate/`](code/alignment_fate/) | per-read fates for a chosen transcript |
| [`panels/`](code/panels/) | one program per paper panel |

## Install

```bash
pip install -r requirements.txt
```



## Inputs

**This repository starts from BAMs.** the FASTQ-to-BAM step is a separate pipeline, described under
[Upstream processing](#upstream-processing) and not re-implemented here.


| input | flag | notes |
|---|---|---|
| Ribo-seq genome BAM (STAR) | `--genome-bam` | coordinate-sorted and indexed |
| Ribo-seq transcriptome BAM (Bowtie2) | `--transcriptome-bam` | reference names must be the APPRIS headers |
| RNA-seq genome + transcriptome BAMs | via `--bams` | the ribo-vs-RNA stage only |
| GENCODE GTF | `--gtf` | v34 for the published cohort; gzipped is fine |
| APPRIS transcript lengths | `--appris` | one line per transcript: the full transcriptome reference header (`ENST…|gene|…|UTR5:a-b|CDS:c-d|UTR3:e-f|`), a TAB, and the transcript length |
| actual-regions BED | `--regions` | *optional*. A coordinate cross-check only: region boundaries are taken from the reference header either way, and a disagreement fails the run rather than overriding them |
| read-length window + offsets | `--qc-genome` / `--qc-txome` | shipped for the published cohort under `data/ribo_seq_qc/` |

## Running it

```bash
# BAMs -> the shared-coordinate coverage HDF5
python code/coverage/build_shared_coverage.py --sample HeLa \
    --genome-bam … --transcriptome-bam … --gtf GTF --appris APPRIS \
    --qc-genome data/ribo_seq_qc/genome/tables/readlen_window_qc.csv \
    --qc-txome  data/ribo_seq_qc/transcriptome/tables/readlen_window_qc.csv \
    --annotation-cache results/.cache/annotation/coverage_annotation.pkl \
    --output results/coverage

# coverage -> concordance tables
python code/coverage/compute_coverage_concordance.py \
    --coverage results/coverage --output results/coverage/concordance

# coverage + a gene ID -> a figure, for any transcript
python code/panels/plot_transcript_coverage.py \
    --coverage-h5 results/coverage/HeLa.shared_coverage.h5 \
    --gene-id ENSG00000111640 --output results/plots/HeLa.GAPDH

# everything
python code/make_tables.py --bams DIR --gtf GTF --appris APPRIS --all
python code/make_panels.py --all
```

## What a run produces
```
results/
├── annotation/orf_catalog.tsv
├── coverage/<sample>.shared_coverage.h5, coverage_checksums.tsv, concordance/*.tsv[.gz]
├── ribo_seq_qc/{genome,transcriptome}/tables/*.csv
├── ribo_rna/*.tsv   read_taxonomy/**/*.tsv   alignment_fate/*.tsv
├── panels/*.pdf     the 14 panel assets
└── .cache/          rebuildable annotation caches and intermediates
```

## What is generated and not

A regenerated table and its shipped counterpart share the same relative path, so promoting
one is a prefix swap:

```
results/coverage/concordance/region_concordance_per_sample.tsv     regenerated, gitignored
data/   coverage/concordance/region_concordance_per_sample.tsv     shipped
```

- **[`results/`](results/)** is where outputs are stored.
- **[`data/`](data/)** is a snapshot of those same tables, panels can be reproduced with these alone.
- **Figures 3A and 3B are the one exception**: they read the complete HeLa coverage HDF5
  under `results/coverage/`, which must be regenerated.
- **Panels are always written to `results/panels/`**

| Workflow | Products | Panels |
|---|---|---|
| [Ribo-seq QC](code/ribo_seq_qc/) | read-length selection, offsets, CDS periodicity | 2A–B |
| [Shared coverage](code/coverage/) | genome/transcriptome coverage HDF5 | 3A–B |
| [Coverage concordance](code/coverage/) | per-transcript and pooled concordance | 3C–D |
| [Ribo/RNA comparison](code/ribo_rna/) | CDS-assigned counts per transcript and the route summary | 4A–C |
| [Read taxonomy](code/read_taxonomy/) | route counts, read unions, biotypes, reach | 5A–D |
| [Alignment fate](code/alignment_fate/) | transcript-specific alignment outcomes | 5E |

[`config/panel_manifest.yaml`](config/panel_manifest.yaml) is the only place figure and panel
letters appear. Figure 1 is an author-drawn schematic with no generator.

## Verify

```bash
python code/make_panels.py --all --verify    # 14 panels against stored references
python benchmark/summarize_benchmarks.py --check   # every value in Table 1
python code/coverage/coverage_schema.py --validate results/coverage/HeLa.shared_coverage.h5
```

## Upstream processing

Ribo-seq and RNA-seq reads were processed with **[RiboFlow_genome](https://github.com/DanielNguyener/riboflow_genome)**, a Nextflow DSL2 pipeline
extending [RiboFlow](https://github.com/ribosomeprofiling/riboflow) to align to the genome as
well as the transcriptome. This repository does not re-implement that step.

The 24 per-sample parameter files are at
[`config/published_cohort/riboflow_configs/`](config/published_cohort/riboflow_configs/) with
the [run script](config/published_cohort/run_riboflow_cohort.sh).

## References

| | |
|---|---|
| **RiboFlow-genome** | <https://github.com/DanielNguyener/riboflow_genome> — processed the Ribo-seq and RNA-seq data |
| **references_for_riboflow** | <https://github.com/ribosomeprofiling/references_for_riboflow> — `transcriptome/human/v2`, The transcriptome reference and region definitions, and the basis of the shared coordinate |
| **GENCODE release 34** | <https://www.gencodegenes.org/human/release_34.html> — mirrored in the reference set above |
| **ribobaser** |  <https://github.com/CenikLab/ribobaser> supplies the QC table used to select the sample panel |

Accessions: [`docs/accessions.tsv`](docs/accessions.tsv).