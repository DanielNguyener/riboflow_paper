# riboflow_paper

[![DOI](https://zenodo.org/badge/1321237103.svg)](https://doi.org/10.5281/zenodo.22102431)

Code and tables for the RiboFlow_v2 manuscript: ribosome profiling and matched RNA-seq
from 24 human cell lines, aligned to the genome and to the transcriptome, with the two
alignment routes compared. This repository turns the RiboFlow_v2 alignments into the
analysis tables, Figures 2–6, S1 Fig, and S1 Table.

Read processing is done by the separate
[RiboFlow_v2](https://github.com/ribosomeprofiling/riboflow) pipeline, with the
configurations in
[`config/published_cohort/`](config/published_cohort/riboflow_configs/README.md).

## What is where

| | |
|---|---|
| `code/make_tables.py` | BAMs → analysis tables (`results/`; `--into-data` copies them to `data/`) |
| `code/make_panels.py` | tables → panel PDFs (`results/panels/`) |
| `code/assemble_figures.py` | panels → `figures/published/FigN.tif` |
| `code/make_figures.py` | runs all three |
| `code/ribo_seq_qc/` | read-length selection and P-site QC (S1 Fig) |
| `code/coverage/` | per-transcript coverage on both routes (Figure 2) |
| `code/ribo_rna/` + `code/te_route/` | CDS counts and translation-efficiency statistics, R (Figure 3) |
| `code/read_categories/` | the five read categories (Figures 4–6); see its [README](code/read_categories/README.md) |
| `code/panels/` | one plotting script per panel |
| `code/common/` | shared input handling and annotation tables |
| `config/panel_manifest.yaml` | which panel is drawn from which table, and how figures are composed |
| `data/` | the shipped tables, one directory per analysis |
| `results/` | regenerated output (not tracked) |
| `figures/` | per-panel reference PDFs and the published figures |
| `docs/` | Figure 3 math, the coverage file format, and every published number with its source |
| `supporting_information/S1_Table/` | the sample table and its generator |

## The five read categories (Figures 4–6)

Every read ID that aligns on either route falls into one of five categories, always in
this order:

| key | meaning |
|---|---|
| SH-U | on both routes, maps uniquely to the genome |
| SH-M | on both routes, multimaps in the genome |
| GO-U | genome only, unique |
| GO-M | genome only, multimapping |
| TO | transcriptome only |

`code/read_categories/categories.py` defines the keys, names and colours once, and every
figure uses it. For these figures a read is "on the transcriptome route" when it has a
primary alignment in the post-dedup BAM (RiboFlow_v2's MAPQ ≥ 10 filter). S1 Fig and
Figures 2–3 use a stricter rule instead: MAPQ ≥ 42.

## Setup

Python 3.9, R ≥ 4 (base only), and the Arial font.

```bash
pip install -r requirements.txt
```

## Reproducing everything

Two steps. Step 1 builds the analysis tables from the BAM files. Step 2 builds the
figures from the tables.

**Step 1: BAMs to tables.** Download the RiboFlow_v2 alignments (Zenodo link below) and
run:

```bash
python code/make_tables.py --bams DIR --gtf GTF --appris APPRIS --all --into-data
```

This writes every analysis table to `data/` and the per-sample coverage files to
`results/coverage/` (one HDF5 per sample; format in
[`docs/hdf5_schema.md`](docs/hdf5_schema.md)). It takes hours. If you already have the
coverage HDF5 files, or you only want the figures, skip this step: the repository ships
every table the figures read, and only Figure 2A/2B needs an HDF5
(`results/coverage/HeLa.shared_coverage.h5`).

**Step 2: tables to figures.**

```bash
python code/make_figures.py --all --check
```

This draws every panel from `data/` and writes
`figures/published/{Fig2,Fig3,Fig4,Fig5,Fig6,S1_Fig}.tif`. To compare the panels against
the shipped references, run `python code/make_panels.py --all --verify`.

## External inputs

| input | identifier | use |
|---|---|---|
| RiboFlow_v2 alignments, 24 cell lines | Zenodo [10.5281/zenodo.22083992](https://doi.org/10.5281/zenodo.22083992) | `--bams DIR`; layout in `config/cohort_manifest.tsv` |
| Sequencing data | GEO accessions in [`docs/accessions.tsv`](docs/accessions.tsv) | input to RiboFlow_v2 |
| GENCODE annotation | release 34 (GRCh38) | `--gtf` |
| APPRIS principal-isoform transcriptome | [`references_for_riboflow`](https://github.com/ribosomeprofiling/references_for_riboflow), `transcriptome/human/v2` | RiboFlow_v2 reference; its transcript-lengths table is `--appris` |
| Housekeeping gene lists | HRT Atlas v1.0 (`data/te_route/housekeeping/`) | Figure 3C labels |
| Consensus pseudogenes | GENCODE release 34, Yale-UCSC 2-way consensus set (`data/clustering/gencode.v34.2wayconspseudos.gtf.gz`) | Figure 6E |
| Sample QC table | [ribobaser](https://github.com/CenikLab/ribobaser) | S1 Table generator (not redistributed) |

## Code and data availability

Code: this repository (MIT), archived at Zenodo
[10.5281/zenodo.22102431](https://doi.org/10.5281/zenodo.22102431). Tables: `data/`.
Alignments: Zenodo [10.5281/zenodo.22083992](https://doi.org/10.5281/zenodo.22083992).
Raw reads: GEO ([`docs/accessions.tsv`](docs/accessions.tsv)). Every number in the
manuscript is listed in [`docs/numeric_claims.tsv`](docs/numeric_claims.tsv) with the
table and rule that produces it.

## Citation and license

Citation information will be added at release. License: [`LICENSE`](LICENSE).
