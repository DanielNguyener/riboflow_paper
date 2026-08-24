# riboflow_paper

Analysis and figure code for the RiboFlow_v2 manuscript: the same Ribo-seq libraries aligned
to the genome and to the transcriptome, compared route against route over a 24-cell-line
panel. Starts from RiboFlow_v2 BAMs; produces every analysis table, every panel and the five
published figures (Fig2–Fig6), Table 1 and S1 Table.

| | |
|---|---|
| [`code/`](code/README.md) | the programs, one directory per scientific function |
| [`config/`](config/README.md) | panel/figure manifest, cohort manifest, RiboFlow_v2 run configuration |
| [`data/`](data/README.md) | the shipped analysis tables the panels read (mirrors `results/`) |
| `results/` | regenerated output — gitignored |
| [`figures/`](figures/README.md) | panel references and the published figures |
| [`docs/`](docs/README.md) | methods, figure index, numeric claims, accessions |
| [`tests/`](tests/) | `python -m pytest tests -q` |
| [`supporting_information/`](supporting_information/S1_Table/README.md) | S1 Table and its generator |
| [`benchmark/`](benchmark/README.md) | the runtime traces behind Table 1 |
| [`get_coverage/`](get_coverage/README.md) | standalone full-read coverage tool (HPC) |

## Install

Python 3.9 (run as `python`), R ≥ 4 (base R only), Arial available to matplotlib.

```bash
pip install -r requirements.txt -r requirements-dev.txt
```

## Reproduce the figures from the shipped tables (no BAM)

```bash
python code/make_figures.py --all --check      # figures/published/Fig{2..6}.tif + _plos.pdf
python -m pytest tests -q                      # includes a clean-copy, byte-for-byte rebuild
```

Figure 3A/3B are the one exception: they read `results/coverage/HeLa.shared_coverage.h5`,
a ~25 MB product built from the HeLa BAMs (`code/coverage/build_shared_coverage.py`; see
[`docs/figures.md`](docs/figures.md)). Without it, `--figure 3` stops with that command.

## Reproduce the tables from BAMs

The RiboFlow_v2 BAMs for the 24-cell-line panel are deposited on Zenodo:
<https://doi.org/10.5281/zenodo.22083992>. Download them into a directory and pass it as `--bams DIR`.

```bash
python code/make_tables.py --bams DIR --validate                         # what would run
python code/make_tables.py --bams DIR --gtf GTF --appris APPRIS --all    # -> results/
python code/make_tables.py --bams DIR --gtf GTF --appris APPRIS --all --into-data   # -> data/
python code/make_figures.py --all --bams DIR --gtf GTF --appris APPRIS --into-data --check
```

| input | flag | notes |
|---|---|---|
| RiboFlow_v2 output tree | `--bams` | Zenodo [10.5281/zenodo.22083992](https://doi.org/10.5281/zenodo.22083992); `{s}/genome/alignment_ribo/merged/{s}.post_dedup.bam` etc.; paths may also come from `config/cohort_manifest.tsv` |
| GENCODE GTF | `--gtf` | v34 for the published cohort |
| APPRIS transcript lengths | `--appris` | the transcriptome reference headers (`references_for_riboflow`, `transcriptome/human/v2`), not redistributed |
| read-length window + offsets | `--qc-genome` / `--qc-txome` | shipped under `data/ribo_seq_qc/` |

Raw inputs have no default. Programs that open BAMs directly take flags, then
`RIBOFLOW_PAPER_{BAMS,GTF,APPRIS}`, then `config/local.yaml` (copy `config/inputs.example.yaml`).

## Figures

| figure | panels | tables from | generators |
|---|---|---|---|
| 2 | A–B read-length window, P-site offsets, CDS periodicity | `code/ribo_seq_qc/` | `code/panels/plot_readlen_psite_selection.py`, `plot_cds_periodicity_difference.py` |
| 3 | A–B COMT/GAPDH coverage, C–D concordance | `code/coverage/` | `code/panels/plot_transcript_coverage.py`, `plot_per_transcript_concordance.py`, `plot_pooled_concordance.py` |
| 4 | A–C translation efficiency by route | `code/ribo_rna/` + `code/te_route/` (R) | `code/te_route/plot_te_route_panels.py` |
| 5 | A–D read taxonomy over the cohort | `code/read_taxonomy/` | `code/panels/plot_fig05_plos_panels.py` |
| 6 | A gene read partition, B LRRFIP1 locus | `code/alignment_fate/` | `code/panels/plot_gene_partition.py`, `plot_locus_coverage.py` |

`config/panel_manifest.yaml` is the only place figure numbers appear. Figure 1 is an
author-drawn schematic with no generator. Details, geometry and checksums:
[`docs/figures.md`](docs/figures.md).

## Verify

```bash
python code/make_panels.py --all --verify          # 20 panels against figures/panel_references/
python code/assemble_figures.py --all --check      # PLOS spec: size, fonts, TIFF mode/dpi/LZW
python benchmark/summarize_benchmarks.py --check   # every value in Table 1
```

## Upstream processing

Reads were processed with **RiboFlow_v2**, a Nextflow DSL2 pipeline extending
[RiboFlow](https://github.com/ribosomeprofiling/riboflow) to align to the genome as well as
the transcriptome. The 24 per-sample parameter files and run script are in
[`config/published_cohort/`](config/published_cohort/riboflow_configs/README.md). The
transcriptome reference is `references_for_riboflow` `transcriptome/human/v2`
(<https://github.com/ribosomeprofiling/references_for_riboflow>), mirroring GENCODE
release 34. Sample selection used the QC table from
[ribobaser](https://github.com/CenikLab/ribobaser). Accessions: [`docs/accessions.tsv`](docs/accessions.tsv).

## Citation

See [`CITATION.cff`](CITATION.cff). License: [`LICENSE`](LICENSE).
