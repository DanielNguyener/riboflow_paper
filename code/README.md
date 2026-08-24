# code

One directory per scientific function; `data/` and `results/` use the same names.

| | |
|---|---|
| `make_tables.py` | every analysis table from BAMs, as ordered stages (`--validate`, `--stages`, `--into-data`) |
| `make_panels.py` | every panel from `config/panel_manifest.yaml` (`--verify`, `--accept`; `fit:` panels are re-rendered to width) |
| `assemble_figures.py` | panels → `figures/published/FigN.tif` + `_plos.pdf`, per the manifest's `figures:` block (`--check`) |
| `make_figures.py` | tables (optional) → panels → figures, end to end |
| `common/` | BAM discovery and templates, raw-input resolution (`inputs.py`), annotation cache, region classification, ORF catalog |
| `ribo_seq_qc/` | read-length window and P-site offsets, both routes (Figure 2) |
| `coverage/` | BAMs → shared-coordinate HDF5 (schema 3, `docs/hdf5_schema.md`) → concordance tables (Figure 3) |
| `ribo_rna/` | CDS-assigned counts per transcript; the cohort count matrices (Figure 4 input) |
| `te_route/` | CPM gate, shared median-of-ratios, ΔTE statistics (R); the Figure 4 plotter |
| `read_taxonomy/` | read taxonomy, alignment concordance, reach, multimapper tie biotypes (Figure 5) |
| `alignment_fate/` | gene read partition and the LRRFIP1 locus artifact (Figure 6) |
| `panels/` | one program per panel, `panel_style.py`, `figure_io.py` (ink crop, fit loop, composition, TIFF) |

Stages (`make_tables.py`): `annotation → qc → offsets, orf_catalog, coverage → concordance,
te_counts → te_normalize → te_stats, taxonomy, alignment_concordance → reach,
multimap_biotype, gene_partition, locus`. Run Python 3.9 as `python`.
