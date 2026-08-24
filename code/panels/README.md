# panels

One program per panel; each takes its inputs as flags, writes `--output.<format>` and refuses
to overwrite without `--force`. `config/panel_manifest.yaml` is the only caller that matters.

| | |
|---|---|
| `panel_style.py` | the shared type scale (Arial 11/11/9/9 pt), colours, `save()`, `legend_below()` |
| `figure_io.py` | ink box, fit loop, 1:1 composition, PyMuPDF TIFF export (Figures 5, 6) |
| `_fig02_common.py`, `fig05_common.py`, `_fig05_side_panel.py` | helpers shared within a figure |
| `plot_readlen_psite_selection.py`, `plot_cds_periodicity_difference.py` | Figure 2 |
| `plot_transcript_coverage.py`, `plot_per_transcript_concordance.py`, `plot_pooled_concordance.py` | Figure 3 (`plot_transcript_coverage.py --gene-id` works for any transcript in a coverage HDF5) |
| `plot_route_read_counts.py`, `plot_read_id_union.py`, `plot_multimap_biotype.py`, `plot_nonselected_isoform_reach.py` | Figure 5, natural size |
| `plot_fig05_plos_panels.py` | the same four at page size: 8 pt type, shrunken shared plot box, rebound in-process |
| `plot_gene_read_partition.py`, `plot_gene_partition.py` | Figure 6A (the fold and `draw()`; the compact-table front end) |
| `plot_locus_coverage.py` | Figure 6B |

Figure 4's plotter lives in `code/te_route/` with its own style module.
