# te_route — does the alignment route change a transcript's translation efficiency?

```bash
Rscript code/te_route/normalization.R    # data/ribo_rna/counts/ -> results/te_route/normalized/
Rscript code/te_route/te_statistics.R    # + data/annotation/orf_catalog.tsv -> results/te_route/tables/
python  code/make_panels.py fig04        # results/panels/fig04_te_route_combined.pdf
```

| | |
|---|---|
| `normalization.R` | CPM gate (> 1 CPM in all four matrices in ≥ 12 of 24 lines), then one median-of-ratios size factor per (assay, library) shared by both routes; asserts 11,589 gated / 7,864 estimation rows |
| `te_statistics.R` | ΔRNA, ΔRibo, ΔTE per transcript × line; per-transcript mean, SD, 95 % t CI, t-test, BH; per-line route correlation |
| `te_panel_style.py` | PLOS type scale (Arial 10–12 pt), page caps, the TIFF writer — deliberately separate from `panels/panel_style.py`, whose sizes differ |
| `plot_te_route_panels.py` | the combined page and each single panel from one solved geometry (`--panel combined|A|B|C`) |

Inputs: `data/ribo_rna/counts/{ribo,rna}_counts_{genome,txome}.csv` (19,736 × 24, from
`code/ribo_rna/build_count_matrices.py`), `data/annotation/orf_catalog.tsv`,
`data/te_route/housekeeping/*.csv` (HRT Atlas v1.0). Base R, no packages. Expected numbers
and mathematics: [`docs/methods_te_route.md`](../../docs/methods_te_route.md),
`docs/numeric_claims.tsv` C13–C18.
