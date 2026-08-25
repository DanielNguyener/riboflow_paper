# S1 Table

S1 Table. Human Ribo-seq libraries and matched RNA-seq datasets included in the analysis
panel. `samples.csv`: 24 rows, 23 columns.

Regenerate from the ribobaser QC table (`.rda`) and the supplementary table of
Liu et al. 2025 (`41587_2025_2718_MOESM3_ESM.xlsx`), neither of which is redistributed:

```bash
python supporting_information/S1_Table/build_s1_table.py \
    --rda "$RIBOBASER_RDA" --xlsx /path/to/41587_2025_2718_MOESM3_ESM.xlsx \
    --output samples.csv
```
