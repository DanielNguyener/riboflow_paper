# S1 Table — the balanced analysis panel

> **S1 Table. Human Ribo-seq libraries and matched RNA-seq datasets included in the balanced
> analysis panel.**

`samples.csv` — 24 rows, 23 columns

```bash
python supporting_information/S1_Table/build_s1_table.py \
    --rda "$RIBOBASER_RDA" --xlsx /path/to/41587_2025_2718_MOESM3_ESM.xlsx \
    --output /tmp/samples.csv
```