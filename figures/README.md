# figures

| | |
|---|---|
| `published/Fig{2..6}.tif`, `Fig{2..6}_plos.pdf` | the submitted figures; rebuilt by `python code/assemble_figures.py --all --check` and compared byte for byte by `tests/test_clean_copy.py` |
| `panel_references/*.pdf` | one reference per manifest panel; `python code/make_panels.py --all --verify` compares, `--accept` promotes |

Nothing here is written by the pipeline except through those two explicit commands.
Figure 1 is an author-drawn schematic and has no generator.
