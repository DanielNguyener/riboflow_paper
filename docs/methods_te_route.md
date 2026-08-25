# Figure 4 methods

The mathematics behind `code/te_route/normalization.R`, `te_statistics.R` and
`plot_te_route_panels.py`. Notation: `c_{g,i}` is the raw count for transcript `g` in cell line
`i`; there are four count matrices, one per (assay, alignment route) with assay ∈ {Ribo, RNA}
and route ∈ {genome, transcriptome}.

---

## 1. CPM filtering

Counts are converted to counts per million **within each matrix**:

```
CPM_{g,i} = 1e6 · c_{g,i} / Σ_g c_{g,i}
```

The denominator is that matrix's own column sum. The four libraries differ in depth for
different reasons, and one shared denominator would make the gate stricter on the shallower
route as an artefact of depth rather than of expression.

A cell line **supports** a transcript when it clears the threshold in **all four matrices in
that same cell line**:

```
support_{g,i} = 1[ CPM^{Ribo,gen}_{g,i} > 1 ] · 1[ CPM^{RNA,gen}_{g,i} > 1 ]
              · 1[ CPM^{Ribo,txo}_{g,i} > 1 ] · 1[ CPM^{RNA,txo}_{g,i} > 1 ]
```

A transcript passes when `Σ_i support_{g,i} ≥ 12` — half the 24-line panel. This yields
**11,589 of 19,736** transcripts.

Requiring the four together *per line* is stricter than requiring each matrix to reach the
line count on its own, which would let a transcript qualify on Ribo in one set of lines and
on RNA in a disjoint set.

**The line requirement does the work, not the threshold.** At these depths 1 CPM is roughly
1–14 raw reads, so the threshold is a gate against the near-empty tail; moving it from 1 to 0
changes the passing set by ~1,354 transcripts, whereas raising the line requirement from 12
to 24 changes it by ~4,347.

CPM is used **only** to decide membership. No CPM value enters any plotted or tested quantity.

---

## 2. Median-of-ratios normalization, shared across routes

Size factors follow Anders & Huber (2010), computed in log space. For a count matrix `m`:

```
log2 ref_g   = mean_i log2 m_{g,i}          # per-transcript geometric mean
s_i          = 2 ^ median_g ( log2 m_{g,i} − log2 ref_g )
```

Working in logs is the same number as `median_g (m_{g,i} / ref_g)` — `2^x` is strictly
increasing and so commutes with the median — but never materialises the ratio matrix.

**The pseudo-reference needs non-zero rows.** A geometric mean is zero the moment any one
entry is zero, and a ratio against a zero reference is undefined. Estimation is therefore
restricted to transcripts non-zero in all four matrices and all 24 lines: **7,864** of the
11,589. Factors are then applied to all 11,589 — the filter chooses which rows *vote* on the
factor, not which rows appear in the output.

**One factor per (assay, library), shared by both routes.** For assay `a`, the estimator is
run on the geometric mean of the two routes,

```
G^a_{g,i} = sqrt( m^{a,genome}_{g,i} · m^{a,txome}_{g,i} )
```

giving `s^a_i`, which is applied to **both** routes of that assay.

### Why the choice of `G` cannot affect ΔTE

Under the shared scheme the same factor divides both routes, so inside any within-library
route comparison it cancels exactly:

```
log2( R^gen_{g,i} / s^a_i ) − log2( R^txo_{g,i} / s^a_i ) = log2 R^gen_{g,i} − log2 R^txo_{g,i}
```

Every quantity here — ΔRNA, ΔRibo and ΔTE alike — is a within-library difference of
routes, so all three are **raw log-ratios carrying no normalization at all**. The geometric
mean is therefore a presentational choice for the scaled tables, not a modelling one: it
cannot move any ΔTE result. It matters only for per-route TE taken on its own.

Under the alternative of four independently estimated factor sets this cancellation does not
occur, and ΔTE picks up a per-cell-line constant
`D_i = log2(s^{Ribo,gen}_i / s^{Ribo,txo}_i) − log2(s^{RNA,gen}_i / s^{RNA,txo}_i)`
that is a property of the estimation rather than of the alignment. That is the reason for the
shared scheme.

---

## 3. Log transformation and translation efficiency

```
TE^route_{g,i} = log2 R^route_{g,i} − log2 N^route_{g,i}     (R = Ribo, N = RNA)
```

equivalently `log2(RPF / mRNA)`. Differences of logs are used throughout rather than a ratio
of pooled totals: the difference is formed **within one library** and every library then
counts once, whereas pooling first would let the deepest libraries dominate.

**No pseudocount.** A zero is a missing measurement, not a small one; adding a constant would
invent a finite log-ratio out of it and would bias exactly the untranslated transcripts. A
(transcript, cell line) cell is used only where **all four** counts are strictly positive.
Consequently a transcript is averaged over **12–24** cell lines (median 24; 7,864 of 11,589
have all 24), not always 24.

---

## 4. ΔTE and its two assay halves

Per transcript and cell line, each as genome minus transcriptome:

```
ΔRNA_{g,i}  = log2 N^gen_{g,i} − log2 N^txo_{g,i}
ΔRibo_{g,i} = log2 R^gen_{g,i} − log2 R^txo_{g,i}
ΔTE_{g,i}   = TE^gen_{g,i} − TE^txo_{g,i} = ΔRibo_{g,i} − ΔRNA_{g,i}
```

The pairing is within cell line; the mean over lines is taken afterwards. Two consequences
used by the figure:

- `ΔTE = ΔRibo − ΔRNA` exactly, so in panel C the line **`y = x` is `ΔTE = 0`** and the
  vertical distance from it is a transcript's ΔTE.
- Negative ΔTE means **lower TE under genome alignment**.

---

## 5. Test, confidence interval, and multiple-testing correction

Each transcript's `n_g ∈ [12, 24]` per-line values are one sample. A **two-sided one-sample
t-test against zero** asks whether the route effect is consistently non-zero *across cell
lines*:

```
mean_g = (1/n_g) Σ_i ΔTE_{g,i}
SE_g   = sd_g / sqrt(n_g)
t_g    = mean_g / SE_g                       on n_g − 1 degrees of freedom
CI_g   = mean_g ± t_{0.975, n_g−1} · SE_g
```

The 95 % band in panel B is that same interval, so a band excluding zero and a nominal
p < 0.05 are the same statement.

Collapsing each cell line to one value before testing is what makes the replicates
independent: the two routes score the *same underlying reads*, and a model fitted jointly
over the four tables would treat that shared signal as independent observations. Collapsing
first puts the dependence inside the statistic rather than across the observations. The 24
cell lines are the replication; the test says nothing about significance within any one line.

**Benjamini–Hochberg** across all transcripts with a defined p. With p-values sorted
ascending and `m` tests:

```
padj_(k) = min_{j ≥ k} ( m · p_(j) / j ),   capped at 1
```

`padj` is an **adjusted p-value**. It is not the realised false-discovery rate, and it is not
a Storey q-value; those are different quantities. Here BH takes the count from 7,037 nominally
significant to 5,956.

The figure's highlighted set is `padj < 0.05` **and** `|mean ΔTE| > 1` — the second is an
effect-size gate, and 1 log₂ unit is a two-fold difference in estimated TE.

---

## 6. Route agreement, and what it cannot see

Panel A correlates `TE^gen_{·,i}` against `TE^txo_{·,i}` across transcripts, separately within
each cell line, on that line's usable transcripts (median 11,130; range 9,640–11,494).

A size factor enters log₂ TE as a term that is **constant across transcripts** within one
(cell line, route):

```
TE^route_{g,i} = log2( R^route_{g,i} / N^route_{g,i} ) − log2( s^{Ribo,route}_i / s^{RNA,route}_i )
```

Both Pearson and Spearman are invariant to adding a constant to either variable, so **these
correlations cannot change with the normalization scheme** — verified numerically at
8 × 10⁻¹⁶ for Pearson. Panel A therefore measures route agreement and cannot be used as
evidence for or against a normalization choice; that is also the sense in which a high
correlation is insensitive to a systematic offset (Bland & Altman, 1986), and why panels B
and C are needed alongside it.

---

## 7. Figure conventions

- **PLOS page.** The figure is one page inside PLOS Computational Biology's 7.5 × 8.75 in
  cap (border included): two rows, A | B above and C with its colorbar below. The axes box is
  solved as the largest square that fits both caps once the margins are taken out (B and C
  ≈3.4 in square; A to fig03D's proportions but never under 2.1 in, so its tick labels do
  not touch); a single-panel page reuses the identical box. Type is Arial at 10–11 pt with
  12 pt bold panel letters, inside PLOS's 8–12 pt window. Panel C's key sits inside the plane's
  empty upper-left corner, framed so its sample dot is not read as data. Besides PDF/PNG the program writes a flattened RGB TIFF at
  300 dpi with LZW compression and a 2-pt white border, and does not write one that exceeds
  2250 × 2625 px or 10 MB.
- **Panel B's y axis is symmetric-log**, linear within ±0.02 and logarithmic outside. The
  crossover sits *below* the interquartile range so the bulk falls in the expanded part; note
  that this exaggerates apparent steepness near zero.
- **Panel C colours points by local transcript density** — a smoothed 2-D histogram read back
  at each point. One point is one transcript; the colour counts its neighbours.
- **Labels are chosen by role, not rank**: the extreme on each axis, the only positive
  transcript, and the three largest RNA-driven housekeeping transcripts. The largest effects
  are mostly *not* in the housekeeping atlas, so ranking alone would show none of them.
  The manuscript's example genes (GAPDH, COMT, LRRFIP1) are always labelled, in red.

---

## References

- Anders S, Huber W. Differential expression analysis for sequence count data.
  *Genome Biology* 11:R106 (2010). — median-of-ratios
- Benjamini Y, Hochberg Y. Controlling the false discovery rate. *JRSS B* 57:289–300 (1995).
- Bland JM, Altman DG. Statistical methods for assessing agreement between two methods of
  clinical measurement. *Lancet* 327:307–310 (1986).
- Hounkpe BW *et al.* HRT Atlas v1.0. *Nucleic Acids Research* 49:D947–D955 (2021).
  — the housekeeping transcript lists in `data/te_route/housekeeping/`
