# Figure 3 methods

The math behind `code/te_route/normalization.R`, `te_statistics.R` and
`plot_te_route_panels.py`. `c_{g,i}` is the raw count for transcript g in cell line i.
Four count matrices, one per assay and route pair (Ribo or RNA, genome or transcriptome).

## 1. CPM filtering

```
CPM_{g,i} = 1e6 * c_{g,i} / sum_g c_{g,i}
```

Each matrix uses its own column sum: the four libraries differ in depth for different
reasons. A cell line supports a transcript when it clears 1 CPM in all four matrices in
that line:

```
support_{g,i} = 1[ CPM^{Ribo,gen}_{g,i} > 1 ] * 1[ CPM^{RNA,gen}_{g,i} > 1 ]
              * 1[ CPM^{Ribo,txo}_{g,i} > 1 ] * 1[ CPM^{RNA,txo}_{g,i} > 1 ]
```

A transcript passes with support in at least 12 of 24 cell lines: 11,589 of 19,736
transcripts. The line requirement does the work, not the threshold (1 CPM is ~1-14 raw
reads at these depths; threshold 1 -> 0 moves the passing set by ~1,354 transcripts,
line requirement 12 -> 24 by ~4,347). CPM decides membership only; no CPM value enters
any plotted or tested number.

## 2. Median-of-ratios normalization, shared across routes

Size factors follow Anders and Huber (2010), computed in log space:

```
log2 ref_g = mean_i log2 m_{g,i}          (per-transcript geometric mean)
s_i        = 2 ^ median_g ( log2 m_{g,i} - log2 ref_g )
```

Logs give the same number as the median of ratios (2^x is increasing) without building
the ratio matrix. The reference needs rows with no zeros, so factors are estimated on
the transcripts non-zero in all four matrices and all 24 lines (7,864 of 11,589) and
applied to all 11,589: the filter chooses which rows vote, not which rows appear.

One factor per assay and library, shared by both routes, estimated on the geometric mean
of the two routes:

```
G^a_{g,i} = sqrt( m^{a,genome}_{g,i} * m^{a,txome}_{g,i} )
```

The shared factor cancels exactly in any within-library route comparison:

```
log2( R^gen_{g,i} / s^a_i ) - log2( R^txo_{g,i} / s^a_i ) = log2 R^gen_{g,i} - log2 R^txo_{g,i}
```

so delta RNA, delta Ribo and delta TE carry no normalization at all. Four independently
estimated factor sets would inject a per-cell-line constant into delta TE that comes
from the estimation, not the alignment; that is the reason for the shared scheme.

## 3. Log transformation and translation efficiency

```
TE^route_{g,i} = log2 R^route_{g,i} - log2 N^route_{g,i}     (R = Ribo, N = RNA)
```

log2 of RPF over mRNA, formed within one library so every library counts once. No
pseudocount: a zero is a missing measurement, and a constant would invent a finite log
ratio for exactly the untranslated transcripts. A transcript/cell-line pair is used only
when all four counts are positive, so a transcript averages over 12-24 cell lines
(median 24; 7,864 of 11,589 have all 24).

## 4. Delta TE and its two assay halves

Per transcript and cell line, each as genome minus transcriptome:

```
dRNA_{g,i}  = log2 N^gen_{g,i} - log2 N^txo_{g,i}
dRibo_{g,i} = log2 R^gen_{g,i} - log2 R^txo_{g,i}
dTE_{g,i}   = TE^gen_{g,i} - TE^txo_{g,i} = dRibo_{g,i} - dRNA_{g,i}
```

Pairing is within cell line; the mean over lines comes afterwards. In panel C the line
y = x is exactly dTE = 0 and the vertical distance from it is a transcript's dTE.
Negative dTE means lower TE under genome alignment.

## 5. Test, confidence interval, and multiple testing

Each transcript's 12-24 per-line values are one sample; a two-sided one-sample t-test
against zero:

```
mean_g = (1/n_g) sum_i dTE_{g,i}
SE_g   = sd_g / sqrt(n_g)
t_g    = mean_g / SE_g                    on n_g - 1 degrees of freedom
CI_g   = mean_g +/- t_{0.975, n_g-1} * SE_g
```

Panel B's 95% band is that interval, so a band excluding zero and nominal p < 0.05 are
the same statement. Collapsing each cell line to one value first is what makes the
replicates independent: both routes score the same underlying reads. The 24 cell lines
are the replication; the test says nothing about any one line.

Benjamini-Hochberg runs across all transcripts with a defined p (p sorted ascending,
m tests):

```
padj_(k) = min over j >= k of ( m * p_(j) / j ),  capped at 1
```

padj is an adjusted p-value, not a realized FDR or q-value. BH takes 7,037 nominally
significant to 5,956. The figure highlights padj < 0.05 AND |mean dTE| > 1 (an
effect-size gate; 1 log2 unit = two-fold TE).

## 6. Route agreement, and what it cannot see

Panel A correlates genome-route TE against transcriptome-route TE across transcripts,
within each cell line, on that line's usable transcripts (median 11,130, range
9,640-11,494). A size factor enters log2 TE as a constant across transcripts within one
cell line and route; Pearson and Spearman are unchanged by adding a constant, so these
correlations cannot depend on the normalization scheme (verified numerically to 8e-16
for Pearson). A high correlation is also insensitive to a systematic offset (Bland and
Altman, 1986), which is why panels B and C exist.

## 7. Figure conventions

- One page inside PLOS Computational Biology's 7.5 x 8.75 inch cap: A and B on top,
  C with its colorbar below. The axes box is the largest square that fits (B and C
  ~3.4 inches; A follows fig02D's proportions, never under 2.1 inches). Arial 10-11 pt,
  12 pt bold panel letters. Panel C's key sits framed in the empty upper-left corner.
  Besides PDF and PNG the program writes a flattened RGB TIFF at 300 dpi, LZW, 2 pt
  white border, capped at 2250 x 2625 px and 10 MB.
- Panel B's y axis is symmetric-log: linear within 0.02 of zero, logarithmic outside,
  crossover below the interquartile range. This exaggerates apparent steepness near zero.
- Panel C colors points by local transcript density (smoothed 2-D histogram). One point
  is one transcript.
- Labels by role, not rank: the extreme on each axis, the only positive transcript, the
  three largest RNA-driven housekeeping transcripts. GAPDH, COMT and LRRFIP1 are always
  labeled, in red.

## References

- Anders S, Huber W. Genome Biology 11:R106 (2010). Median-of-ratios.
- Benjamini Y, Hochberg Y. JRSS B 57:289-300 (1995).
- Bland JM, Altman DG. Lancet 327:307-310 (1986).
- Hounkpe BW et al. HRT Atlas v1.0. Nucleic Acids Research 49:D947-D955 (2021).
  The housekeeping lists in data/te_route/housekeeping/.
