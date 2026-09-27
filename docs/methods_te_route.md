# Figure 3 methods

The math behind `code/te_route/normalization.R`, `te_statistics.R` and
`plot_te_route_panels.py`. Notation: `c_{g,i}` is the raw count for transcript g in cell
line i. There are four count matrices, one per assay and route pair (Ribo or RNA, genome
or transcriptome).

## 1. CPM filtering

Counts are converted to counts per million within each matrix:

```
CPM_{g,i} = 1e6 * c_{g,i} / sum_g c_{g,i}
```

Each matrix uses its own column sum. The four libraries differ in depth for different
reasons, and one shared denominator would punish the shallower route for its depth
instead of its expression.

A cell line supports a transcript when the transcript clears 1 CPM in all four matrices
in that same cell line:

```
support_{g,i} = 1[ CPM^{Ribo,gen}_{g,i} > 1 ] * 1[ CPM^{RNA,gen}_{g,i} > 1 ]
              * 1[ CPM^{Ribo,txo}_{g,i} > 1 ] * 1[ CPM^{RNA,txo}_{g,i} > 1 ]
```

A transcript passes when at least 12 of the 24 cell lines support it. This keeps 11,589
of 19,736 transcripts. Requiring all four matrices per line is stricter than letting a
transcript qualify on Ribo in one set of lines and on RNA in another.

The line requirement does the work, not the threshold. At these depths 1 CPM is roughly
1 to 14 raw reads. Moving the threshold from 1 to 0 changes the passing set by about
1,354 transcripts; raising the line requirement from 12 to 24 changes it by about 4,347.

CPM is used only to decide membership. No CPM value enters any plotted or tested number.

## 2. Median-of-ratios normalization, shared across routes

Size factors follow Anders and Huber (2010), computed in log space. For a count matrix m:

```
log2 ref_g = mean_i log2 m_{g,i}          (per-transcript geometric mean)
s_i        = 2 ^ median_g ( log2 m_{g,i} - log2 ref_g )
```

Working in logs gives the same number as the median of the ratios, because 2^x is
increasing and so passes through the median, but it never builds the ratio matrix.

The reference needs rows with no zeros. A geometric mean is zero as soon as one entry is
zero, and a ratio against a zero reference is undefined. So the factors are estimated on
the transcripts that are non-zero in all four matrices and all 24 lines: 7,864 of the
11,589. The factors are then applied to all 11,589. The filter chooses which rows vote on
the factor, not which rows appear in the output.

There is one factor per assay and library, shared by both routes. For assay a the
estimator runs on the geometric mean of the two routes:

```
G^a_{g,i} = sqrt( m^{a,genome}_{g,i} * m^{a,txome}_{g,i} )
```

The resulting s^a_i is applied to both routes of that assay.

Why this choice cannot affect delta TE: the same factor divides both routes, so it
cancels exactly in any within-library route comparison.

```
log2( R^gen_{g,i} / s^a_i ) - log2( R^txo_{g,i} / s^a_i ) = log2 R^gen_{g,i} - log2 R^txo_{g,i}
```

Delta RNA, delta Ribo and delta TE are all within-library differences of routes, so all
three carry no normalization at all. The shared factor matters only for per-route TE
taken on its own. With four independently estimated factor sets the cancellation would
not happen, and delta TE would pick up a per-cell-line constant that comes from the
estimation rather than from the alignment. That is the reason for the shared scheme.

## 3. Log transformation and translation efficiency

```
TE^route_{g,i} = log2 R^route_{g,i} - log2 N^route_{g,i}     (R = Ribo, N = RNA)
```

This is log2 of RPF over mRNA. Differences of logs are used throughout, formed within one
library, so every library counts once; pooling totals first would let the deepest
libraries dominate.

No pseudocount is added. A zero is a missing measurement, not a small one, and a constant
would invent a finite log ratio for exactly the untranslated transcripts. A transcript
and cell line pair is used only when all four counts are positive. As a result a
transcript is averaged over 12 to 24 cell lines (median 24; 7,864 of 11,589 have all 24).

## 4. Delta TE and its two assay halves

Per transcript and cell line, each as genome minus transcriptome:

```
dRNA_{g,i}  = log2 N^gen_{g,i} - log2 N^txo_{g,i}
dRibo_{g,i} = log2 R^gen_{g,i} - log2 R^txo_{g,i}
dTE_{g,i}   = TE^gen_{g,i} - TE^txo_{g,i} = dRibo_{g,i} - dRNA_{g,i}
```

The pairing is within cell line; the mean over lines is taken afterwards. Two facts the
figure uses: in panel C the line y = x is exactly dTE = 0, and the vertical distance from
it is a transcript's dTE. Negative dTE means lower TE under genome alignment.

## 5. Test, confidence interval, and multiple testing

Each transcript's 12 to 24 per-line values are one sample. A two-sided one-sample t-test
against zero asks whether the route effect is consistently non-zero across cell lines:

```
mean_g = (1/n_g) sum_i dTE_{g,i}
SE_g   = sd_g / sqrt(n_g)
t_g    = mean_g / SE_g                    on n_g - 1 degrees of freedom
CI_g   = mean_g +/- t_{0.975, n_g-1} * SE_g
```

The 95 percent band in panel B is that same interval, so a band excluding zero and a
nominal p below 0.05 are the same statement.

Collapsing each cell line to one value before testing is what makes the replicates
independent. The two routes score the same underlying reads, and a model fitted jointly
over the four tables would treat that shared signal as independent observations. The 24
cell lines are the replication; the test says nothing about any one line.

Benjamini-Hochberg runs across all transcripts with a defined p. With p-values sorted
ascending and m tests:

```
padj_(k) = min over j >= k of ( m * p_(j) / j ),  capped at 1
```

padj is an adjusted p-value, not a realized false-discovery rate and not a q-value. Here
BH takes the count from 7,037 nominally significant to 5,956. The figure highlights
transcripts with padj below 0.05 and absolute mean dTE above 1; the second condition is
an effect-size gate, and 1 log2 unit is a two-fold difference in estimated TE.

## 6. Route agreement, and what it cannot see

Panel A correlates TE on the genome route against TE on the transcriptome route across
transcripts, within each cell line, on that line's usable transcripts (median 11,130,
range 9,640 to 11,494).

A size factor enters log2 TE as a term that is constant across transcripts within one
cell line and route. Pearson and Spearman are both unchanged by adding a constant to
either variable, so these correlations cannot change with the normalization scheme
(verified numerically to 8e-16 for Pearson). Panel A therefore measures route agreement
only. A high correlation is also insensitive to a systematic offset (Bland and Altman,
1986), which is why panels B and C are needed alongside it.

## 7. Figure conventions

- The figure is one page inside PLOS Computational Biology's 7.5 x 8.75 inch cap: A and B
  on the top row, C with its colorbar below. The axes box is the largest square that fits
  both caps after margins (B and C about 3.4 inches square; A follows fig02D's
  proportions but never drops under 2.1 inches). Type is Arial at 10 to 11 pt with 12 pt
  bold panel letters. Panel C's key sits inside the empty upper-left corner of the plane,
  framed so its sample dot is not read as data. Besides PDF and PNG the program writes a
  flattened RGB TIFF at 300 dpi with LZW compression and a 2 pt white border, capped at
  2250 x 2625 px and 10 MB.
- Panel B's y axis is symmetric-log: linear within 0.02 of zero and logarithmic outside.
  The crossover sits below the interquartile range so the bulk falls in the expanded
  part. This exaggerates apparent steepness near zero.
- Panel C colors points by local transcript density, a smoothed 2-D histogram read back
  at each point. One point is one transcript.
- Labels are chosen by role, not rank: the extreme on each axis, the only positive
  transcript, and the three largest RNA-driven housekeeping transcripts. The manuscript's
  example genes (GAPDH, COMT, LRRFIP1) are always labeled, in red.

## References

- Anders S, Huber W. Differential expression analysis for sequence count data.
  Genome Biology 11:R106 (2010). Median-of-ratios.
- Benjamini Y, Hochberg Y. Controlling the false discovery rate. JRSS B 57:289-300 (1995).
- Bland JM, Altman DG. Statistical methods for assessing agreement between two methods of
  clinical measurement. Lancet 327:307-310 (1986).
- Hounkpe BW et al. HRT Atlas v1.0. Nucleic Acids Research 49:D947-D955 (2021).
  The housekeeping transcript lists in data/te_route/housekeeping/.
