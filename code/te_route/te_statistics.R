#!/usr/bin/env Rscript
# Per-transcript delta TE (and its two assay halves) and per-cell-line route agreement.
# Reads the scaled matrices and the ORF catalog, writes the two statistics tables. Base R only.
# Mathematics: docs/methods_te_route.md.
#
#   Rscript code/te_route/te_statistics.R [--normalized DIR] [--orf-catalog FILE] [--output DIR]

# Arguments: --key value pairs; every default is repository-relative.
here <- dirname(normalizePath(sub("^--file=", "",
                                  grep("^--file=", commandArgs(FALSE), value = TRUE))[1]))
root <- dirname(dirname(here))
parse_args <- function(defaults) {
  raw <- commandArgs(TRUE)
  if (length(raw) && raw[1] %in% c("--help", "-h")) {
    cat("usage: Rscript", basename(sub("^--file=", "",
        grep("^--file=", commandArgs(FALSE), value = TRUE))[1]),
        paste(sprintf("[--%s PATH]", names(defaults)), collapse = " "), "\n")
    for (k in names(defaults)) cat(sprintf("  --%-14s default %s\n", k, defaults[[k]]))
    quit(status = 0)
  }
  if (length(raw) %% 2 != 0) stop("arguments come in --key value pairs")
  for (i in seq_len(length(raw) %/% 2) * 2 - 1) {
    key <- sub("^--", "", raw[i])
    if (!key %in% names(defaults)) stop("unknown argument --", key)
    defaults[[key]] <- raw[i + 1]
  }
  defaults
}
opts <- parse_args(list(
  normalized = file.path(root, "results", "te_route", "normalized"),
  "orf-catalog" = file.path(root, "data", "annotation", "orf_catalog.tsv"),
  output = file.path(root, "results", "te_route", "tables")))

CONF <- 0.95

read_scaled <- function(name) {
  frame <- utils::read.csv(file.path(opts$normalized, name),
                           check.names = FALSE, stringsAsFactors = FALSE)
  m <- as.matrix(frame[, -1, drop = FALSE])
  storage.mode(m) <- "double"
  rownames(m) <- as.character(frame[[1]])
  m
}

tab <- list(genome_ribo = read_scaled("ribo_scaled_genome.csv"),
            genome_rna  = read_scaled("rna_scaled_genome.csv"),
            txome_ribo  = read_scaled("ribo_scaled_txome.csv"),
            txome_rna   = read_scaled("rna_scaled_txome.csv"))
ids <- rownames(tab[[1]])
samples <- colnames(tab[[1]])

# All four must be strictly positive. A zero is a missing measurement; a pseudocount would
# invent a finite delta out of one, so the cell is dropped and the transcript's n falls.
usable <- (tab$genome_ribo > 0) & (tab$genome_rna > 0) &
          (tab$txome_ribo > 0) & (tab$txome_rna > 0)

d_rna  <- log2(tab$genome_rna)  - log2(tab$txome_rna)
d_ribo <- log2(tab$genome_ribo) - log2(tab$txome_ribo)
d_te   <- d_ribo - d_rna
for (nm in c("d_rna", "d_ribo", "d_te")) {
  m <- get(nm); m[!usable] <- NA_real_; assign(nm, m)
}

n_lines <- rowSums(usable)
cat(sprintf("%s transcripts x %d cell lines; n per transcript %d-%d (median %g)\n",
            format(length(ids), big.mark = ","), length(samples),
            min(n_lines), max(n_lines), stats::median(n_lines)))

# ── per transcript: mean, SD, 95% t interval ───────────────────────────────────────────────

summarise <- function(d) {
  present <- !is.na(d)
  n <- rowSums(present)
  d0 <- d; d0[!present] <- 0
  mean_d <- rowSums(d0) / n
  sd_d <- sqrt(pmax(rowSums(d0^2) - n * mean_d^2, 0) / (n - 1))
  se <- sd_d / sqrt(n)
  tcrit <- stats::qt(1 - (1 - CONF) / 2, df = n - 1)
  list(mean = mean_d, sd = sd_d, ci_low = mean_d - tcrit * se, ci_high = mean_d + tcrit * se,
       t = mean_d / se, p = 2 * stats::pt(-abs(mean_d / se), df = n - 1))
}

s_rna  <- summarise(d_rna)
s_ribo <- summarise(d_ribo)
s_te   <- summarise(d_te)

# Benjamini-Hochberg across every transcript with a defined p.
padj <- stats::p.adjust(s_te$p, method = "BH")

catalog <- utils::read.delim(opts[["orf-catalog"]],
                             stringsAsFactors = FALSE)
catalog <- catalog[!duplicated(catalog$transcript_id), c("transcript_id", "gene_name")]

out <- opts$output
dir.create(out, recursive = TRUE, showWarnings = FALSE)

per_gene <- data.frame(
  transcript_id = ids,
  gene_name = catalog$gene_name[match(ids, catalog$transcript_id)],
  n_lines = as.integer(n_lines),
  drna_mean = s_rna$mean,   drna_sd = s_rna$sd,
  dribo_mean = s_ribo$mean, dribo_sd = s_ribo$sd,
  dte_mean = s_te$mean,     dte_sd = s_te$sd,
  dte_ci_low = s_te$ci_low, dte_ci_high = s_te$ci_high,
  dte_t = s_te$t, dte_p = s_te$p, dte_padj = padj,
  stringsAsFactors = FALSE)
op <- options(digits = 17L)
utils::write.table(per_gene, file.path(out, "per_gene_delta.tsv"), sep = "\t",
                   quote = FALSE, row.names = FALSE, na = "NA")
options(op)

# ── per cell line: do the two routes agree on TE itself ────────────────────────────────────
# Correlated on log2 TE over that line's usable transcripts.

rows <- lapply(samples, function(s) {
  k <- usable[, s]
  lg <- log2(tab$genome_ribo[k, s]) - log2(tab$genome_rna[k, s])
  lt <- log2(tab$txome_ribo[k, s])  - log2(tab$txome_rna[k, s])
  data.frame(sample = s, n_transcripts = sum(k),
             spearman_rho = stats::cor(lg, lt, method = "spearman"),
             pearson_r = stats::cor(lg, lt), stringsAsFactors = FALSE)
})
correlation <- do.call(rbind, rows)
op <- options(digits = 17L)
utils::write.table(correlation, file.path(out, "route_correlation.tsv"), sep = "\t",
                   quote = FALSE, row.names = FALSE)
options(op)

# ── summary ────────────────────────────────────────────────────────────────────────────────

q <- stats::quantile(s_te$mean, c(0.25, 0.75))
big <- !is.na(padj) & padj < 0.05 & abs(s_te$mean) > 1
cat(sprintf("delta TE: median %+.4f  IQR %+.4f to %+.4f\n",
            stats::median(s_te$mean), q[1], q[2]))
cat(sprintf("  lower under genome %s   higher %s\n",
            format(sum(s_te$mean < 0), big.mark = ","),
            format(sum(s_te$mean > 0), big.mark = ",")))
cat(sprintf("  padj < 0.05 and |delta TE| > 1: %d  (%d negative, %d positive)\n",
            sum(big), sum(big & s_te$mean < 0), sum(big & s_te$mean > 0)))
cat(sprintf("route agreement: median Spearman %.4f  median Pearson %.4f\n",
            stats::median(correlation$spearman_rho), stats::median(correlation$pearson_r)))
cat(sprintf("wrote %s\n", out))
