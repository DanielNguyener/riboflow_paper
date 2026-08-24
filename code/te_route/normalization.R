#!/usr/bin/env Rscript
# CPM gate + assay-shared median-of-ratios normalization.
# Reads the four count matrices, writes the scaled matrices + size factors. Base R only.
# Mathematics: docs/methods_te_route.md.
#
#   Rscript code/te_route/normalization.R [--counts DIR] [--output DIR]

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
opts <- parse_args(list(counts = file.path(root, "data", "ribo_rna", "counts"),
                        output = file.path(root, "results", "te_route", "normalized")))

MIN_CPM      <- 1
MIN_LINES    <- 12L
N_GATED      <- 11589L   # asserted, not assumed: this folder must not drift from the analysis
N_ESTIMATION <- 7864L

MATRICES <- c(genome_ribo = "ribo_counts_genome.csv", genome_rna = "rna_counts_genome.csv",
              txome_ribo  = "ribo_counts_txome.csv",  txome_rna  = "rna_counts_txome.csv")

read_counts <- function(path) {
  frame <- utils::read.csv(path, check.names = FALSE, stringsAsFactors = FALSE)
  ids <- as.character(frame[[1]])
  if (anyDuplicated(ids)) stop("duplicate transcript ids in ", basename(path))
  m <- as.matrix(frame[, -1, drop = FALSE])
  storage.mode(m) <- "double"
  rownames(m) <- ids
  m
}

# ── inputs ─────────────────────────────────────────────────────────────────────────────────

raw <- lapply(MATRICES, function(f) read_counts(file.path(opts$counts, f)))
names(raw) <- names(MATRICES)

index <- rownames(raw[[1]])
for (nm in names(raw)) {
  if (!identical(rownames(raw[[nm]]), index)) stop(nm, " has a different transcript index")
}
# Columns are intersected BY NAME; pairing positionally would divide one library's ribo by
# another's rna the first time a matrix is rewritten.
samples <- Reduce(intersect, lapply(raw, colnames))
raw <- lapply(raw, function(m) m[, samples, drop = FALSE])

cat(sprintf("input: %s transcripts x %d cell lines\n",
            format(length(index), big.mark = ","), length(samples)))

# ── CPM gate ───────────────────────────────────────────────────────────────────────────────
# CPM against each matrix's OWN column sum; a line supports a transcript only when all four
# matrices clear the threshold in that same line.

supported <- matrix(TRUE, nrow = length(index), ncol = length(samples))
for (nm in names(raw)) {
  cpm <- sweep(raw[[nm]], 2, colSums(raw[[nm]]), "/") * 1e6
  supported <- supported & (cpm > MIN_CPM)
}
keep <- rowSums(supported) >= MIN_LINES
gated <- index[keep]
cat(sprintf("CPM gate (> %g CPM in all four matrices, >= %d lines): %s of %s\n",
            MIN_CPM, MIN_LINES, format(length(gated), big.mark = ","),
            format(length(index), big.mark = ",")))
if (length(gated) != N_GATED) {
  stop("expected ", N_GATED, " gated transcripts, got ", length(gated))
}

raw <- lapply(raw, function(m) m[gated, , drop = FALSE])

# ── size factors ───────────────────────────────────────────────────────────────────────────
# One factor per (assay, library), shared by both routes, estimated from the geometric mean
# of the two routes. Restricted to rows non-zero everywhere, because a geometric-mean
# reference is undefined the moment any entry is zero.

zero_free <- Reduce(`&`, lapply(raw, function(m) apply(m > 0, 1L, all)))
cat(sprintf("estimation set (non-zero in all four matrices and all lines): %s\n",
            format(sum(zero_free), big.mark = ",")))
if (sum(zero_free) != N_ESTIMATION) {
  stop("expected ", N_ESTIMATION, " estimation rows, got ", sum(zero_free))
}

median_of_ratios <- function(counts) {
  log_counts <- log2(counts)
  log_ref <- rowMeans(log_counts)                       # per-feature geometric mean
  2 ^ apply(log_counts - log_ref, 2L, stats::median)    # median ratio to that reference
}

size_factors <- list()
scaled <- list()
for (assay in c("ribo", "rna")) {
  g <- raw[[paste0("genome_", assay)]][zero_free, , drop = FALSE]
  t <- raw[[paste0("txome_", assay)]][zero_free, , drop = FALSE]
  sf <- median_of_ratios(2 ^ ((log2(g) + log2(t)) / 2))
  if (!all(is.finite(sf)) || any(sf <= 0)) stop("non-finite size factor for ", assay)
  size_factors[[assay]] <- sf
  for (route in c("genome", "txome")) {
    key <- paste(route, assay, sep = "_")
    scaled[[key]] <- sweep(raw[[key]], 2L, sf, "/")
  }
  cat(sprintf("  %-4s size factors: min %.4f  median %.4f  max %.4f\n",
              assay, min(sf), stats::median(sf), max(sf)))
}

# ── outputs ────────────────────────────────────────────────────────────────────────────────

out <- opts$output
dir.create(out, recursive = TRUE, showWarnings = FALSE)

for (route in c("genome", "txome")) {
  for (assay in c("ribo", "rna")) {
    m <- scaled[[paste(route, assay, sep = "_")]]
    frame <- data.frame(transcript_id = rownames(m), m, check.names = FALSE)
    utils::write.csv(frame, file.path(out, sprintf("%s_scaled_%s.csv", assay, route)),
                     row.names = FALSE, quote = FALSE)
  }
}
utils::write.csv(
  data.frame(sample = samples,
             ribo = as.numeric(size_factors$ribo[samples]),
             rna = as.numeric(size_factors$rna[samples])),
  file.path(out, "size_factors.csv"), row.names = FALSE, quote = FALSE)

cat(sprintf("wrote %s\n", out))
