#!/usr/bin/env Rscript
# Ward (ward.D2) hierarchical clustering of the genes by read-fate composition.
#
#   Rscript code/clustering/ward_cluster.R --input <stem>.gene_counts_filtered.tsv \
#       --output results/clustering --stem <stem> [--k 4] [--k-min 2] [--k-max 8] \
#       [--seed 1] [--silhouette-sample 4000]
#
# Each gene is its five fate counts divided by n_union (a composition), the distance is
# Euclidean between those shares, the tree is stats::hclust(method = "ward.D2"). For every
# k in k-min..k-max the cut's within-SS, between/total SS, average silhouette (on a seeded
# subsample: the full distance matrix at ~12,000 genes is ~1 GB for nothing) and the merge
# height that turns k clusters into k - 1 are recorded. The cut at --k is then relabelled so
# that cluster 1 has the most concordant centroid (highest shared/unique share) and written
# per gene with the composition in %, plus one centroid row per cluster (member mean).
#
# Writes, under --output, all named by the stem:
#   <stem>.tree.rds              the hclust object
#   <stem>.tree_merge.tsv        its merge matrix and heights (what the panel draws)
#   <stem>.tree.tsv              linkage, n, cophenetic correlation
#   <stem>.hclust_sweep.tsv      one row per k, the chosen one marked
#   <stem>.clusters_k<K>.tsv     per gene: cluster and the five shares in %
#   <stem>.cluster_centroids.tsv per cluster: n and the centroid composition in %
#
# Base R only; `here` is derived from Rscript's --file= so it runs from any directory.

options(digits = 15)   # so the heights round-trip through the TSV unchanged

# ── arguments ──────────────────────────────────────────────────────────────────────────
here <- dirname(normalizePath(sub("^--file=", "",
                                  grep("^--file=", commandArgs(FALSE), value = TRUE))[1]))
root <- dirname(dirname(here))
source(file.path(root, "code", "common", "cli_args.R"))
log_line <- function(fmt, ...) cat(sprintf("[clustering/ward] %s\n", sprintf(fmt, ...)))

# ── the five fates ─────────────────────────────────────────────────────────────────────
#: Column order is load-bearing: index 1 is the concordant fate, which is how clusters
#: get relabelled (cluster 1 = most concordant). Same order as the count table.
COMPONENTS <- c("shared_genome_unique", "shared_genome_multimapped",
                "genome_only_unique", "genome_only_multimapped", "transcriptome_only")
ID_COLS <- c("gene", "gene_id", "transcript_id")
n_cols <- paste0("n_", COMPONENTS)
pct_cols <- paste0("pct_", COMPONENTS)

read_tsv <- function(path, required = character()) {
  if (!file.exists(path)) stop("no such file: ", path)
  d <- utils::read.delim(path, stringsAsFactors = FALSE, check.names = FALSE)
  missing <- setdiff(required, names(d))
  if (length(missing)) stop(basename(path), " lacks column(s): ", paste(missing, collapse = ", "))
  d
}
write_tsv <- function(frame, path) {
  dir.create(dirname(path), showWarnings = FALSE, recursive = TRUE)
  utils::write.table(frame, path, sep = "\t", quote = FALSE, row.names = FALSE)
  log_line("wrote %s", path)
  invisible(path)
}

# Average silhouette on a seeded subsample, Euclidean in X (the geometry Ward minimises).
silhouette_mean <- function(X, cluster, n_sample, seed) {
  set.seed(seed)
  idx <- if (nrow(X) > n_sample) sample.int(nrow(X), n_sample) else seq_len(nrow(X))
  lab <- cluster[idx]
  if (length(unique(lab)) < 2) return(NA_real_)
  dm <- as.matrix(stats::dist(X[idx, , drop = FALSE], method = "euclidean"))
  sizes <- table(lab)
  s <- vapply(seq_along(lab), function(i) {
    own <- lab[i]
    if (sizes[[as.character(own)]] < 2) return(0)
    means <- vapply(split(dm[i, -i], lab[-i]), mean, numeric(1))
    a <- means[[as.character(own)]]
    b <- min(means[names(means) != as.character(own)])
    (b - a) / max(a, b)
  }, numeric(1))
  mean(s)
}

# ── run ────────────────────────────────────────────────────────────────────────────────
opts <- parse_args(list(input = "", output = "", stem = "", k = "4", `k-min` = "2",
                        `k-max` = "8", seed = "1", `silhouette-sample` = "4000"))
input <- require_arg(opts, "input"); out_dir <- require_arg(opts, "output")
stem <- require_arg(opts, "stem")
chosen <- as.integer(opts$k)
k_range <- as.integer(opts[["k-min"]]):as.integer(opts[["k-max"]])
seed <- as.integer(opts$seed); sil_n <- as.integer(opts[["silhouette-sample"]])
if (!chosen %in% k_range) stop("--k must lie in --k-min..--k-max")

d <- read_tsv(input, c(ID_COLS, "n_union", n_cols))
counts <- as.matrix(d[, n_cols]); storage.mode(counts) <- "double"
if (!all(rowSums(counts) == d$n_union)) stop("the five counts do not sum to n_union")
if (any(counts < 0)) stop("negative count")
X <- counts / rowSums(counts)
stopifnot(all(abs(rowSums(X) - 1) < 1e-12), all(X >= 0))
log_line("%s: %d genes x %d shares; ward.D2, k = %d..%d, cut at k = %d",
         basename(input), nrow(X), ncol(X), min(k_range), max(k_range), chosen)

D <- stats::dist(X, method = "euclidean")
tree <- stats::hclust(D, method = "ward.D2")
coph <- stats::cor(D, stats::cophenetic(tree))
rm(D)
log_line("cophenetic correlation %.4f", coph)

dir.create(out_dir, showWarnings = FALSE, recursive = TRUE)
path <- function(suffix) file.path(out_dir, sprintf("%s.%s", stem, suffix))
saveRDS(tree, path("tree.rds"))
write_tsv(data.frame(merge1 = tree$merge[, 1], merge2 = tree$merge[, 2], height = tree$height),
          path("tree_merge.tsv"))
write_tsv(data.frame(linkage = "ward.D2", n = nrow(X), cophenetic_correlation = coph),
          path("tree.tsv"))

totss <- sum(scale(X, scale = FALSE)^2)
rows <- list()
for (k in k_range) {
  cluster <- stats::cutree(tree, k = k)
  centres <- t(sapply(seq_len(k), function(i) colMeans(X[cluster == i, , drop = FALSE])))
  sizes <- as.integer(table(factor(cluster, levels = seq_len(k))))
  withinss <- sum((X - centres[cluster, , drop = FALSE])^2)
  sil <- silhouette_mean(X, cluster, sil_n, seed)
  # heights are increasing; the (n - k + 1)-th merge turns the k-cluster cut into k - 1
  height <- tree$height[nrow(X) - k + 1]
  rows[[length(rows) + 1]] <- data.frame(
    k = k, tot_withinss = withinss, between_over_total = 1 - withinss / totss,
    avg_silhouette = sil, smallest_cluster = min(sizes),
    sizes = paste(sizes, collapse = "/"), merge_height = height, chosen = k == chosen,
    stringsAsFactors = FALSE)
  log_line("k = %d  withinSS %9.2f  between/total %.4f  silhouette %6.3f  height %.3f  sizes %s",
           k, withinss, 1 - withinss / totss, sil, height, paste(sizes, collapse = "/"))
  if (k == chosen) chosen_cut <- list(cluster = cluster, centres = centres, sizes = sizes)
}
write_tsv(do.call(rbind, rows), path("hclust_sweep.tsv"))

# the chosen cut: cluster 1 = the most concordant centroid, 2 the next, ...
centres_P <- chosen_cut$centres / rowSums(chosen_cut$centres)
colnames(centres_P) <- COMPONENTS
ord <- order(-centres_P[, 1])
cluster <- match(chosen_cut$cluster, ord)
centres_P <- centres_P[ord, , drop = FALSE]
sizes <- as.integer(table(factor(cluster, levels = seq_len(chosen))))

log_line("centroid composition (%%) at k = %d:", chosen)
for (i in seq_len(chosen))
  cat(sprintf("  %d  n = %6d  %s\n", i, sizes[i],
              paste(sprintf("%5.1f", 100 * centres_P[i, ]), collapse = "  ")))

pct <- 100 * counts / d$n_union
colnames(pct) <- pct_cols
write_tsv(data.frame(d[, c(ID_COLS, "n_union")], cluster = cluster, round(pct, 6),
                     check.names = FALSE),
          path(sprintf("clusters_k%d.tsv", chosen)))
write_tsv(data.frame(cluster = seq_len(chosen), n = sizes,
                     `colnames<-`(round(100 * centres_P, 4), pct_cols), check.names = FALSE),
          path("cluster_centroids.tsv"))
