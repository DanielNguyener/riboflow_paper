# Shared --key value argument parsing for the base-R analysis scripts.
# Sourced with each script's own `here`; base R only.

parse_args <- function(defaults) {
  raw <- commandArgs(TRUE)
  if (length(raw) && raw[1] %in% c("--help", "-h")) {
    cat("usage: Rscript", basename(sub("^--file=", "",
        grep("^--file=", commandArgs(FALSE), value = TRUE))[1]),
        paste(sprintf("[--%s VALUE]", names(defaults)), collapse = " "), "\n")
    for (k in names(defaults)) cat(sprintf("  --%-18s default %s\n", k,
                                           if (nzchar(defaults[[k]])) defaults[[k]] else "(none)"))
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

require_arg <- function(opts, key) {
  if (!nzchar(opts[[key]])) stop("--", key, " is required")
  opts[[key]]
}
