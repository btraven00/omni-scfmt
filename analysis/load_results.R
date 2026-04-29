#!/usr/bin/env Rscript
# Shared results loader for omni-scfmt analysis scripts.
#
# Walks an omnibenchmark output tree, finds every *.timing.tsv, attaches
# (dataset, module, tool, format, access_mode, converter, role) classifiers
# parsed from the module name, and (optionally) joins inner-denet RSS/CPU
# peaks. Returns one long-format tibble the plot scripts consume.
#
# Usage:
#   source("analysis/load_results.R")
#   results <- load_results("out")

suppressPackageStartupMessages({
  library(dplyr)
  library(jsonlite)
  library(tidyr)
})

# Parse a module name (e.g. "seurat_inmemory_rds", "scx_h5ad_to_h5seurat",
# "anndatar_h5ad_to_rds_sce") into a classifier row. Everything in the plot
# scripts keys off these columns, so keep it consistent.
classify_module <- function(module) {
  # Method modules:  {tool}_{access}_{format}
  if (grepl("^(seurat|scanpy)_", module)) {
    parts       <- strsplit(module, "_", fixed = TRUE)[[1]]
    tool        <- parts[1]
    access_mode <- parts[2]
    fmt         <- paste(parts[-c(1, 2)], collapse = "_")
    native <- switch(tool,
      seurat = fmt %in% c("rds", "bpcells"),
      scanpy = fmt %in% c("h5ad")
    )
    return(tibble::tibble(
      role        = "method",
      tool        = tool,
      access_mode = access_mode,
      format      = fmt,
      native      = native,
      converter   = NA_character_
    ))
  }

  # Convert modules: {converter}_{from}_to_{to}
  if (grepl("^(scx|anndatar)_", module)) {
    parts     <- strsplit(module, "_", fixed = TRUE)[[1]]
    converter <- parts[1]
    to_idx    <- which(parts == "to")
    from_fmt  <- paste(parts[2:(to_idx - 1)], collapse = "_")
    to_fmt    <- paste(parts[(to_idx + 1):length(parts)], collapse = "_")
    return(tibble::tibble(
      role        = "convert",
      tool        = NA_character_,
      access_mode = NA_character_,
      format      = paste(from_fmt, "->", to_fmt),
      native      = NA,
      converter   = converter
    ))
  }

  # Data modules: {dataset}
  tibble::tibble(
    role = "data", tool = NA_character_, access_mode = NA_character_,
    format = "h5ad", native = NA, converter = NA_character_
  )
}

# Extract dataset name from an omnibenchmark output path:
#   .../out/data/<dataset>/.../  →  <dataset>
extract_dataset <- function(path) {
  parts    <- strsplit(normalizePath(path, mustWork = FALSE),
                       .Platform$file.sep)[[1]]
  data_idx <- which(parts == "data")
  if (length(data_idx) && data_idx < length(parts)) parts[data_idx + 1] else NA_character_
}

# Load inner-denet JSONL (one sample per line) if present next to a timing
# file, and reduce it to a single-row peak/mean summary for the rule.
summarize_inner_denet <- function(rule_dir) {
  files <- list.files(rule_dir, pattern = "\\.inner_denet\\.jsonl$",
                      full.names = TRUE)
  if (!length(files)) {
    return(tibble::tibble(peak_rss_mb = NA_real_, mean_cpu_pct = NA_real_))
  }
  `%||%` <- function(a, b) if (!is.null(a)) a else b
  samples <- do.call(rbind, lapply(files, function(f) {
    lines <- readLines(f, warn = FALSE)
    if (length(lines) < 2) return(NULL)
    rows <- lapply(lines[-1], function(l) {
      tryCatch({
        d <- jsonlite::fromJSON(l)
        agg <- if (!is.null(d$aggregated)) d$aggregated else d
        data.frame(mem_rss_kb = agg$mem_rss_kb %||% NA_real_,
                   cpu_usage  = agg$cpu_usage  %||% NA_real_)
      }, error = function(e) NULL)
    })
    do.call(rbind, Filter(Negate(is.null), rows))
  }))
  if (is.null(samples) || !nrow(samples)) {
    return(tibble::tibble(peak_rss_mb = NA_real_, mean_cpu_pct = NA_real_))
  }
  rss  <- samples$mem_rss_kb / 1024
  cpu  <- samples$cpu_usage * 100
  tibble::tibble(
    peak_rss_mb  = suppressWarnings(max(rss,  na.rm = TRUE)),
    mean_cpu_pct = suppressWarnings(mean(cpu, na.rm = TRUE))
  )
}

load_results <- function(out_dir = "out") {
  tsv_files <- unique(normalizePath(list.files(
    out_dir, pattern = "\\.timing\\.tsv$",
    recursive = TRUE, full.names = TRUE, all.files = TRUE
  )))
  if (!length(tsv_files)) stop("No *.timing.tsv under: ", out_dir)

  do.call(bind_rows, lapply(tsv_files, function(f) {
    df      <- read.delim(f)
    module  <- basename(dirname(dirname(f)))
    dataset <- extract_dataset(f)
    klass   <- classify_module(module)
    inner   <- summarize_inner_denet(dirname(f))
    df |>
      mutate(module = module, dataset = dataset) |>
      bind_cols(klass[rep(1, nrow(df)), ], inner[rep(1, nrow(df)), ])
  }))
}

# When sourced from a plot script, don't run anything — just export loader.
invisible(NULL)
