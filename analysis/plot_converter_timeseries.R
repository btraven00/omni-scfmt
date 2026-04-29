#!/usr/bin/env Rscript
# CPU and RSS time series for every converter run, one PNG per dataset.
# Reads *.inner_denet.jsonl from converter output dirs (not method dirs).
# Usage: Rscript plot_converter_timeseries.R [--out_dir out] [--plot_dir plots]

suppressPackageStartupMessages({
  library(argparse); library(dplyr); library(ggplot2); library(jsonlite)
  library(tidyr)
})
parser <- ArgumentParser()
parser$add_argument("--out_dir",  default = "out")
parser$add_argument("--plot_dir", default = "plots")
args <- parser$parse_args()
dir.create(args$plot_dir, showWarnings = FALSE, recursive = TRUE)

denet_files <- list.files(
  args$out_dir, pattern = "\\.inner_denet\\.jsonl$",
  recursive = TRUE, full.names = TRUE, all.files = TRUE
)
# keep only converter dirs (contain /convert_/ but not /method/)
denet_files <- denet_files[grepl("/convert[_/]", denet_files) &
                           !grepl("/method/", denet_files)]

if (!length(denet_files)) {
  message("No converter inner_denet files found.")
  quit(status = 0)
}

read_ts <- function(f) {
  lines <- readLines(f, warn = FALSE)
  # first line is header (pid/cmd/t0_ms), rest are samples
  if (length(lines) < 2) return(NULL)
  t0 <- tryCatch(fromJSON(lines[1])$t0_ms, error = function(e) NA_real_)
  rows <- lapply(lines[-1], function(l) {
    tryCatch({
      d <- fromJSON(l)
      agg <- if (!is.null(d$aggregated)) d$aggregated else d
      data.frame(
        ts_ms      = agg$ts_ms,
        cpu_pct    = agg$cpu_usage * 100,
        rss_mb     = agg$mem_rss_kb / 1024
      )
    }, error = function(e) NULL)
  })
  rows <- do.call(rbind, Filter(Negate(is.null), rows))
  if (is.null(rows) || !nrow(rows)) return(NULL)

  parts    <- strsplit(normalizePath(f), .Platform$file.sep)[[1]]
  data_idx <- which(parts == "data")
  dataset  <- if (length(data_idx)) parts[data_idx + 1] else NA_character_
  module   <- basename(dirname(dirname(f)))  # parent of .default

  rows$t_s     <- (rows$ts_ms - min(rows$ts_ms, na.rm = TRUE)) / 1000
  rows$dataset <- dataset
  rows$module  <- module
  rows
}

ts <- do.call(rbind, lapply(denet_files, read_ts))
if (is.null(ts) || !nrow(ts)) {
  message("No time series data parsed.")
  quit(status = 0)
}

# long format for faceting on metric
ts_long <- ts |>
  pivot_longer(c(cpu_pct, rss_mb), names_to = "metric", values_to = "value") |>
  mutate(metric = recode(metric, cpu_pct = "CPU (%)", rss_mb = "RSS (MB)"))

for (ds in unique(ts_long$dataset)) {
  sub <- filter(ts_long, dataset == ds)
  p <- ggplot(sub, aes(x = t_s, y = value, colour = module)) +
    geom_line(linewidth = 0.7) +
    facet_wrap(~ metric, ncol = 1, scales = "free_y") +
    labs(title = paste("Converter CPU & RSS —", ds),
         x = "Elapsed (s)", y = NULL, colour = "Module") +
    theme_bw() +
    theme(legend.position = "bottom",
          legend.text = element_text(size = 7))
  out <- file.path(args$plot_dir, paste0("converter_timeseries_", ds, ".png"))
  ggsave(out, p, width = 8, height = 6, dpi = 150)
  message("Wrote: ", out)
}
