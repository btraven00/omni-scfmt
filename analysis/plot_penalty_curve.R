#!/usr/bin/env Rscript
# Contribution 1 primary figure: foreign/native wall-time ratio vs N cells,
# one line per (tool, foreign_format). Log-log axes.
# Usage: Rscript plot_penalty_curve.R [--out_dir out] [--plot_dir plots]

suppressPackageStartupMessages({
  library(argparse); library(dplyr); library(ggplot2)
})
parser <- ArgumentParser()
parser$add_argument("--out_dir",  default = "out")
parser$add_argument("--plot_dir", default = "plots")
args <- parser$parse_args()
dir.create(args$plot_dir, showWarnings = FALSE, recursive = TRUE)

HERE <- dirname(normalizePath(sub("--file=", "", grep("--file=", commandArgs(FALSE), value = TRUE)[1])))
source(file.path(HERE, "load_results.R"))
r <- load_results(args$out_dir) |>
  filter(role == "method", access_mode == "inmemory",
         stage %in% c("load", "access", "write"))

# Map each (tool, dataset) to its native-format total time; then compute
# ratio for every foreign-format run.
totals <- r |>
  group_by(tool, format, native, dataset) |>
  summarise(total_s = sum(elapsed_s), .groups = "drop")
natives <- totals |> filter(native) |>
  select(tool, dataset, native_s = total_s)
penalty <- totals |> filter(!native) |>
  inner_join(natives, by = c("tool", "dataset")) |>
  mutate(ratio = total_s / native_s)

# N cells per dataset — parse from dataset name when possible, else fall
# back to a small lookup so the plot still works for pbmc3k/pbmc68k.
n_cells <- function(d) {
  m <- regmatches(d, regexpr("\\d+[km]?", d))
  if (!length(m) || !nzchar(m)) return(NA_real_)
  if (endsWith(m, "k")) as.numeric(sub("k$", "", m)) * 1e3
  else if (endsWith(m, "m")) as.numeric(sub("m$", "", m)) * 1e6
  else as.numeric(m)
}
penalty$n_cells <- vapply(penalty$dataset, n_cells, numeric(1))

p <- ggplot(penalty, aes(x = n_cells, y = ratio, colour = format,
                          linetype = tool, shape = tool)) +
  geom_hline(yintercept = 1, linetype = "dashed", colour = "grey50") +
  geom_line() + geom_point(size = 2.5) +
  scale_x_log10() + scale_y_log10() +
  labs(
    title  = "Cross-tool format penalty",
    x      = "N cells (log)",
    y      = "Foreign / native wall-time ratio (log)",
    colour = "Foreign format",
    linetype = "Tool", shape = "Tool"
  ) +
  theme_bw()
ggsave(file.path(args$plot_dir, "penalty_curve.png"), p,
       width = 7, height = 5, dpi = 150)
message("Wrote: ", file.path(args$plot_dir, "penalty_curve.png"))
