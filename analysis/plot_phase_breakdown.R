#!/usr/bin/env Rscript
# Contribution 1 supporting figure: stacked wall-time per (tool, format)
# showing where the cross-tool penalty lands (LOAD vs ACCESS vs WRITE vs
# compute). One panel per dataset.
# Usage: Rscript plot_phase_breakdown.R [--out_dir out] [--plot_dir plots]

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
all_r <- load_results(args$out_dir) |> filter(role == "method")

# peak RSS per (dataset, module) — one value per bar
rss_labels <- all_r |>
  group_by(dataset, module, tool, format, access_mode) |>
  summarise(peak_rss_mb = first(peak_rss_mb), .groups = "drop") |>
  mutate(cell = paste(tool, access_mode, format, sep = "/"),
         rss_label = sprintf("%.0f MB", peak_rss_mb))

r <- all_r |>
  mutate(phase = case_when(
    stage == "load"      ~ "LOAD",
    stage == "access"    ~ "ACCESS",
    stage == "write"     ~ "WRITE",
    stage == "to_memory" ~ "TO_MEMORY",
    TRUE                 ~ "compute"
  )) |>
  group_by(dataset, tool, format, access_mode, phase) |>
  summarise(elapsed_s = sum(elapsed_s), .groups = "drop")

r$cell <- paste(r$tool, r$access_mode, r$format, sep = "/")
# ggplot2 4.x stacks the last data row at xmin=0 (leftmost). Reverse the
# factor levels so compute sorts first and LOAD sorts last → LOAD lands left.
r$phase <- factor(r$phase, levels = c("compute", "TO_MEMORY", "WRITE", "ACCESS", "LOAD"))
r <- r |> arrange(dataset, cell, phase)

# total elapsed per bar — for positioning the RSS label just past the bar end
bar_totals <- r |>
  group_by(dataset, cell) |>
  summarise(total_s = sum(elapsed_s), .groups = "drop") |>
  left_join(rss_labels, by = c("dataset", "cell"))

for (ds in unique(r$dataset)) {
  p <- ggplot(filter(r, dataset == ds),
              aes(y = cell, x = elapsed_s, fill = phase)) +
    geom_col() +
    geom_text(data = filter(bar_totals, dataset == ds),
              aes(y = cell, x = total_s, label = rss_label),
              inherit.aes = FALSE, hjust = -0.1, size = 3) +
    scale_x_continuous(expand = expansion(mult = c(0, 0.2))) +
    guides(fill = guide_legend(reverse = TRUE)) +
    labs(title = paste("Phase breakdown —", ds),
         y = NULL, x = "Elapsed (s)", fill = "Phase") +
    theme_bw()
  out <- file.path(args$plot_dir, paste0("phase_breakdown_", ds, ".png"))
  ggsave(out, p, width = 8, height = 4 + 0.25 * length(unique(r$cell)),
         dpi = 150)
  message("Wrote: ", out)
}
