#!/usr/bin/env Rscript
# Contribution 2: converter throughput and phase breakdown from omnibench-events.jsonl.
# Writes:
#   converter_throughput.png  — total elapsed per converter/pair, faceted by dataset
#   converter_phases.png      — load vs write breakdown per converter/pair
#
# Usage: Rscript plot_converter_compare.R [--out_dir out] [--plot_dir plots]

suppressPackageStartupMessages({
  library(argparse); library(dplyr); library(ggplot2); library(jsonlite)
})
parser <- ArgumentParser()
parser$add_argument("--out_dir",  default = "out")
parser$add_argument("--plot_dir", default = "plots")
args <- parser$parse_args()
dir.create(args$plot_dir, showWarnings = FALSE, recursive = TRUE)

# Parse all converter omnibench-events.jsonl files.
# Converter dirs sit at: out/data/<ds>/.default/convert_<pair>/<module>/.default/
# We identify them by: path contains a convert_* segment and no /method/ segment.
event_files <- list.files(args$out_dir, pattern = "omnibench-events\\.jsonl$",
                          recursive = TRUE, full.names = TRUE, all.files = TRUE)
event_files <- event_files[grepl("/convert[_/]", event_files) &
                           !grepl("/method/", event_files)]

if (!length(event_files)) {
  message("No converter omnibench-events.jsonl found — skipping converter plots.")
  quit(status = 0)
}

parse_events <- function(f) {
  lines <- readLines(f, warn = FALSE)
  rows  <- lapply(lines, function(l) tryCatch(fromJSON(l), error = function(e) NULL))
  rows  <- Filter(Negate(is.null), rows)
  ends  <- Filter(function(x) !is.null(x$phase) && x$phase == "end" &&
                               !is.null(x$attrs$elapsed_s), rows)
  if (!length(ends)) return(NULL)

  parts   <- strsplit(normalizePath(f), .Platform$file.sep)[[1]]
  data_i  <- which(parts == "data")
  dataset <- if (length(data_i)) parts[data_i + 1] else NA_character_

  # module dir is two levels up from the file (parent of .default)
  module  <- basename(dirname(dirname(f)))
  converter <- sub("_.*", "", module)

  do.call(rbind, lapply(ends, function(x) {
    data.frame(dataset = dataset, module = module, converter = converter,
               stage = x$event, elapsed_s = x$attrs$elapsed_s,
               stringsAsFactors = FALSE)
  }))
}

r <- do.call(rbind, lapply(event_files, parse_events))
if (is.null(r) || !nrow(r)) {
  message("Converter events parsed but empty — skipping converter plots.")
  quit(status = 0)
}

n_cells <- function(d) {
  m <- regmatches(d, regexpr("\\d+[km]?", d, ignore.case = TRUE))
  if (!length(m) || !nzchar(m)) return(NA_real_)
  if (endsWith(m, "k")) as.numeric(sub("k$", "", m)) * 1e3
  else if (endsWith(m, "m")) as.numeric(sub("m$", "", m)) * 1e6
  else as.numeric(m)
}

totals <- r |>
  group_by(dataset, module, converter) |>
  summarise(total_s = sum(elapsed_s), .groups = "drop") |>
  mutate(n = vapply(dataset, n_cells, numeric(1)),
         cells_per_s = n / total_s)

# ---- throughput (total elapsed) ---------------------------------------------
p_thru <- ggplot(totals, aes(x = module, y = total_s, fill = converter)) +
  geom_col(position = position_dodge(0.8), width = 0.75) +
  facet_wrap(~ dataset, scales = "free_y") +
  labs(title = "Converter total elapsed time", x = NULL, y = "Elapsed (s)",
       fill = "Converter") +
  theme_bw() +
  theme(axis.text.x = element_text(angle = 30, hjust = 1))
ggsave(file.path(args$plot_dir, "converter_throughput.png"),
       p_thru, width = 9, height = 5, dpi = 150)
message("Wrote: ", file.path(args$plot_dir, "converter_throughput.png"))

# ---- load vs write phase breakdown ------------------------------------------
phases <- r |>
  filter(stage %in% c("load", "write")) |>
  mutate(stage = toupper(stage),
         stage = factor(stage, levels = c("WRITE", "LOAD")))

p_phases <- ggplot(phases, aes(x = module, y = elapsed_s, fill = stage)) +
  geom_col(position = position_stack(reverse = TRUE)) +
  facet_wrap(~ dataset, scales = "free_y") +
  coord_flip() +
  labs(title = "Converter phase breakdown (LOAD vs WRITE)",
       x = NULL, y = "Elapsed (s)", fill = "Phase") +
  theme_bw()
ggsave(file.path(args$plot_dir, "converter_phases.png"),
       p_phases, width = 9, height = 5, dpi = 150)
message("Wrote: ", file.path(args$plot_dir, "converter_phases.png"))

message("Converter plots written to: ", args$plot_dir)
