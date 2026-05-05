#!/usr/bin/env Rscript
# CPU+RSS time-series profiles per concurrency level. For each `runs/N*`
# directory, picks one representative replicate, calls
# obkit::plot_annotated_profiles() (CPU on left axis, memory on right,
# phase rectangles overlaid), and tiles the resulting plots into one PDF.
#
# Usage:
#   Rscript analysis/scalability_profiles.R               # auto-detect runs/N*
#   Rscript analysis/scalability_profiles.R runs/N4 ...   # explicit paths
#
# Output: plots/scalability_profiles.pdf

suppressPackageStartupMessages({
  library(obkit)
  library(ggplot2)
  library(patchwork)
})

args <- commandArgs(trailingOnly = TRUE)
if (!length(args)) {
  candidates <- list.dirs("runs", full.names = TRUE, recursive = FALSE)
  candidates <- candidates[grepl("/N[0-9]+$", candidates)]
  ord <- order(as.integer(sub(".*/N", "", candidates)))
  args <- candidates[ord]
  if (!length(args)) stop("No runs/N* directories found; pass paths explicitly.")
}
cat("Reading runs:", paste(args, collapse = ", "), "\n")

# Pick one replicate dir per run. Prefer rep 1 (lowest pid → earliest
# launch under most schedulers) for representativeness; fall back to first.
pick_rep_dir <- function(run_dir) {
  events <- list.files(run_dir,
                       pattern = "obkit-events.jsonl|omnibench-events.jsonl",
                       recursive = TRUE, full.names = TRUE, all.files = TRUE)
  events <- events[grepl("/scale/scanpy_scale/", events, fixed = TRUE)]
  events <- unique(normalizePath(events))
  if (!length(events)) return(NULL)
  # Sort by replicate number (parsed from sibling denet filename).
  rep_dirs <- dirname(events)
  reps <- vapply(rep_dirs, function(d) {
    f <- list.files(d, pattern = "scanpy_scale\\.rep[0-9]+\\.denet\\.jsonl",
                    full.names = TRUE)
    if (!length(f)) return(NA_integer_)
    as.integer(sub(".*\\.rep([0-9]+)\\.denet\\.jsonl$", "\\1", f[1]))
  }, integer(1))
  rep_dirs[order(reps, na.last = TRUE)][1]
}

build_panel <- function(run_dir) {
  rep_dir <- pick_rep_dir(run_dir)
  if (is.null(rep_dir)) {
    warning(run_dir, ": no scanpy_scale replicate found")
    return(NULL)
  }
  denet <- list.files(rep_dir, pattern = "\\.denet\\.jsonl$", full.names = TRUE)[1]
  evt   <- list.files(rep_dir,
                      pattern = "obkit-events.jsonl|omnibench-events.jsonl",
                      full.names = TRUE, all.files = TRUE)[1]

  samples <- obkit::read_aggregates(denet)
  # Workaround for obkit 0.0.3: read_aggregates() returns columns named
  # `ts_ms`/`cpu_usage`/..., but plot_annotated_profiles() expects them
  # prefixed with `aggregated.`. Rename here until upstream fixes it.
  names(samples) <- paste0("aggregated.", names(samples))
  ranges  <- obkit::read_phase_ranges(evt)
  if (!nrow(samples)) return(NULL)

  N <- basename(run_dir)
  mx <- max(110, ceiling(max(samples$aggregated.cpu_usage) / 100) * 100)
  obkit::plot_annotated_profiles(
    samples = samples,
    ranges  = ranges,
    title   = paste0(N, " (rep 1)"),
    mx      = mx
  )
}

panels <- Filter(Negate(is.null), lapply(args, build_panel))
if (!length(panels)) stop("No panels could be built — check `runs/`.")

# Assemble into a tidy grid. Patchwork lays out into rows of <=3 panels.
ncol <- min(3, length(panels))
nrow <- ceiling(length(panels) / ncol)
combined <- patchwork::wrap_plots(panels, ncol = ncol) +
  patchwork::plot_annotation(
    title = "scanpy_scale: CPU% (left) + RSS MB (right) with phase overlay",
    subtitle = "one representative replicate per concurrency level"
  )

dir.create("plots", showWarnings = FALSE)
ggsave("plots/scalability_profiles.pdf", combined,
       width = 5 * ncol, height = 3.2 * nrow, limitsize = FALSE)
cat("\nWrote plots/scalability_profiles.pdf (", length(panels), "panels)\n")
