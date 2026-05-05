#!/usr/bin/env Rscript
# Scalability plot: I/O overhead (load + writing) vs compute, across the
# replicates of a scanpy_scale run. Reads obkit logger events under the
# scale stage, one event log per replicate.
#
# Concurrency level = number of replicate dirs in the run (the parameter
# expansion in scalability.yaml). To compare across N, run the benchmark
# at different replicate counts and pass each `out/` path:
#
#   Rscript analysis/scalability.R out/                # single run
#   Rscript analysis/scalability.R out_N3 out_N8 ...   # multiple runs
#
# Outputs:
#   plots/scalability_io_vs_compute.pdf

suppressPackageStartupMessages({
  library(dplyr)
  library(tidyr)
  library(ggplot2)
  library(obkit)        # read_phase_ranges()
})

# obkit's logger filename changed from "omnibench-events.jsonl" (<= 0.0.2)
# to "obkit-events.jsonl" on main. Match either.
EVENTS_GLOB <- c("omnibench-events.jsonl", "obkit-events.jsonl")

find_event_files <- function(root) {
  pat <- paste(EVENTS_GLOB, collapse = "|")
  # all.files=TRUE because ob writes outputs into dot-prefixed dirs
  # (e.g. `.default/`, `.0546af9d/`) which list.files() skips by default.
  files <- list.files(root, pattern = pat, recursive = TRUE,
                      full.names = TRUE, all.files = TRUE)
  files <- files[grepl("/scale/scanpy_scale/", files, fixed = TRUE)]
  # ob writes both a hash dir (`.0546af9d`) and a human-readable symlink
  # (`replicate-1`) — same file reachable via multiple paths. Dedupe.
  unique(normalizePath(files))
}

load_run <- function(root) {
  files <- find_event_files(root)
  if (!length(files)) stop("No scanpy_scale event files under ", root)

  rows <- lapply(files, function(f) {
    df <- obkit::read_phase_ranges(f)
    df$elapsed_s <- (df$xmax - df$xmin) / 1000
    df$replicate_dir <- basename(dirname(f))
    df
  })
  out <- do.call(rbind, rows)
  out$run <- basename(normalizePath(root))
  out$concurrency <- length(unique(out$replicate_dir))
  out
}

args <- commandArgs(trailingOnly = TRUE)
if (!length(args)) {
  # Default: every runs/N* directory, sorted numerically.
  candidates <- list.dirs("runs", full.names = TRUE, recursive = FALSE)
  candidates <- candidates[grepl("/N[0-9]+$", candidates)]
  ord <- order(as.integer(sub(".*/N", "", candidates)))
  args <- candidates[ord]
  if (!length(args)) stop("No runs/N* directories found; pass paths explicitly.")
}
cat("Reading runs:", paste(args, collapse = ", "), "\n")

events <- do.call(rbind, lapply(args, load_run))

# Wide form: one row per (run, replicate), columns = stages.
wide <- events %>%
  select(run, concurrency, replicate_dir, event, elapsed_s) %>%
  pivot_wider(names_from = event, values_from = elapsed_s) %>%
  mutate(
    io_s      = load + writing,
    compute_s = compute,
    io_frac   = io_s / (io_s + compute_s)
  )

cat("\nPer-replicate timings:\n")
print(as.data.frame(wide))

cat("\nPer-concurrency summary (mean ± sd):\n")
summary_tbl <- wide %>%
  group_by(run, concurrency) %>%
  summarise(
    load_mean    = mean(load),    load_sd    = sd(load),
    compute_mean = mean(compute), compute_sd = sd(compute),
    write_mean   = mean(writing), write_sd   = sd(writing),
    io_frac_mean = mean(io_frac),
    .groups = "drop"
  )
print(as.data.frame(summary_tbl))

# Aggregate across replicates within each (run, concurrency, stage):
# mean elapsed time and sd for the error bars.
long <- wide %>%
  select(run, concurrency, replicate_dir, load, compute, writing) %>%
  pivot_longer(cols = c(load, compute, writing),
               names_to = "stage", values_to = "elapsed_s") %>%
  mutate(stage = factor(stage, levels = c("load", "compute", "writing")))

agg <- long %>%
  group_by(concurrency, stage) %>%
  summarise(
    mean_s = mean(elapsed_s),
    sd_s   = sd(elapsed_s),
    n      = dplyr::n(),
    .groups = "drop"
  ) %>%
  mutate(
    # Cap the error bar at zero on the log scale: clip lower bound away
    # from <=0, otherwise log10 fails.
    ymin = pmax(mean_s - sd_s, mean_s * 0.5, 1e-4),
    ymax = mean_s + ifelse(is.na(sd_s), 0, sd_s)
  )

dir.create("plots", showWarnings = FALSE)

p <- ggplot(agg, aes(x = concurrency, y = mean_s, color = stage, group = stage)) +
  geom_line(linewidth = 0.7) +
  geom_point(size = 2.5) +
  geom_errorbar(aes(ymin = ymin, ymax = ymax), width = 0.15, linewidth = 0.5) +
  scale_x_continuous(breaks = sort(unique(agg$concurrency))) +
  scale_y_log10() +
  scale_color_manual(values = c(load = "#4C72B0",
                                compute = "#55A868",
                                writing = "#C44E52")) +
  labs(
    title = "scanpy_scale: stage time vs concurrency",
    subtitle = "mean ± sd across replicates within each N",
    x = "concurrency level (N)",
    y = "elapsed (s, log scale)",
    color = NULL
  ) +
  theme_minimal(base_size = 11)

ggsave("plots/scalability_io_vs_compute.pdf", p,
       width = 6, height = 4)
cat("\nWrote plots/scalability_io_vs_compute.pdf\n")
