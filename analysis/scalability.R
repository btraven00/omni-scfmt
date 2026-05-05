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
  files <- list.files(root, pattern = pat, recursive = TRUE, full.names = TRUE)
  files[grepl("/scale/scanpy_scale/", files, fixed = TRUE)]
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
if (!length(args)) args <- "out"

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

# Plot: stacked elapsed (load + compute + writing) per replicate, faceted
# by concurrency. The bar height shows total wall, the read+write segments
# show the I/O overhead this run is asking us to characterize.
long <- wide %>%
  select(run, concurrency, replicate_dir, load, compute, writing) %>%
  pivot_longer(cols = c(load, compute, writing),
               names_to = "stage", values_to = "elapsed_s") %>%
  mutate(stage = factor(stage, levels = c("load", "compute", "writing")))

dir.create("plots", showWarnings = FALSE)

p <- ggplot(long, aes(x = replicate_dir, y = elapsed_s, fill = stage)) +
  geom_col() +
  facet_wrap(~ paste0("N=", concurrency, "  (", run, ")"),
             scales = "free_x") +
  scale_fill_manual(values = c(load = "#4C72B0",
                               compute = "#55A868",
                               writing = "#C44E52")) +
  labs(
    title = "scanpy_scale: I/O vs compute under concurrent h5ad reads",
    x = "replicate (param hash)",
    y = "elapsed (s)",
    fill = NULL
  ) +
  theme_minimal(base_size = 11) +
  theme(axis.text.x = element_text(angle = 30, hjust = 1))

ggsave("plots/scalability_io_vs_compute.pdf", p,
       width = 7, height = 4)
cat("\nWrote plots/scalability_io_vs_compute.pdf\n")
