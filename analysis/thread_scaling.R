#!/usr/bin/env Rscript
# Thread-scaling speedup plot.  Reads obkit logger events from the
# scale/scanpy_thread_scale output directories and plots wall time and
# speedup (relative to threads=1) per stage.
#
# Usage:
#   Rscript analysis/thread_scaling.R              # reads out/ in cwd
#   Rscript analysis/thread_scaling.R /path/to/out
#
# Outputs:
#   plots/thread_scaling_time.pdf    — wall time per stage vs thread count
#   plots/thread_scaling_speedup.pdf — speedup (T1 / Tn) per stage

suppressPackageStartupMessages({
  library(dplyr)
  library(tidyr)
  library(ggplot2)
  library(obkit)
})

EVENTS_GLOB <- c("omnibench-events.jsonl", "obkit-events.jsonl")

find_event_files <- function(root) {
  pat <- paste(EVENTS_GLOB, collapse = "|")
  files <- list.files(root, pattern = pat, recursive = TRUE,
                      full.names = TRUE, all.files = TRUE)
  files <- files[grepl("/scale/scanpy_thread_scale/", files, fixed = TRUE)]
  unique(normalizePath(files))
}

# Extract the threads value from the parameter dir name or parameters.json.
threads_from_dir <- function(path) {
  # ob names the human-readable symlink "threads-N"
  dirs <- strsplit(path, .Platform$file.sep)[[1]]
  m <- regmatches(dirs, regexpr("^threads-([0-9]+)$", dirs))
  m <- m[nchar(m) > 0]
  if (length(m)) return(as.integer(sub("threads-", "", m[1])))
  # fallback: read parameters.json
  pjson <- file.path(dirname(path), "parameters.json")
  if (file.exists(pjson)) {
    p <- jsonlite::fromJSON(pjson)
    if (!is.null(p$threads)) return(as.integer(p$threads))
  }
  NA_integer_
}

args <- commandArgs(trailingOnly = TRUE)
if (!length(args)) args <- "out"
root <- args[1]

files <- find_event_files(root)
if (!length(files)) stop("No scanpy_thread_scale event files under ", root)
cat("Found", length(files), "event file(s)\n")

rows <- lapply(files, function(f) {
  df <- obkit::read_phase_ranges(f)
  df$elapsed_s <- (df$xmax - df$xmin) / 1000
  df$threads <- threads_from_dir(f)
  df
})
events <- do.call(rbind, rows)

wide <- events %>%
  select(threads, event, elapsed_s) %>%
  pivot_wider(names_from = event, values_from = elapsed_s) %>%
  arrange(threads)

cat("\nRaw timings (s):\n")
print(as.data.frame(wide))

# Speedup relative to threads=1
baseline <- wide %>% filter(threads == 1) %>%
  select(load_1 = load, compute_1 = compute, writing_1 = writing)

speedup <- wide %>%
  mutate(
    load_speedup    = baseline$load_1    / load,
    compute_speedup = baseline$compute_1 / compute,
    writing_speedup = baseline$writing_1 / writing,
  )

cat("\nSpeedup vs 1 thread:\n")
print(as.data.frame(speedup %>% select(threads, ends_with("speedup"))))

long_time <- wide %>%
  pivot_longer(c(load, compute, writing), names_to = "stage", values_to = "elapsed_s") %>%
  mutate(stage = factor(stage, levels = c("load", "compute", "writing")))

long_speedup <- speedup %>%
  select(threads, load_speedup, compute_speedup, writing_speedup) %>%
  pivot_longer(-threads, names_to = "stage", values_to = "speedup") %>%
  mutate(stage = sub("_speedup", "", stage),
         stage = factor(stage, levels = c("load", "compute", "writing")))

stage_colors <- c(load = "#4C72B0", compute = "#55A868", writing = "#C44E52")

dir.create("plots", showWarnings = FALSE)

p_time <- ggplot(long_time, aes(x = threads, y = elapsed_s, color = stage, group = stage)) +
  geom_line(linewidth = 0.8) +
  geom_point(size = 2.5) +
  scale_x_continuous(breaks = unique(long_time$threads)) +
  scale_color_manual(values = stage_colors) +
  labs(title = "scanpy thread scaling: wall time per stage",
       x = "threads", y = "elapsed (s)", color = NULL) +
  theme_minimal(base_size = 11)

p_speedup <- ggplot(long_speedup, aes(x = threads, y = speedup, color = stage, group = stage)) +
  geom_line(linewidth = 0.8) +
  geom_point(size = 2.5) +
  geom_abline(slope = 1, intercept = 0, linetype = "dashed", color = "grey60") +
  scale_x_continuous(breaks = unique(long_speedup$threads)) +
  scale_color_manual(values = stage_colors) +
  labs(title = "scanpy thread scaling: speedup vs 1 thread",
       subtitle = "dashed = ideal linear speedup",
       x = "threads", y = "speedup (T1 / Tn)", color = NULL) +
  theme_minimal(base_size = 11)

ggsave("plots/thread_scaling_time.pdf",    p_time,    width = 6, height = 4)
ggsave("plots/thread_scaling_speedup.pdf", p_speedup, width = 6, height = 4)
cat("\nWrote plots/thread_scaling_time.pdf and plots/thread_scaling_speedup.pdf\n")
