#!/usr/bin/env Rscript
# Thread-scaling speedup plot comparing h5ad (HDF5) vs zarr (blosc) backends.
# Reads obkit logger events from scale/scanpy_thread_scale* output dirs.
#
# Usage:
#   Rscript analysis/thread_scaling.R              # reads out/ in cwd
#   Rscript analysis/thread_scaling.R /path/to/out
#
# Outputs:
#   plots/thread_scaling_time.pdf    — wall time per stage and backend
#   plots/thread_scaling_speedup.pdf — speedup (T1 / Tn) per stage and backend

suppressPackageStartupMessages({
  library(dplyr)
  library(tidyr)
  library(ggplot2)
  library(jsonlite)
  library(obkit)
})

EVENTS_GLOB <- c("omnibench-events.jsonl", "obkit-events.jsonl")

BACKENDS <- list(
  scanpy_thread_scale      = "HDF5 (h5ad)",
  scanpy_thread_scale_zarr = "zarr (blosc)"
)

find_event_files <- function(root, module_id) {
  pat <- paste(EVENTS_GLOB, collapse = "|")
  files <- list.files(root, pattern = pat, recursive = TRUE,
                      full.names = TRUE, all.files = TRUE)
  needle <- paste0("/scale/", module_id, "/")
  files <- files[grepl(needle, files, fixed = TRUE)]
  unique(normalizePath(files))
}

threads_from_path <- function(path) {
  dirs <- strsplit(path, .Platform$file.sep)[[1]]
  m <- regmatches(dirs, regexpr("^threads-([0-9]+)$", dirs))
  m <- m[nchar(m) > 0]
  if (length(m)) return(as.integer(sub("threads-", "", m[1])))
  pjson <- file.path(dirname(path), "parameters.json")
  if (file.exists(pjson)) {
    p <- tryCatch(jsonlite::fromJSON(pjson), error = function(e) NULL)
    if (!is.null(p$threads)) return(as.integer(p$threads))
  }
  NA_integer_
}

args <- commandArgs(trailingOnly = TRUE)
root <- if (length(args)) args[1] else "out"

all_events <- lapply(names(BACKENDS), function(mod) {
  label <- BACKENDS[[mod]]
  files <- find_event_files(root, mod)
  if (!length(files)) {
    message("No event files for ", mod, " under ", root)
    return(NULL)
  }
  cat(sprintf("%-40s %d file(s)\n", label, length(files)))
  rows <- lapply(files, function(f) {
    df <- obkit::read_phase_ranges(f)
    df$elapsed_s <- (df$xmax - df$xmin) / 1000
    df$threads   <- threads_from_path(f)
    df
  })
  df <- do.call(rbind, rows)
  df$backend <- label
  df
})
all_events <- do.call(rbind, Filter(Negate(is.null), all_events))

if (is.null(all_events) || !nrow(all_events)) stop("No events found under ", root)

# One row per (backend, threads, stage)
wide <- all_events %>%
  select(backend, threads, event, elapsed_s) %>%
  pivot_wider(names_from = event, values_from = elapsed_s) %>%
  arrange(backend, threads)

cat("\nRaw timings (s):\n")
print(as.data.frame(wide))

# Speedup relative to threads=1 within each backend
speedup <- wide %>%
  group_by(backend) %>%
  mutate(across(where(is.numeric) & !threads,
                ~ first(.x[threads == 1]) / .x,
                .names = "{.col}_speedup")) %>%
  ungroup()

# Long form for plotting — only stages present in both backends
stage_cols <- intersect(
  c("convert", "load", "compute", "writing"),
  names(wide)
)

long_time <- wide %>%
  pivot_longer(all_of(stage_cols), names_to = "stage", values_to = "elapsed_s") %>%
  filter(!is.na(elapsed_s)) %>%
  mutate(stage = factor(stage, levels = c("convert", "load", "compute", "writing")))

speedup_cols <- paste0(stage_cols, "_speedup")
long_speedup <- speedup %>%
  select(backend, threads, all_of(speedup_cols)) %>%
  pivot_longer(all_of(speedup_cols), names_to = "stage", values_to = "speedup") %>%
  filter(!is.na(speedup)) %>%
  mutate(stage = sub("_speedup", "", stage),
         stage = factor(stage, levels = c("convert", "load", "compute", "writing")))

stage_colors <- c(convert = "#DD8452", load = "#4C72B0",
                  compute = "#55A868", writing = "#C44E52")
backend_ltys  <- c("HDF5 (h5ad)" = "solid", "zarr (blosc)" = "dashed")

dir.create("plots", showWarnings = FALSE)

p_time <- ggplot(long_time,
                 aes(x = threads, y = elapsed_s, color = stage,
                     linetype = backend, group = interaction(stage, backend))) +
  geom_line(linewidth = 0.8) +
  geom_point(size = 2.5) +
  scale_x_continuous(breaks = unique(long_time$threads)) +
  scale_color_manual(values = stage_colors) +
  scale_linetype_manual(values = backend_ltys) +
  labs(title = "scanpy thread scaling: wall time per stage",
       subtitle = "solid = HDF5 (h5ad)  |  dashed = zarr (blosc)",
       x = "threads", y = "elapsed (s)",
       color = "stage", linetype = "backend") +
  theme_minimal(base_size = 11)

p_speedup <- ggplot(long_speedup,
                    aes(x = threads, y = speedup, color = stage,
                        linetype = backend, group = interaction(stage, backend))) +
  geom_line(linewidth = 0.8) +
  geom_point(size = 2.5) +
  geom_abline(slope = 1, intercept = 0, linetype = "dotted", color = "grey50") +
  scale_x_continuous(breaks = unique(long_speedup$threads)) +
  scale_color_manual(values = stage_colors) +
  scale_linetype_manual(values = backend_ltys) +
  labs(title = "scanpy thread scaling: speedup vs 1 thread",
       subtitle = "solid = HDF5 (h5ad)  |  dashed = zarr (blosc)  |  dotted = ideal",
       x = "threads", y = "speedup (T1 / Tn)",
       color = "stage", linetype = "backend") +
  theme_minimal(base_size = 11)

ggsave("plots/thread_scaling_time.pdf",    p_time,    width = 7, height = 4)
ggsave("plots/thread_scaling_speedup.pdf", p_speedup, width = 7, height = 4)
cat("\nWrote plots/thread_scaling_time.pdf and plots/thread_scaling_speedup.pdf\n")
