#!/usr/bin/env Rscript
# Convert h5seurat → h5ad using scx.
# CLI: --h5seurat <input> --output_dir <dir> [--name <str>]

suppressPackageStartupMessages(library(argparse))

parser <- ArgumentParser()
parser$add_argument("--h5seurat",   required = TRUE)
parser$add_argument("--output_dir", required = TRUE)
parser$add_argument("--name",       default  = NULL)  # ignored; derived
args <- parser$parse_args()
# Derive dataset name from input so the output matches snakemake's pattern
# (the backend passes --name <module_id>).
args$name <- paste0(sub("\\.h5seurat$", "", basename(args$h5seurat)), ".rt_scx")
dir.create(args$output_dir, showWarnings = FALSE, recursive = TRUE)

suppressPackageStartupMessages(library(omnibench.logger))
init_logger(args$output_dir)
source(file.path(Sys.getenv("SCFMT_ROOT", "/workspace"), "scripts", "stage.R"))

out <- file.path(args$output_dir, paste0(args$name, ".h5ad"))

stage("convert", {
  rc <- system2("scx", c("convert", args$h5seurat, out))
  if (rc != 0) stop("scx convert failed: exit ", rc)
})

message("Wrote: ", out)
