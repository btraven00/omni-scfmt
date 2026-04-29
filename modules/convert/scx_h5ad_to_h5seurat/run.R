#!/usr/bin/env Rscript
# Convert h5ad → h5seurat using scx.
# CLI: --h5ad <input> --output_dir <dir> [--name <str>]

suppressPackageStartupMessages(library(argparse))

parser <- ArgumentParser()
parser$add_argument("--h5ad",       required = TRUE)
parser$add_argument("--output_dir", required = TRUE)
parser$add_argument("--name",       default  = NULL)  # ignored; derived
args <- parser$parse_args()
# The omnibenchmark backend passes --name <module_id>, but the dataset name
# must come from the input filename for snakemake's output pattern to match.
args$name <- sub("\\.h5ad$", "", basename(args$h5ad))
dir.create(args$output_dir, showWarnings = FALSE, recursive = TRUE)

suppressPackageStartupMessages(library(omnibench.logger))
init_logger(args$output_dir)
source(file.path(Sys.getenv("SCFMT_ROOT", "/workspace"), "scripts", "stage.R"))

out <- file.path(args$output_dir, paste0(args$name, ".h5seurat"))

stage("convert", {
  # scx streams read→write internally, so we can't split load/write here.
  # Inner denet timeline will still show the single convert block.
  rc <- system2("scx", c("convert", args$h5ad, out))
  if (rc != 0) stop("scx convert failed: exit ", rc)
})

message("Wrote: ", out)
