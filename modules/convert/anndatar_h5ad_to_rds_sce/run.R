#!/usr/bin/env Rscript
# Convert h5ad → rds(SingleCellExperiment) using anndataR.
# CLI: --h5ad <input> --output_dir <dir> [--name <str>]

suppressPackageStartupMessages(library(argparse))

parser <- ArgumentParser()
parser$add_argument("--h5ad",       required = TRUE)
parser$add_argument("--output_dir", required = TRUE)
parser$add_argument("--name",       default  = NULL)  # ignored; derived
args <- parser$parse_args()
args$name <- sub("\\.h5ad$", "", basename(args$h5ad))
dir.create(args$output_dir, showWarnings = FALSE, recursive = TRUE)

suppressPackageStartupMessages({
  library(omnibench.logger)
  library(anndataR)
  library(SingleCellExperiment)
})
init_logger(args$output_dir)
source(file.path(Sys.getenv("SCFMT_ROOT", "/workspace"), "scripts", "stage.R"))

out <- file.path(args$output_dir, paste0(args$name, ".sce.rds"))

stage("load", {
  ad  <- anndataR::read_h5ad(args$h5ad)
  sce <- anndataR:::as_SingleCellExperiment(ad)  # not yet exported
})

stage("write", {
  saveRDS(sce, out)
})

message("Wrote: ", out)
