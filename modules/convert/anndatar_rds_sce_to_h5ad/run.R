#!/usr/bin/env Rscript
# Convert rds(SingleCellExperiment) → h5ad using anndataR.
# CLI: --rds_sce <input> --output_dir <dir> [--name <str>]

suppressPackageStartupMessages(library(argparse))

parser <- ArgumentParser()
parser$add_argument("--rds_sce",    required = TRUE)
parser$add_argument("--output_dir", required = TRUE)
parser$add_argument("--name",       default  = NULL)  # ignored; derived
args <- parser$parse_args()
args$name <- paste0(sub("\\.sce\\.rds$", "", basename(args$rds_sce)), ".rt_anndatar")
dir.create(args$output_dir, showWarnings = FALSE, recursive = TRUE)

suppressPackageStartupMessages({
  library(omnibench.logger)
  library(anndataR)
  library(SingleCellExperiment)
})
init_logger(args$output_dir)
source(file.path(Sys.getenv("SCFMT_ROOT", "/workspace"), "scripts", "stage.R"))

out <- file.path(args$output_dir, paste0(args$name, ".h5ad"))

stage("load", {
  sce <- readRDS(args$rds_sce)
})

stage("write", {
  ad <- anndataR::as_AnnData(sce, x_mapping = "counts")
  anndataR::write_h5ad(ad, out)
})

message("Wrote: ", out)
