#!/usr/bin/env Rscript
# Load 10x 1.3M mouse brain from TENxBrainData, deterministic cell subsample,
# materialize counts sparse, write h5ad. Shared across all brain1m_Nk
# modules — size is driven by --n_cells.
# CLI: --output_dir <dir> --n_cells <int> [--seed <int>] [--name <str>]

suppressPackageStartupMessages(library(argparse))

parser <- ArgumentParser()
parser$add_argument("--output_dir", required = TRUE)
parser$add_argument("--n_cells",    required = TRUE, type = "integer")
parser$add_argument("--seed",       type = "integer", default = 1)
parser$add_argument("--name",       default = NULL)
args <- parser$parse_args()
if (is.null(args$name)) args$name <- sprintf("brain1m_%dk", args$n_cells %/% 1000)

dir.create(args$output_dir, showWarnings = FALSE, recursive = TRUE)

suppressPackageStartupMessages({
  library(omnibench.logger)
  library(SingleCellExperiment)
  library(TENxBrainData)
  library(DelayedArray)
  library(Matrix)
  library(anndataR)
})
init_logger(args$output_dir)
source(file.path(Sys.getenv("SCFMT_ROOT", "/workspace"), "scripts", "stage.R"))

stage("fetch", {
  sce <- TENxBrainData()   # HDF5-backed
})

stage("subsample", {
  set.seed(args$seed)
  n_total <- ncol(sce)
  n_keep  <- min(args$n_cells, n_total)
  idx     <- sort(sample.int(n_total, n_keep))
  sce     <- sce[, idx]
})

stage("materialize", {
  # HDF5-backed slice → in-memory sparse; the memory peak stage.
  assay(sce, "counts") <- as(as(counts(sce), "CsparseMatrix"), "dgCMatrix")
  if (!is.null(rowData(sce)$Symbol)) {
    syms <- as.character(rowData(sce)$Symbol)
    syms[is.na(syms)] <- rownames(sce)[is.na(syms)]
    rownames(sce) <- make.unique(syms)
  }
})

out <- file.path(args$output_dir, paste0(args$name, ".h5ad"))
stage("write_h5ad", {
  ad <- anndataR::as_AnnData(sce, x_mapping = "counts")
  anndataR::write_h5ad(ad, out)
})

message("Saved: ", out, " (", ncol(sce), " cells, ", nrow(sce), " genes)")
