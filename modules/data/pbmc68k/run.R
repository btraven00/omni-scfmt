#!/usr/bin/env Rscript
# Load PBMC 68k from TENxPBMCData, materialize counts sparse, write h5ad.
# CLI: --output_dir <dir> [--name <str>]

suppressPackageStartupMessages(library(argparse))

parser <- ArgumentParser()
parser$add_argument("--output_dir", required = TRUE)
parser$add_argument("--name",       default  = "pbmc68k")
args <- parser$parse_args()

dir.create(args$output_dir, showWarnings = FALSE, recursive = TRUE)

suppressPackageStartupMessages({
  library(omnibench.logger)
  library(SingleCellExperiment)
  library(TENxPBMCData)
  library(Matrix)
  library(anndataR)
})
init_logger(args$output_dir)
source(file.path(Sys.getenv("SCFMT_ROOT", "/workspace"), "scripts", "stage.R"))

stage("fetch", {
  sce <- TENxPBMCData("pbmc68k")
})

stage("materialize", {
  if (!is.null(rowData(sce)$Symbol_TENx)) {
    syms <- as.character(rowData(sce)$Symbol_TENx)
    syms[is.na(syms)] <- rownames(sce)[is.na(syms)]
    rownames(sce) <- make.unique(syms)
  }
  assay(sce, "counts") <- as(counts(sce), "dgCMatrix")
})

out <- file.path(args$output_dir, paste0(args$name, ".h5ad"))
stage("write_h5ad", {
  ad <- anndataR::as_AnnData(sce, x_mapping = "counts")
  anndataR::write_h5ad(ad, out)
})

message("Saved: ", out, " (", ncol(sce), " cells, ", nrow(sce), " genes)")
