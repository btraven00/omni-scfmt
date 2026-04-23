#!/usr/bin/env Rscript
# Seurat, in-memory, native .rds (serialized Seurat object).
# This is the true Seurat-native on-disk path without SeuratDisk: one
# readRDS() call, no conversion.
# LOAD → ACCESS → WRITE → normalize → HVG → PCA.
# CLI: --rds_seurat <input> --output_dir <dir> [--name <str>]

suppressPackageStartupMessages(library(argparse))
parser <- ArgumentParser()
parser$add_argument("--rds_seurat", required = TRUE)
parser$add_argument("--h5ad",       default  = NULL)  # provided by pipeline, not used
parser$add_argument("--output_dir", required = TRUE)
parser$add_argument("--name",       default  = "seurat_inmemory_rds")
parser$add_argument("--n_hvg",      type = "integer", default = 2000)
parser$add_argument("--n_comp",     type = "integer", default = 50)
parser$add_argument("--access_n",   type = "integer", default = 10000)
parser$add_argument("--seed",       type = "integer", default = 1)
args <- parser$parse_args()
dataset <- sub("\\.seurat\\.rds$", "", basename(args$rds_seurat))
dir.create(args$output_dir, showWarnings = FALSE, recursive = TRUE)

suppressPackageStartupMessages({
  library(omnibench.logger)
  library(Seurat); library(SeuratObject); library(Matrix)
})
init_logger(args$output_dir)
source(file.path(Sys.getenv("SCFMT_ROOT", "/workspace"), "scripts", "stage.R"))

stage("load", {
  obj <- readRDS(args$rds_seurat)
  if ("X" %in% Layers(obj) && !"counts" %in% Layers(obj))
    obj <- CreateSeuratObject(counts = LayerData(obj, layer = "X"), meta.data = obj@meta.data)
})

stage("access", {
  set.seed(args$seed)
  idx <- sample.int(ncol(obj), min(args$access_n, ncol(obj)))
  counts_mat <- LayerData(obj, layer = "counts")
  access_sum <- sum(Matrix::colSums(counts_mat[, idx, drop = FALSE]))
  rm(counts_mat); invisible(gc(full = TRUE, verbose = FALSE))
})

scratch_out <- file.path(args$output_dir, paste0(dataset, ".scratch.rds"))
stage("write", {
  saveRDS(obj, scratch_out)
})

stage("normalize", { obj <- NormalizeData(obj, verbose = FALSE) })
stage("hvg",       { obj <- FindVariableFeatures(obj, nfeatures = args$n_hvg, verbose = FALSE) })
stage("pca", {
  obj <- ScaleData(obj, features = VariableFeatures(obj), verbose = FALSE)
  obj <- RunPCA(obj, npcs = args$n_comp, verbose = FALSE)
})

timing_df <- stage_timings_df(dataset = dataset)
write.table(timing_df,
  file = file.path(args$output_dir, paste0(dataset, ".timing.tsv")),
  sep = "\t", quote = FALSE, row.names = FALSE)
saveRDS(list(obj = obj, timing = timing_df, access_sum = access_sum),
  file.path(args$output_dir, paste0(dataset, ".results.rds")))
message("Done: ", dataset)
