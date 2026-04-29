#!/usr/bin/env Rscript
# Seurat, BPCells-backed (counts live on disk; obs/var in memory).
# LOAD includes one-time BPCells write if the bpcells dir doesn't exist.
# CLI: --h5ad <input> --output_dir <dir> --bpcells_dir <path> [--name <str>]

suppressPackageStartupMessages(library(argparse))
parser <- ArgumentParser()
parser$add_argument("--h5ad",        required = TRUE)
parser$add_argument("--output_dir",  required = TRUE)
parser$add_argument("--bpcells_dir", required = TRUE)
parser$add_argument("--name",        default  = "seurat_backed_bpcells")
parser$add_argument("--n_hvg",       type = "integer", default = 2000)
parser$add_argument("--n_comp",      type = "integer", default = 50)
parser$add_argument("--access_n",    type = "integer", default = 10000)
parser$add_argument("--seed",        type = "integer", default = 1)
args <- parser$parse_args()
dir.create(args$output_dir, showWarnings = FALSE, recursive = TRUE)

suppressPackageStartupMessages({
  library(omnibench.logger)
  library(Seurat); library(SeuratObject); library(Matrix)
  library(BPCells); library(anndataR)
})
init_logger(args$output_dir)
source(file.path(Sys.getenv("SCFMT_ROOT", "/workspace"), "scripts", "stage.R"))

# ----- LOAD: one-time BPCells write, then mount as backed matrix ------------
stage("load", {
  if (!dir.exists(args$bpcells_dir)) {
    ad <- anndataR::read_h5ad(args$h5ad)
    m  <- ad$X   # cells × genes in AnnData; transpose to genes × cells
    m  <- Matrix::t(as(m, "CsparseMatrix"))
    BPCells::write_matrix_dir(BPCells::as(m, "IterableMatrix"),
                               dir = args$bpcells_dir)
    rm(ad, m); invisible(gc(full = TRUE, verbose = FALSE))
  }
  mat <- BPCells::open_matrix_dir(args$bpcells_dir)   # lazy / on-disk
  obj <- CreateSeuratObject(counts = mat)
})

stage("access", {
  set.seed(args$seed)
  idx <- sample.int(ncol(obj), min(args$access_n, ncol(obj)))
  access_sum <- sum(Matrix::colSums(GetAssayData(obj, slot = "counts")[, idx, drop = FALSE]))
})

# WRITE as a fresh bpcells dir (typical user action when subsetting).
scratch_out <- file.path(args$output_dir, paste0(args$name, ".scratch.bpcells"))
stage("write", {
  BPCells::write_matrix_dir(GetAssayData(obj, slot = "counts"),
                             dir = scratch_out, overwrite = TRUE)
})

stage("normalize",   { obj <- NormalizeData(obj, verbose = FALSE) })
stage("hvg",         { obj <- FindVariableFeatures(obj, nfeatures = args$n_hvg, verbose = FALSE) })
stage("pca", {
  obj <- ScaleData(obj, features = VariableFeatures(obj), verbose = FALSE)
  obj <- RunPCA(obj, npcs = args$n_comp, verbose = FALSE)
})

timing_df <- stage_timings_df(dataset = args$name)
write.table(timing_df,
  file = file.path(args$output_dir, paste0(args$name, ".timing.tsv")),
  sep = "\t", quote = FALSE, row.names = FALSE)
saveRDS(list(timing = timing_df, access_sum = access_sum),
  file.path(args$output_dir, paste0(args$name, ".results.rds")))
message("Done: ", args$name)
