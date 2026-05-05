#!/usr/bin/env python3
"""scanpy thread-scaling benchmark — zarr storage backend.

Same compute pipeline as scanpy_thread_scale but reads via zarr, whose blosc
decompressor respects numcodecs.blosc.set_nthreads(n).  Four obkit stages:

  convert  — h5ad → zarr on disk (parallel blosc compression)
  load     — zarr → AnnData (parallel blosc decompression)
  compute  — normalize → log1p → HVG → scale → PCA → neighbors
  writing  — AnnData → h5ad output

CLI: --h5ad INPUT --output_dir DIR --threads N [--name STR]
     [--n_hvg N] [--n_comp N]
"""
from __future__ import annotations

import argparse
import shutil
import time
from contextlib import contextmanager
from pathlib import Path

import anndata
import numcodecs
import scanpy as sc
from obkit.logger import init_logger, emit


@contextmanager
def stage(name: str, threads: int):
    attrs = {"threads": threads}
    emit(name, "start", attrs=attrs)
    t0 = time.monotonic()
    try:
        yield
    finally:
        attrs["elapsed_s"] = time.monotonic() - t0
        emit(name, "end", attrs=attrs)


def main() -> None:
    p = argparse.ArgumentParser()
    p.add_argument("--h5ad", required=True)
    p.add_argument("--output_dir", required=True)
    p.add_argument("--threads", type=int, required=True)
    p.add_argument("--name", default="scanpy_thread_scale_zarr")
    p.add_argument("--n_hvg", type=int, default=2000)
    p.add_argument("--n_comp", type=int, default=50)
    args = p.parse_args()

    out = Path(args.output_dir)
    out.mkdir(parents=True, exist_ok=True)
    init_logger(str(out))

    n = args.threads
    numcodecs.blosc.set_nthreads(n)

    dataset = Path(args.h5ad).stem
    zarr_path = out / f"{dataset}.zarr"

    with stage("convert", n):
        adata_tmp = sc.read_h5ad(args.h5ad)
        adata_tmp.write_zarr(zarr_path)
        del adata_tmp

    with stage("load", n):
        adata = anndata.read_zarr(zarr_path)

    with stage("compute", n):
        sc.pp.normalize_total(adata, target_sum=1e4)
        sc.pp.log1p(adata)
        sc.pp.highly_variable_genes(adata, n_top_genes=args.n_hvg)
        adata = adata[:, adata.var.highly_variable].copy()
        sc.pp.scale(adata, max_value=10)
        sc.tl.pca(adata, n_comps=args.n_comp)
        sc.pp.neighbors(adata, n_neighbors=15)

    out_h5ad = out / f"{dataset}.h5ad"
    with stage("writing", n):
        adata.write_h5ad(out_h5ad)

    shutil.rmtree(zarr_path, ignore_errors=True)
    print(f"threads={n} done (zarr): {out_h5ad}")


if __name__ == "__main__":
    main()
