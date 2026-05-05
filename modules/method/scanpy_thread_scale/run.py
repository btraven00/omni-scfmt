#!/usr/bin/env python3
"""scanpy thread-scaling benchmark.

Measures load/compute/writing wall time at a given thread count. BLAS and
OpenMP thread counts are set by the `scale` wrapper before numpy is imported.
sklearn n_jobs is set explicitly here for operations that accept it.

CLI: --h5ad INPUT --output_dir DIR --threads N [--name STR]
     [--n_hvg N] [--n_comp N]
"""
from __future__ import annotations

import argparse
import time
from contextlib import contextmanager
from pathlib import Path

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
    p.add_argument("--name", default="scanpy_thread_scale")
    p.add_argument("--n_hvg", type=int, default=2000)
    p.add_argument("--n_comp", type=int, default=50)
    args = p.parse_args()

    out = Path(args.output_dir)
    out.mkdir(parents=True, exist_ok=True)
    init_logger(str(out))

    n = args.threads
    dataset = Path(args.h5ad).stem

    with stage("load", n):
        adata = sc.read_h5ad(args.h5ad)

    with stage("compute", n):
        sc.pp.normalize_total(adata, target_sum=1e4)
        sc.pp.log1p(adata)
        sc.pp.highly_variable_genes(adata, n_top_genes=args.n_hvg)
        adata = adata[:, adata.var.highly_variable].copy()
        sc.pp.scale(adata, max_value=10)
        sc.tl.pca(adata, n_comps=args.n_comp)
        sc.pp.neighbors(adata, n_neighbors=15, n_jobs=n)

    out_h5ad = out / f"{dataset}.h5ad"
    with stage("writing", n):
        adata.write_h5ad(out_h5ad)

    print(f"threads={n} done: {out_h5ad}")


if __name__ == "__main__":
    main()
