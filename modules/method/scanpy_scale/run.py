#!/usr/bin/env python3
"""scanpy scalability replicate: load → compute → writing under obkit.logger.

One process per replicate; omnibenchmark fans out N replicates as parameter
expansions, so N copies of this script run concurrently against the same
input h5ad.

CLI: --h5ad INPUT --output_dir DIR --replicate N [--name STR]
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
def stage(name: str, replicate: int):
    attrs = {"replicate": replicate}
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
    p.add_argument("--replicate", type=int, required=True)
    p.add_argument("--name", default="scanpy_scale")
    p.add_argument("--n_hvg", type=int, default=2000)
    p.add_argument("--n_comp", type=int, default=50)
    args = p.parse_args()

    out = Path(args.output_dir)
    out.mkdir(parents=True, exist_ok=True)
    init_logger(str(out))

    rep = args.replicate
    dataset = Path(args.h5ad).stem

    with stage("load", rep):
        adata = sc.read_h5ad(args.h5ad)

    with stage("compute", rep):
        sc.pp.normalize_total(adata, target_sum=1e4)
        sc.pp.log1p(adata)
        sc.pp.highly_variable_genes(adata, n_top_genes=args.n_hvg)
        adata = adata[:, adata.var.highly_variable].copy()
        sc.pp.scale(adata, max_value=10)
        # PCA is deterministic given input, but neighbors + UMAP take a
        # random_state — seed per replicate so we get genuine across-rep
        # variance independent of contention noise. neighbors is the heavy
        # step here (kNN graph build dominates wall time on 68k+).
        sc.tl.pca(adata, n_comps=args.n_comp, random_state=rep * 1000)
        sc.pp.neighbors(adata, n_neighbors=15, random_state=rep * 1000)
        sc.tl.umap(adata, random_state=rep * 1000)

    out_h5ad = out / f"{dataset}.h5ad"
    with stage("writing", rep):
        adata.write_h5ad(out_h5ad)

    print(f"rep={rep} done: {out_h5ad}")


if __name__ == "__main__":
    main()
