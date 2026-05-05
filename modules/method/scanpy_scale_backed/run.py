#!/usr/bin/env python3
"""scanpy scalability replicate (backed/memory-mapped mode).

Same pipeline as scanpy_scale but opens the h5ad with backed='r' so the
kernel memory-maps the file instead of bulk-loading it.  The actual I/O
pressure therefore moves from the `load` stage into `compute` (where
to_memory() forces the read), giving a different contention signature.

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
    p.add_argument("--name", default="scanpy_scale_backed")
    p.add_argument("--n_hvg", type=int, default=2000)
    p.add_argument("--n_comp", type=int, default=50)
    args = p.parse_args()

    out = Path(args.output_dir)
    out.mkdir(parents=True, exist_ok=True)
    init_logger(str(out))

    rep = args.replicate
    dataset = Path(args.h5ad).stem

    # backed='r': mmap the file — open is near-instant, no bulk I/O here
    with stage("load", rep):
        adata = sc.read_h5ad(args.h5ad, backed="r")

    with stage("compute", rep):
        # to_memory() is where actual file read (and contention) happens
        adata = adata.to_memory()
        sc.pp.normalize_total(adata, target_sum=1e4)
        sc.pp.log1p(adata)
        sc.pp.highly_variable_genes(adata, n_top_genes=args.n_hvg)
        adata = adata[:, adata.var.highly_variable].copy()
        sc.pp.scale(adata, max_value=10)
        sc.tl.pca(adata, n_comps=args.n_comp, random_state=rep * 1000)
        sc.pp.neighbors(adata, n_neighbors=15, random_state=rep * 1000)

    out_h5ad = out / f"{dataset}.h5ad"
    with stage("writing", rep):
        adata.write_h5ad(out_h5ad)

    print(f"rep={rep} done (backed): {out_h5ad}")


if __name__ == "__main__":
    main()
