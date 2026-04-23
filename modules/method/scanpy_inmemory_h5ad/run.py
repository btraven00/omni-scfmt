#!/usr/bin/env python3
"""scanpy, in-memory, native h5ad.

LOAD → ACCESS → WRITE → normalize → HVG → PCA.
CLI: --h5ad INPUT --output_dir DIR [--name STR] [--n_hvg N] [--n_comp N]
     [--access_n N] [--seed S]
"""
from __future__ import annotations

import argparse
import os
import subprocess
import sys
from pathlib import Path

import numpy as np
import pandas as pd
import scanpy as sc
import anndata as ad_mod
from omnibench_logger import init_logger

sys.path.insert(0, str(Path(os.environ.get("SCFMT_ROOT", "/workspace")) / "scripts"))
from stage import stage, stage_timings_df  # noqa: E402


def main() -> None:
    p = argparse.ArgumentParser()
    p.add_argument("--h5ad", required=True)
    p.add_argument("--rds_seurat", default=None)  # provided by pipeline, not used
    p.add_argument("--output_dir", required=True)
    p.add_argument("--name", default="scanpy_inmemory_h5ad")
    p.add_argument("--n_hvg", type=int, default=2000)
    p.add_argument("--n_comp", type=int, default=50)
    p.add_argument("--access_n", type=int, default=10_000)
    p.add_argument("--seed", type=int, default=1)
    args = p.parse_args()
    dataset = Path(args.h5ad).stem

    out = Path(args.output_dir); out.mkdir(parents=True, exist_ok=True)
    init_logger(str(out))

    with stage("load"):
        adata = sc.read_h5ad(args.h5ad)

    with stage("access"):
        rng = np.random.default_rng(args.seed)
        n = min(args.access_n, adata.n_obs)
        idx = rng.choice(adata.n_obs, size=n, replace=False)
        access_sum = float(np.asarray(adata.X[idx, :].sum()))

    scratch_out = out / f"{args.name}.scratch.h5ad"
    with stage("write"):
        adata.write_h5ad(scratch_out)

    with stage("normalize"):
        sc.pp.normalize_total(adata, target_sum=1e4)
        sc.pp.log1p(adata)

    with stage("hvg"):
        sc.pp.highly_variable_genes(adata, n_top_genes=args.n_hvg)

    with stage("pca"):
        sc.pp.scale(adata, max_value=10)
        sc.tl.pca(adata, n_comps=args.n_comp)

    timing_df = stage_timings_df(dataset=dataset)
    timing_df.to_csv(out / f"{dataset}.timing.tsv", sep="\t", index=False)

    pd.DataFrame({"access_sum": [access_sum]}).to_csv(
        out / f"{dataset}.access.tsv", sep="\t", index=False
    )
    results_out = str(out / f"{dataset}.results.rds")
    subprocess.run(["Rscript", "-e", f'saveRDS(NULL, "{results_out}")'], check=True)
    print(f"Done: {dataset}")


if __name__ == "__main__":
    main()
