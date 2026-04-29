#!/usr/bin/env python3
"""scanpy, in-memory, foreign h5seurat (convert via scx → h5ad, then load).

LOAD includes the scx conversion step; that's the real user cost of
handing a scanpy pipeline a Seurat-native file.
CLI: --h5seurat INPUT --output_dir DIR [--name STR]
"""
from __future__ import annotations

import argparse
import os
import subprocess
import sys
import tempfile
from pathlib import Path

import numpy as np
import pandas as pd
import scanpy as sc
from omnibench_logger import init_logger

sys.path.insert(0, str(Path(os.environ.get("SCFMT_ROOT", "/workspace")) / "scripts"))
from stage import stage, stage_timings_df  # noqa: E402


def main() -> None:
    p = argparse.ArgumentParser()
    p.add_argument("--h5seurat", required=True)
    p.add_argument("--output_dir", required=True)
    p.add_argument("--name", default="scanpy_inmemory_h5seurat")
    p.add_argument("--n_hvg", type=int, default=2000)
    p.add_argument("--n_comp", type=int, default=50)
    p.add_argument("--access_n", type=int, default=10_000)
    p.add_argument("--seed", type=int, default=1)
    args = p.parse_args()

    out = Path(args.output_dir); out.mkdir(parents=True, exist_ok=True)
    init_logger(str(out))

    tmp_h5ad = Path(tempfile.mkdtemp(prefix="scfmt_")) / "converted.h5ad"
    with stage("load"):
        subprocess.run(
            ["scx", "convert", args.h5seurat, str(tmp_h5ad)], check=True
        )
        adata = sc.read_h5ad(tmp_h5ad)

    with stage("access"):
        rng = np.random.default_rng(args.seed)
        n = min(args.access_n, adata.n_obs)
        idx = rng.choice(adata.n_obs, size=n, replace=False)
        access_sum = float(np.asarray(adata.X[idx, :].sum()))

    # WRITE back to h5seurat via scx in reverse (fair round-trip cost).
    scratch_out = out / f"{args.name}.scratch.h5seurat"
    with stage("write"):
        tmp_h5ad_out = out / f"{args.name}.scratch.h5ad"
        adata.write_h5ad(tmp_h5ad_out)
        subprocess.run(
            ["scx", "convert", str(tmp_h5ad_out), str(scratch_out)], check=True
        )

    with stage("normalize"):
        sc.pp.normalize_total(adata, target_sum=1e4)
        sc.pp.log1p(adata)
    with stage("hvg"):
        sc.pp.highly_variable_genes(adata, n_top_genes=args.n_hvg)
    with stage("pca"):
        sc.pp.scale(adata, max_value=10)
        sc.tl.pca(adata, n_comps=args.n_comp)

    timing_df = stage_timings_df(dataset=args.name)
    timing_df.to_csv(out / f"{args.name}.timing.tsv", sep="\t", index=False)
    pd.DataFrame({"access_sum": [access_sum]}).to_csv(
        out / f"{args.name}.access.tsv", sep="\t", index=False
    )
    print(f"Done: {args.name}")


if __name__ == "__main__":
    main()
