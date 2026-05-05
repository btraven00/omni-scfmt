#!/usr/bin/env python3
"""Convert omni-data's pbmc3k tarball into a cached h5ad.

Inputs:  --tarball <path-to-pbmc3k.tar.gz> (the 10x filtered matrix bundle)
Outputs: <output_dir>/<name>.h5ad
"""
from __future__ import annotations

import argparse
import tarfile
import tempfile
from pathlib import Path

import scanpy as sc


def main() -> None:
    p = argparse.ArgumentParser()
    p.add_argument("--tarball", required=True)
    p.add_argument("--output_dir", required=True)
    # ob passes --name <module_id>; ignore it and derive the basename from
    # the upstream dataset (encoded in the tarball filename) so the output
    # path matches the stage template "{dataset}.h5ad".
    p.add_argument("--name", default=None)
    args = p.parse_args()

    dataset = Path(args.tarball).name.split(".")[0]
    out = Path(args.output_dir)
    out.mkdir(parents=True, exist_ok=True)

    with tempfile.TemporaryDirectory() as tmp:
        with tarfile.open(args.tarball, "r:gz") as tar:
            tar.extractall(tmp)
        # 10x bundle layout: filtered_gene_bc_matrices/hg19/{matrix.mtx,genes.tsv,barcodes.tsv}
        candidates = list(Path(tmp).rglob("matrix.mtx"))
        if not candidates:
            raise SystemExit(f"matrix.mtx not found inside {args.tarball}")
        mtx_dir = candidates[0].parent
        adata = sc.read_10x_mtx(str(mtx_dir), var_names="gene_symbols", make_unique=True)

    h5ad = out / f"{dataset}.h5ad"
    adata.write_h5ad(h5ad)
    print(f"Wrote {h5ad}: {adata.n_obs} cells x {adata.n_vars} genes")


if __name__ == "__main__":
    main()
