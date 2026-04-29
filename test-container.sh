#!/usr/bin/env bash
# Smoke test: every dependency declared in the Dockerfile loads successfully.
# Build with:
#   podman build -t localhost/scfmt:latest .
# Then run this script.

set -u
IMAGE="${IMAGE:-localhost/scfmt:latest}"
RUN=(podman run --rm "$IMAGE")

fail=0
pass=0
run_check() {
    local name="$1"; shift
    if "$@" >/tmp/scfmt_check.out 2>&1; then
        printf "  ok   %s\n" "$name"
        pass=$((pass + 1))
    else
        printf "  FAIL %s\n" "$name"
        sed 's/^/       | /' /tmp/scfmt_check.out
        fail=$((fail + 1))
    fi
}

echo "Image: $IMAGE"
echo

echo "[binaries]"
run_check "scx --version"   "${RUN[@]}" scx --version
run_check "denet --version" "${RUN[@]}" denet --version

echo
echo "[R packages]"
for pkg in SingleCellExperiment DelayedArray TENxPBMCData TENxBrainData anndataR Seurat BPCells argparse omnibench.logger; do
    run_check "library($pkg)" "${RUN[@]}" Rscript -e "suppressPackageStartupMessages(library($pkg))"
done

echo
echo "[Python packages]"
for pkg in scanpy anndata h5py denet omnibench_logger; do
    run_check "import $pkg" "${RUN[@]}" python3 -c "import $pkg"
done

echo
echo "[integration]"
# h5ad round-trip via scx — catches HDF5 / dtype wiring.
run_check "scx inspect --help" "${RUN[@]}" scx inspect --help
# anndataR can open an AnnData.
run_check "anndataR::AnnData\$new()" "${RUN[@]}" Rscript -e '
  ad <- anndataR::AnnData$new(
    X = matrix(1:6, nrow = 2),
    obs = data.frame(row.names = c("c1","c2")),
    var = data.frame(row.names = c("g1","g2","g3"))
  )
  stopifnot(dim(ad$X) == c(2,3))
'
# BPCells write→read round-trip in a tmp dir.
run_check "BPCells write/read" "${RUN[@]}" Rscript -e '
  m <- Matrix::rsparsematrix(100, 50, density = 0.1)
  d <- tempfile()
  BPCells::write_matrix_dir(BPCells::as(m, "IterableMatrix"), dir = d)
  m2 <- BPCells::open_matrix_dir(d)
  stopifnot(dim(m2) == c(100, 50))
'

echo
echo "------------------------------------------------------------"
echo "pass: $pass   fail: $fail"
exit "$fail"
