#!/usr/bin/env bash
# Build and push container images for omni-scfmt.
# Usage: containers/build.sh [--push] [omni-data|scanpy|all]
set -euo pipefail

cd "$(dirname "$0")"

REGISTRY="ghcr.io/btraven00"
PUSH=0

while [[ $# -gt 0 ]]; do
    case "$1" in
        --push) PUSH=1; shift ;;
        --) shift; break ;;
        *) break ;;
    esac
done

TARGETS=("$@")
[[ ${#TARGETS[@]} -eq 0 ]] && TARGETS=(omni-data scanpy)

declare -A IMAGE_NAMES=(
    [omni-data]="$REGISTRY/omni-data:latest"
    [scanpy]="$REGISTRY/scfmt-scanpy:latest"
)

for target in "${TARGETS[@]}"; do
    [[ $target == "all" ]] && { TARGETS=(omni-data scanpy); continue; }
    image="${IMAGE_NAMES[$target]}"
    echo "==> Building $image from $target/"
    podman build -t "$image" -f "$target/Containerfile" "$target/"
    if [[ $PUSH -eq 1 ]]; then
        echo "==> Pushing $image"
        podman push "$image"
    fi
done

echo "Done. To push: containers/build.sh --push"
