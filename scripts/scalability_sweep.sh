#!/usr/bin/env bash
# Sweep the scanpy_scale concurrency level. For each N in the list, rewrite
# the `_replicates:` anchor in scalability.yaml to have N entries, run
# obrun, and stash the output dir as `runs/N<N>/`. Original yaml is
# restored on exit (or interrupt) via trap.
#
# The snakemake worker pool size (`obrun --cores`) is held constant and
# decoupled from N so we measure read contention at fixed worker
# availability. Set it with -c, or fall back to SNAKEMAKE_CORES env var,
# default 64.
#
# Usage:
#   scripts/scalability_sweep.sh                          # N=1..64, cores=64
#   scripts/scalability_sweep.sh -c 16 1 2 4 8            # cores=16, custom Ns
#   SNAKEMAKE_CORES=16 scripts/scalability_sweep.sh       # via env var

set -euo pipefail

cd "$(dirname "$0")/.."
ROOT="$(pwd)"

YAML="$ROOT/scalability.yaml"
BACKUP="$YAML.sweep.bak"
RUNS_DIR="$ROOT/runs"

CORES="${SNAKEMAKE_CORES:-64}"

while [[ $# -gt 0 ]]; do
    case "$1" in
        -c|--cores) CORES="$2"; shift 2 ;;
        --) shift; break ;;
        -*) echo "unknown flag: $1" >&2; exit 2 ;;
        *) break ;;
    esac
done

# Restore on any exit (including Ctrl-C). Idempotent: removes BACKUP at end.
restore() {
    if [[ -f "$BACKUP" ]]; then
        mv "$BACKUP" "$YAML"
    fi
}
trap restore EXIT INT TERM

cp "$YAML" "$BACKUP"

NS=("$@")
[[ ${#NS[@]} -eq 0 ]] && NS=(1 2 4 8 16 32 64)

mkdir -p "$RUNS_DIR"

echo "Sweep config: cores=$CORES, Ns=${NS[*]}"

for N in "${NS[@]}"; do
    echo "=========================================="
    echo "  N=$N readers, --cores $CORES (pool size)"
    echo "=========================================="

    # Build the replicate list, then splice it into the yaml in place of
    # the existing `_replicates: &replicates` anchor block. The block ends
    # at the first blank line.
    {
        echo "_replicates: &replicates"
        for i in $(seq 1 "$N"); do
            echo "  - replicate: $i"
        done
    } > /tmp/scalability_anchor.$$

    awk -v anchor_file="/tmp/scalability_anchor.$$" '
        BEGIN { in_block = 0; injected = 0 }
        /^_replicates: &replicates/ {
            in_block = 1
            if (!injected) {
                while ((getline line < anchor_file) > 0) print line
                close(anchor_file)
                injected = 1
            }
            next
        }
        in_block && /^[[:space:]]*$/ { in_block = 0; print; next }
        in_block { next }
        { print }
    ' "$BACKUP" > "$YAML"

    rm -f /tmp/scalability_anchor.$$

    # Wipe scale-stage outputs so snakemake re-runs them. fetch + prep are
    # cached across Ns — same dataset, same h5ad.
    find "$ROOT/out" -type d -name 'scale' 2>/dev/null \
        | while read -r d; do rm -rf "$d"; done

    obrun --dirty --unpinned --cores "$CORES"

    # Stash this run's outputs.
    DEST="$RUNS_DIR/N$N"
    rm -rf "$DEST"
    cp -r "$ROOT/out" "$DEST"
    echo "  → stashed in $DEST"
done

echo
echo "Sweep complete. Plot with:"
RUNDIRS=""
for N in "${NS[@]}"; do RUNDIRS="$RUNDIRS runs/N$N"; done
echo "  Rscript analysis/scalability.R$RUNDIRS"
