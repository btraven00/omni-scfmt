#!/usr/bin/env bash
# Sweep the scanpy_scale concurrency level. For each N in the list, rewrite
# the `_replicates:` anchor in scalability.yaml to have N entries, run
# obrun, and stash the output dir as `runs/N<N>/`. Original yaml is
# restored on exit (or interrupt) via trap.
#
# Usage:
#   scripts/scalability_sweep.sh                  # default 1 2 3 5 8
#   scripts/scalability_sweep.sh 1 2 4 8 16       # custom

set -euo pipefail

cd "$(dirname "$0")/.."
ROOT="$(pwd)"

YAML="$ROOT/scalability.yaml"
BACKUP="$YAML.sweep.bak"
RUNS_DIR="$ROOT/runs"

# Restore on any exit (including Ctrl-C). Idempotent: removes BACKUP at end.
restore() {
    if [[ -f "$BACKUP" ]]; then
        mv "$BACKUP" "$YAML"
    fi
}
trap restore EXIT INT TERM

cp "$YAML" "$BACKUP"

NS=("$@")
[[ ${#NS[@]} -eq 0 ]] && NS=(1 2 3 5 8)

mkdir -p "$RUNS_DIR"

for N in "${NS[@]}"; do
    echo "=========================================="
    echo "  Concurrency sweep: N=$N"
    echo "=========================================="

    # Build the replicate list as a here-doc, then splice it in by deleting
    # the existing anchor block (`_replicates: &replicates` through the next
    # blank line) and inserting the new one in its place.
    {
        echo "_replicates: &replicates"
        for i in $(seq 1 "$N"); do
            echo "  - replicate: $i"
        done
    } > /tmp/scalability_anchor.$$

    # Replace the anchor block in the yaml. The block ends at the first
    # blank line after the anchor declaration.
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

    # Wipe scale-stage outputs from prior N so snakemake re-runs them.
    rm -rf "$ROOT/out/fetch/pbmc3k/.e024bea6/prep/pbmc3k_h5ad/.default/scale" || true

    obrun --dirty --unpinned --cores "$N"

    # Stash this run's outputs.
    DEST="$RUNS_DIR/N$N"
    rm -rf "$DEST"
    cp -r "$ROOT/out" "$DEST"
    echo "  → stashed in $DEST"
done

echo
echo "Sweep complete. Plot with:"
echo "  Rscript analysis/scalability.R $(printf '%s ' "${NS[@]}" | sed 's@\([0-9]*\) @runs/N\1 @g')"
