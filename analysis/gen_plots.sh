#!/usr/bin/env bash
# Regenerate all omni-scfmt plots from out/.
# Run from the repo root (or set --out_dir / --plot_dir).
set -euo pipefail

HERE="$(cd "$(dirname "$0")" && pwd)"
OUT="${1:-out}"
PLOTS="${2:-plots}"

Rscript "$HERE/plot_penalty_curve.R"     --out_dir "$OUT" --plot_dir "$PLOTS"
Rscript "$HERE/plot_phase_breakdown.R"   --out_dir "$OUT" --plot_dir "$PLOTS"
Rscript "$HERE/plot_converter_compare.R"    --out_dir "$OUT" --plot_dir "$PLOTS"
Rscript "$HERE/plot_converter_timeseries.R" --out_dir "$OUT" --plot_dir "$PLOTS"

echo "Plots in: $PLOTS"
