#!/usr/bin/env bash
# Run the main design under three strengths of the region x age selection term.
# The rstanarm subset is run only for the baseline (set SIM_N_STAN to override).
set -euo pipefail
cd "$(dirname "$0")"
for k in 1 0.5 0; do
  dir=$(printf "output/sel_ia_%.1f" "$k")
  mkdir -p "$dir"
  SIM_SEL_IA=$k SIM_N_STAN=${SIM_N_STAN:-0} Rscript run_simulation.R > "$dir/run_log.txt" 2>&1
  echo "done: $dir"
done
Rscript compare_scenarios.R
