#!/usr/bin/env bash
# Run the main design under three strengths of the region x age selection term,
# for the linear and the non-linear outcome interaction (SHAPES to restrict,
# e.g. SHAPES=nonlinear). The rstanarm subset is off unless SIM_N_STAN is set.
set -euo pipefail
cd "$(dirname "$0")"
for shape in ${SHAPES:-linear nonlinear}; do
  prefix=$([ "$shape" = nonlinear ] && echo "nonlin_" || echo "")
  for k in 1 0.5 0; do
    dir=$(printf "output/%ssel_ia_%.1f" "$prefix" "$k")
    mkdir -p "$dir"
    SIM_IA_SHAPE=$shape SIM_SEL_IA=$k SIM_N_STAN=${SIM_N_STAN:-0} \
      Rscript run_simulation.R > "$dir/run_log.txt" 2>&1
    echo "done: $dir"
  done
done
Rscript compare_scenarios.R
