# MRP vs. raking under sparsity and an omitted region × age interaction

A Monte Carlo study of when multilevel regression and poststratification (MRP)
outperforms raking for a binary outcome (support for party A), varying sample
size, cell sparsity, and whether each estimator accounts for a region × age
interaction that exists in both the outcome and the selection mechanism.

## Running

```sh
Rscript run_simulation.R                       # full design (500 reps, 20 Stan reps)
SIM_R=20 SIM_N_STAN=0 Rscript run_simulation.R # quick check
```

Requires `tidyverse`, `survey`, `lme4`, `furrr`, and optionally `rstanarm`.
The full run takes about 15 minutes on 4 cores for the main design, plus
about 30 minutes for the Stan subset. All randomness (population and all
samples) is fixed by `seed <- 20261018` before estimation; every estimator
receives the identical sample in each replication.

## Design

**Population** (`simulate_population()`). N = 50,000 across age (4) × gender
(2) × region (5) × education (3) = 120 cells. Cell shares combine uneven
region sizes (R2 = 8%, R5 = 10%), region-specific age profiles (R2 young, R5
old), age-dependent education, and lognormal noise (SD 0.5), so several cells
hold fewer than 50 people. True support follows a logistic model with main
effects for all four variables plus a region × age interaction: a flatter age
gradient in R1, a steep age polarisation in R2 (young strongly against A), and
a reversed gradient in R5 (young pro-A). Votes are Bernoulli draws; the truth
is the finite-population mean (national 0.479; R2 0.420; R5 0.530).

**Selection** (`draw_sample_ids()`). Unequal-probability sampling without
replacement, with log-propensity additive in age (young oversampled), region
(R2, R3 oversampled; R5 undersampled), education and gender, *plus* a region ×
age term aligned with the outcome interaction: anti-A young in R2 are further
oversampled, pro-A young in R5 further undersampled. Within the full 120-cell
cross-classification selection is ignorable; within single-variable margins it
is not.

**Estimators.**
(1) raw mean;
(2) `survey::rake` on the four single-variable margins;
(3) `survey::rake` on the region × age joint margin plus gender and education;
(4) `lme4::glmer` with varying intercepts for age, gender, region, education;
(5) as (4) plus a varying intercept for region × age.
Both raking estimators use `maxit = 100`, `epsilon = 1e-6` (relative), and
strict upper trimming at 5 × the mean weight. MRP is fitted to the
cell-aggregated binomial data, predicted for all 120 cells
(`allow.new.levels = TRUE`), and poststratified with the true joint counts.
The rstanarm refit of (4) and (5) uses 4 chains × 2,000 iterations,
`adapt_delta = 0.99`, default priors, on the first 20 replications of each
condition.

**Empty region × age cells in (3).** `postStratify()` errors whenever a
population stratum is absent from the sample, so estimator (3) as specified
is infeasible for such samples. Rather than dropping those replications, each
empty cell is collapsed into the adjacent age group within the same region
(`collapse_region_age()`), and the replication is flagged
(`pct_collapsed`, `mean_empty_ra`). Any hard failure would be retained with
`NA` estimates and counted in `pct_failed`; none occurred.

## Outputs (`output/`)

| File | Content |
|---|---|
| `results_table.csv` | Tidy results: bias, empirical SD, RMSE per n × estimator × estimand, joined with estimator diagnostics |
| `results_diagnostics.csv` | Kish deff, nonconvergence, failure and collapse rates (raking); singular-fit and convergence-warning rates (MRP) |
| `results_sparsity.csv` | Empty cells, cells with 1–2 respondents, empty region × age cells |
| `results_stan_check.csv` | glmer vs. stan_glmer on the identical replication subset |
| `fig_rmse.png`, `fig_bias.png` | RMSE and bias by estimator and n, faceted by estimand |
| `pattern_check.txt` | RMSE grid and the expected-pattern checks |
| `replication_*.csv*` | Replication-level estimates, diagnostics, sparsity |
| `population_cells.csv`, `truth.csv` | Joint cell counts, true cell probabilities, true support |

## Results

RMSE in percentage points (500 replications per condition; bias in
`results_table.csv`, figure in `output/fig_rmse.png`):

| Estimand | n | (1) Raw | (2) Rake margins | (3) Rake region × age | (4) MRP int. | (5) MRP + region × age |
|---|---|---|---|---|---|---|
| National | 200 | 10.03 | 4.80 | **4.71** | 5.27 | 5.28 |
| National | 500 | 9.76 | 3.08 | **2.88** | 3.16 | 3.21 |
| National | 2000 | 9.27 | 1.78 | **1.48** | 1.54 | 1.56 |
| R2 (small, oversampled) | 200 | 18.30 | 12.98 | 11.36 | 9.21 | **8.52** |
| R2 (small, oversampled) | 500 | 16.98 | 9.98 | **6.25** | 8.31 | 6.83 |
| R2 (small, oversampled) | 2000 | 16.38 | 8.67 | **3.10** | 8.53 | 4.44 |
| R5 (small, undersampled) | 200 | 20.30 | 22.04 | 21.90 | **9.42** | 9.80 |
| R5 (small, undersampled) | 500 | 12.42 | 13.47 | 13.25 | **6.86** | 7.34 |
| R5 (small, undersampled) | 2000 | 5.75 | 6.29 | 5.91 | **4.29** | 4.65 |

Sparsity: at n = 200 a sample leaves on average 49.9 of the 120 cells empty
and 42.3 with 1–2 respondents; 84.6% of samples have at least one empty
region × age cell (mean 1.5 of 20). At n = 500 these figures are 27.0, 36.0
and 22.4%; at n = 2,000, 6.2, 17.1 and 0%. No raking run failed or failed to
converge; Kish deff was 1.74–1.85 for (2) and 1.76–1.91 for (3). Singular
fits occurred in 77.6% / 49.6% / 16.2% of (4) fits and 92.8% / 74.4% / 41.2%
of (5) fits at n = 200 / 500 / 2,000; no glmer convergence warnings.

**Against the expected pattern.**

*Interaction-aware raking (3) matches or beats MRP at n = 2,000*: confirmed
nationally and in R2, where it is the best estimator outright (3.10 vs. 4.44
for MRP (5)). It does not hold in R5: an undersampled region still yields
few respondents at n = 2,000, and the weighting variance (empirical SD 5.9 pp
vs. 4.3 for MRP) dominates.

*(3) degrades at n = 200*: confirmed for the regional estimands (R2 11.4 vs.
8.5 for MRP (5); R5 21.9 vs. 9.4–9.8), but not nationally, where (3) remains
best. The national result needs qualifying: strictly as specified, (3) is
not computable in 84.6% of n = 200 samples, and the figure reflects the
collapsed fallback, which in sparse regions reverts towards margin raking.
Degradation at n = 200 is therefore primarily a feasibility problem and a
regional variance problem, not a national RMSE problem.

*Intercept-only MRP (4) biased where the interaction matters*: confirmed for
R2, where bias is about −6 to −8 pp and does not shrink with n (RMSE 8.5 at
n = 2,000, no better than margin raking). It is not visible in R5 or
nationally, where (4) and (5) are practically indistinguishable; in R5 the
omitted interaction produces little region-level bias. I have not
decomposed why; a plausible reading is that, after poststratification with
the true R5 age composition, the misfit of the age gradient partly cancels
across R5's age groups.

*MRP with the interaction (5) best or near best across conditions*: only
partly. (5) is best or within ~10% of the best for regional estimands at
n ≤ 500 and in R5 throughout, but it is beaten clearly by (3) in R2 at
n ≥ 500 and is the worst non-raw estimator nationally at n = 200. The
reason is over-smoothing of the region × age effects. In this DGP, the
interaction is concentrated in 8 of 20 region × age cells with large
deviations (true RMS deviation ≈ 0.45 on the logit scale), and the most
important cell (R2 × 60+, true support 0.83) has about five respondents even
at n = 2,000 because three selection terms work against it. A single
exchangeable normal variance for 20 cells is a poor prior for this
structure: on the same samples, the estimated region × age SD has a median
of 0.24 at n = 200 and 0.38 at n = 2,000, and it is estimated at exactly
zero in 35% of n = 200 fits, which reduces (5) to (4). The resulting bias
(−3.6 pp in R2 and −0.6 pp nationally at n = 2,000) persists as n grows.
Most singular fits come from the two-level gender intercept (variance zero
in 55–67% of fits at n ≤ 500).

**rstanarm robustness check** (20 replications per condition,
`results_stan_check.csv`). `stan_glmer` reproduces the glmer results
closely: RMSE differs by at most 0.8 pp, and the residual biases of (4) in R2
and of (5) in R2 are virtually identical (e.g. −3.72 vs. −3.61 pp for (5) in
R2 at n = 2,000). Integrating over the variance components therefore does not
remove the shrinkage bias; it is a property of the model, not of
glmer's boundary estimates. Sampler diagnostics were imperfect: 15–50% of
fits had at least one divergent transition despite `adapt_delta = 0.99`, and
one (5) fit at n = 500 had R-hat 1.31.

**Implications and extensions.** The data-generating process is close to a
worst case for exchangeable random intercepts: a sparse, heavy-tailed
interaction, aligned with selection, in the cells the sample covers worst.
Natural next conditions are (a) a structured prior for the interaction,
e.g. region-specific age slopes `(1 + age_num | region)` or the
autoregressive / random-walk priors over ordered age of Gao et al. (2021,
*Bayesian Analysis*); (b) gender (and possibly education) as fixed effects,
which removes most singular fits; (c) a milder selection interaction, to
separate the effect of the interaction's size from its overlap with
selection; and (d) coverage of interval estimates, which this design does
not evaluate.
