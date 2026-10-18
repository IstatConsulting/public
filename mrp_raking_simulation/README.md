# MRP vs. raking under sparsity and an omitted region × age interaction

A Monte Carlo study of when multilevel regression and poststratification (MRP)
outperforms raking for a binary outcome (support for party A), varying sample
size, cell sparsity, and whether each estimator accounts for a region × age
interaction that exists in both the outcome and the selection mechanism.

## Running

```sh
./run_scenarios.sh                             # three selection scenarios + comparison
Rscript run_simulation.R                       # one scenario (baseline, 20 Stan reps)
SIM_SEL_IA=0.5 SIM_N_STAN=0 Rscript run_simulation.R
SIM_R=20 SIM_N_STAN=0 Rscript run_simulation.R # quick check
```

Requires `tidyverse`, `survey`, `lme4`, `furrr`, and optionally `rstanarm`.
One scenario takes about 7 minutes on 4 cores for the main design, plus
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
is not. In the scenario runs, the region × age selection term is multiplied
by k ∈ {1, 0.5, 0} (`SIM_SEL_IA`); the population, its votes and the truth
are identical across scenarios. At k = 0 selection is log-additive in the
four margins.

**Estimators.**
(1) raw mean;
(2) `survey::rake` on the four single-variable margins;
(3) `survey::rake` on the region × age joint margin plus gender and education;
(4) `lme4::glmer` with varying intercepts for age, gender, region, education;
(5) as (4) plus a varying intercept for region × age;
(6) as (5) with gender as a fixed effect;
(7) as (6) plus region-specific linear age slopes around a fixed common
slope, `gender + age_c + (1 | age) + (1 | educ) + (1 + age_c || region) +
(1 | region_age)`, with `age_c` the centred age-group index. This is a
frequentist stand-in for a structured prior over ordered age (Gao et al.
2021): the interaction is pooled along a linear age trend within region,
and `(1 | region_age)` absorbs residual non-linearity.
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

## Outputs

Each scenario writes to `output/sel_ia_<k>/`; `compare_scenarios.R` writes
`output/results_by_scenario.csv`, `output/rmse_grid_by_scenario.csv` and
`output/fig_{rmse,bias}_by_scenario.png`. Per scenario:

| File | Content |
|---|---|
| `results_table.csv` | Tidy results: bias, empirical SD, RMSE per n × estimator × estimand, joined with estimator diagnostics |
| `results_diagnostics.csv` | Kish deff, nonconvergence, failure and collapse rates (raking); singular-fit and convergence-warning rates (MRP) |
| `results_sparsity.csv` | Empty cells, cells with 1–2 respondents, empty region × age cells |
| `results_stan_check.csv` | glmer vs. stan_glmer on the identical replication subset (baseline only) |
| `fig_rmse.png`, `fig_bias.png` | RMSE and bias by estimator and n, faceted by estimand |
| `pattern_check.txt` | RMSE grid and the expected-pattern checks |
| `replication_*.csv*` | Replication-level estimates, diagnostics, sparsity |
| `population_cells.csv`, `truth.csv` | Joint cell counts, true cell probabilities, true support |

## Results, round 1: estimators (1)–(5), baseline selection (k = 1)

RMSE in percentage points (500 replications per condition; bias in
`results_table.csv`, figure in `output/sel_ia_1.0/fig_rmse.png`):

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
small bias of (4) turns out to be a cancellation of two opposing biases
(see round 2: it becomes +4.8 pp once the selection interaction is
removed).

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

## Results, round 2: structured MRP and the strength of the selection interaction

Estimators (6) and (7) were added and the full design was rerun under
k = 1, 0.5 and 0 (500 replications each, same seed; estimators (1)–(5)
reproduce round 1 exactly at k = 1). RMSE in percentage points, raw mean
omitted (full grid in `output/rmse_grid_by_scenario.csv`, figure in
`output/fig_rmse_by_scenario.png`):

| k | Estimand | n | (2) Rake margins | (3) Rake region × age | (4) MRP int. | (5) MRP + r × a | (6) (5), gender FE | (7) Structured |
|---|---|---|---|---|---|---|---|---|
| 1 | National | 200 | 4.80 | 4.71 | 5.27 | 5.28 | 5.37 | **4.42** |
| 1 | National | 2000 | 1.78 | 1.48 | 1.54 | 1.56 | 1.58 | **1.45** |
| 1 | R2 | 200 | 12.95 | 11.40 | 9.21 | 8.52 | 8.58 | **6.45** |
| 1 | R2 | 2000 | 8.67 | 3.10 | 8.53 | 4.44 | 4.48 | **2.89** |
| 1 | R5 | 200 | 21.98 | 21.88 | 9.42 | 9.80 | 9.82 | **7.85** |
| 1 | R5 | 2000 | 6.29 | 5.91 | **4.29** | 4.65 | 4.65 | 4.86 |
| 0.5 | R2 | 2000 | 5.40 | 2.87 | 6.87 | 3.52 | 3.56 | **2.48** |
| 0.5 | R5 | 2000 | 6.19 | 6.00 | **4.64** | 4.66 | 4.67 | 5.14 |
| 0 | R2 | 200 | 10.24 | 9.75 | 6.44 | 6.27 | 6.30 | **5.25** |
| 0 | R2 | 2000 | 3.31 | 2.78 | 4.99 | 2.81 | 2.84 | **2.28** |
| 0 | R5 | 2000 | 6.69 | 6.74 | 6.92 | **5.40** | 5.40 | 5.87 |

**Gender as a fixed effect (6) changes nothing of substance.** RMSE differs
from (5) by at most 0.1 pp in every cell; the singular-fit rate falls (e.g.
74% → 51% at n = 500, k = 1). The zero gender variance in round 1 was a
cosmetic boundary estimate, not a source of bias.

**The structured interaction (7) is the best estimator in 22 of 27
scenario × n × estimand cells (ties counted as best).** Pooling the region × age deviations along a
region-specific age trend roughly halves the shrinkage bias of (5) in R2
(k = 1, n = 2,000: −1.7 vs. −3.6 pp) and lowers its variance, so (7) beats
interaction-aware raking (3) in R2 even at n = 2,000, and is best (once tied)
nationally at every n and k. All five exceptions are R5 at n ≥ 500, where (7) retains a
bias of −0.9 to −2.3 pp: R5 is undersampled and its young respondents are rare, so
its age slope, estimated from only five regions, is still pulled towards the
common slope (median over 30 baseline samples at n = 2,000: −0.45 against
a true −0.61). (7) is
singular in 90–100% of fits; this is expected rather than a symptom, because
in this DGP the age main effect and most of the interaction are linear in
ordered age, so the true residual `(1 | age)` and `(1 | region_age)` variances
are close to zero.

**Raking and MRP fail under complementary conditions.** Scaling the selection
interaction isolates the mechanism. Margin raking (2) is consistent whenever
the selection propensity is log-linear in the raked margins, whatever the
outcome model: its R2 bias at n = 2,000 goes from −8.3 (k = 1) to −4.5
(k = 0.5) to +0.2 pp (k = 0). Intercept-only MRP (4) is consistent whenever
the outcome model is correct, whatever the selection: its R2 bias stays at
−8.2, −6.4 and −4.4 pp, because even purely marginal over-selection of the
young leaves R2's within-region age composition unbalanced, and the omitted
interaction then biases the regional prediction. In R5 the two sources of
bias in (4) have opposite signs, so its small baseline bias (−0.2 pp) is a
cancellation that becomes +4.8 pp at k = 0. The practical implication is that
(4)'s good showing in R5 at k ≥ 0.5 should not be read as robustness.

**Caveats.** The DGP favours (7) by construction: its interaction is exactly
linear in age in R1 and nearly linear in R2 and R5, which is precisely the
structure (7) pools towards. A non-linear interaction (e.g. a mid-age
peak in one region) is the condition that would test whether the structured
model's advantage survives misspecification of the structure. Interval
coverage is still not evaluated; it requires design-based variances for
raking (e.g. `svymean` on the raked design, with and without accounting for
trimming) and posterior intervals for MRP on all replications, i.e. Stan at
roughly 30 s per fit.
