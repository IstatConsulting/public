# MRP vs. raking under sparsity and an omitted region × age interaction

A Monte Carlo study of when multilevel regression and poststratification (MRP)
outperforms raking for a binary outcome (support for party A), varying sample
size, cell sparsity, and whether each estimator accounts for a region × age
interaction that exists in both the outcome and the selection mechanism.

Working paper (short write-up for colleagues): `paper/MRP_vs_raking_working_paper.docx`.

## Running

```sh
./run_scenarios.sh                             # 2 shapes x 3 selection scenarios + comparison
SHAPES=nonlinear ./run_scenarios.sh            # non-linear outcome interaction only
Rscript run_simulation.R                       # one scenario (baseline, 20 Stan reps)
SIM_SEL_IA=0.5 SIM_IA_SHAPE=nonlinear SIM_N_STAN=0 Rscript run_simulation.R
SIM_R=20 SIM_N_STAN=0 Rscript run_simulation.R # quick check
```

Requires `tidyverse`, `survey`, `lme4`, `furrr`, and optionally `rstanarm`.
One scenario takes about 8 minutes on 4 cores for the main design, plus
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

**Non-linear outcome interaction** (`SIM_IA_SHAPE=nonlinear`). The region ×
age rows of the outcome model are replaced by shapes that are exactly
orthogonal to the linear age contrast (−1.5, −0.5, 0.5, 1.5), with the same
per-region RMS as the linear rows: R2 a mid-age peak (−0.70, 1.10, −0.10,
−0.30), R5 a U-shape with young and old pro-A (0.69, −0.69, −0.69, 0.69), R1 a
small cubic wiggle (−0.10, 0.30, −0.30, 0.10). Region-specific linear age
slopes therefore have nothing to recover. The young cells keep their sign, so
the (unchanged) selection mechanism stays aligned with the interaction where
it was. Votes are redrawn from the same seed; truth: national 0.498, R2 0.500,
R5 0.626.

**Estimators.**
(1) raw mean;
(2) `survey::rake` on the four single-variable margins;
(3) `survey::rake` on the region × age joint margin plus gender and education;
(4) `lme4::glmer` with varying intercepts for age, gender, region, education;
(5) as (4) plus a varying intercept for region × age;
(6) as (5) with gender as a fixed effect;
(6b) as (6) plus a fixed linear trend over ordered age, `gender + age_c +
(1 | age) + (1 | region) + (1 | educ) + (1 | region_age)`, with `age_c` the
centred age-group index, so that `(1 | age)` only absorbs departures from
linearity;
(7) as (6b) plus region-specific linear age slopes, `gender + age_c +
(1 | age) + (1 | educ) + (1 + age_c || region) + (1 | region_age)`. This is a
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

Each scenario writes to `output/sel_ia_<k>/` (linear interaction) or
`output/nonlin_sel_ia_<k>/` (non-linear). `compare_scenarios.R` writes
`output/results_by_scenario.csv`, `output/rmse_grid_by_scenario.csv`,
`output/fig_{rmse,bias}_by_scenario_{linear,nonlinear}.png` and
`output/fig_rmse_by_shape_k1.png`. Per scenario:

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
`output/fig_rmse_by_scenario_linear.png`). Round 3 adds (6b), which changes
the interpretation of (7)'s advantage:

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

**Among (1)–(7), the structured model (7) is the best estimator in 22 of 27
scenario × n × estimand cells (ties counted as best).** Round 3 shows that
most of this gain comes from its fixed linear age trend, not from the
region-specific slopes. Relative to (5), (7) roughly halves the shrinkage bias of (5) in R2
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

## Results, round 3: non-linear interaction and decomposition of (7)

The full design (estimators (1)–(7) plus (6b), three selection strengths) was
run for both outcome-interaction shapes; all previously reported numbers
reproduce exactly. RMSE in percentage points at baseline selection (k = 1);
the full 54-row grid is in `output/rmse_grid_by_scenario.csv` and the
comparison figure is `output/fig_rmse_by_shape_k1.png`:

| Shape | Estimand | n | (3) Rake region × age | (4) MRP int. | (5) MRP + r × a | (6) gender FE | (6b) + age trend | (7) + region slopes |
|---|---|---|---|---|---|---|---|---|
| Linear | National | 200 | 4.71 | 5.27 | 5.28 | 5.37 | 4.43 | **4.42** |
| Linear | R2 | 200 | 11.40 | 9.21 | 8.52 | 8.58 | 6.82 | **6.45** |
| Linear | R2 | 2000 | 3.10 | 8.53 | 4.44 | 4.48 | 4.70 | **2.89** |
| Linear | R5 | 200 | 21.88 | 9.42 | 9.80 | 9.82 | **7.48** | 7.85 |
| Linear | R5 | 2000 | 5.91 | 4.29 | 4.65 | 4.65 | **4.24** | 4.86 |
| Non-linear | National | 200 | 4.50 | 4.66 | 4.79 | 4.86 | **4.22** | 4.23 |
| Non-linear | R2 | 200 | 11.55 | 7.35 | 6.67 | 6.74 | 6.14 | **6.12** |
| Non-linear | R2 | 2000 | 3.26 | 3.22 | 2.16 | 2.16 | **2.08** | **2.08** |
| Non-linear | R5 | 200 | 20.34 | 11.35 | 11.52 | 11.63 | 7.97 | **7.92** |
| Non-linear | R5 | 2000 | 5.11 | 4.99 | 4.14 | 4.15 | **3.95** | **3.95** |

**The fixed age trend does most of the work, under both shapes.** Moving
from (6) to (6b) lowers RMSE in 51 of 54 cells, by 0.4 pp nationally and by
1.2–1.9 pp in R5 on average. The mechanism is the age *main* effect, which
is strong and linear in both DGPs: as a four-level random intercept it is
shrunk towards zero, which biases every cell's prediction, most visibly at
n = 200 (national bias at k = 1 falls from −2.8 to −0.9 pp in the linear DGP
and from −2.1 to 0.0 pp in the non-linear one). The three exceptions are R2
at n = 2,000 under the linear interaction (e.g. 4.48 → 4.70), where the fixed
trend leaves the steep R2 gradient to the exchangeable `(1 | region_age)`
term alone.

**Region-specific slopes help only when the interaction is linear and the
region is well covered.** Under the non-linear interaction, (7) and (6b)
differ by at most 0.08 pp in any cell: the slope variance is estimated near
zero and the slopes cost nothing. Under the linear interaction they improve
R2 in all nine cells (by 0.4, 1.0 and 1.2 pp on average at n = 200, 500 and
2,000) but worsen R5 in all nine (by 0.5–0.9 pp) and leave the national
estimate unchanged. I have not decomposed the R5 penalty further; it is
consistent with R5's slope being estimated from few young respondents and
shrunk towards the common trend, while the slope term takes over variance
that `(1 | region_age)` would otherwise allot to R5's cells.

**The non-linear interaction does not rescue raking.** Across the 54 cells
the best estimator is always an MRP variant: (6b) or (7) in 53, and (4) once
(non-linear, national, n = 2,000, k = 1: 1.33 against 1.34). Interaction-aware
raking (3) has small bias throughout (|bias| ≤ 2.8 pp at k = 1, and ≤ 1.3 pp
at n = 2,000) but
its variance in the undersampled R5 keeps it well behind MRP up to n = 500.
Its closest approach is R2 at n = 2,000 under the linear interaction (3.10
against 2.89 for (7)). The complementary-failure pattern of round 2 also
holds under the non-linear interaction: margin raking (2) is consistent at
k = 0 (R2 bias −0.2 pp at n = 2,000) and biased at k = 1 (−4.0 pp), while
intercept-only MRP (4) is biased in R5 at every k (−3.5 to −4.5 pp at
n = 2,000).

**Revised reading of round 2.** The earlier conclusion that a structured
interaction prior is what makes MRP competitive was too strong. The robust
improvement is to model ordered age with a fixed trend, so that a strong main
effect is not shrunk through a four-level random intercept. Region-specific
slopes are a second-order refinement that pays off only when the interaction
really is a regional age gradient and the region is not undersampled.

**Remaining caveat.** Interval coverage is still not evaluated; it requires
design-based variances for raking (e.g. `svymean` on the raked design, with
and without accounting for trimming) and posterior intervals for MRP on all
replications, i.e. Stan at roughly 30 s per fit.
