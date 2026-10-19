# =============================================================================
# MRP vs. raking simulation: functions
# -----------------------------------------------------------------------------
# Population generation, biased sampling, the five estimators (+ optional
# rstanarm refits), and evaluation helpers. Sourced by run_simulation.R.
# =============================================================================

suppressPackageStartupMessages({
  library(tidyverse)
  library(survey)
  library(lme4)
})

# -----------------------------------------------------------------------------
# 0. Design constants
# -----------------------------------------------------------------------------

age_lv    <- c("18-29", "30-44", "45-59", "60+")
gender_lv <- c("Male", "Female")
region_lv <- paste0("R", 1:5)
educ_lv   <- c("Low", "Mid", "High")

# R2 = small, oversampled region; R5 = small, undersampled region.
focus_regions <- c(R2 = "R2 (small, oversampled)", R5 = "R5 (small, undersampled)")

# Population composition -------------------------------------------------------
pop_params <- list(
  region_share = c(R1 = .30, R2 = .08, R3 = .27, R4 = .25, R5 = .10),
  # age distribution within region (rows = region): R2 young, R5 old
  age_by_region = rbind(
    R1 = c(.22, .30, .25, .23),
    R2 = c(.34, .28, .22, .16),
    R3 = c(.18, .26, .28, .28),
    R4 = c(.17, .25, .28, .30),
    R5 = c(.10, .20, .30, .40)
  ),
  gender_share = c(.48, .52),
  # education distribution within age group (rows = age)
  educ_by_age = rbind(
    c(.20, .50, .30),
    c(.20, .45, .35),
    c(.30, .48, .22),
    c(.42, .43, .15)
  ),
  cell_noise_sd = 0.5 # lognormal noise on cell shares -> some very rare cells
)

# True outcome model (logit scale) ---------------------------------------------
outcome_params <- list(
  intercept = -0.10,
  age       = c(-0.60, -0.20, 0.20, 0.60),
  gender    = c(-0.10, 0.10),
  region    = c(-0.40, 0.10, 0.10, 0.20, 0.30),
  educ      = c(0.40, 0.00, -0.50),
  # region x age interaction (rows = region): R1 flatter age gradient,
  # R2 steeply age-polarised, R5 reversed gradient (young pro-A)
  region_age = rbind(
    R1 = c(0.30, 0.10, -0.10, -0.30),
    R2 = c(-0.90, -0.40, 0.30, 0.90),
    R3 = c(0, 0, 0, 0),
    R4 = c(0, 0, 0, 0),
    R5 = c(1.00, 0.40, -0.30, -0.80)
  )
)

#' Non-linear alternative for the region x age interaction
#'
#' Each row is exactly orthogonal to the linear age contrast (-1.5, -0.5, 0.5,
#' 1.5), so region-specific linear age slopes cannot capture any of it, and
#' has (approximately) the same RMS as the corresponding linear row: R2 a
#' mid-age peak at 30-44, R5 a U-shape (young and old pro-A), R1 a small
#' cubic wiggle. The young cells keep their sign, so selection stays aligned.
region_age_nonlinear <- rbind(
  R1 = c(-0.10, 0.30, -0.30, 0.10),
  R2 = c(-0.70, 1.10, -0.10, -0.30),
  R3 = c(0, 0, 0, 0),
  R4 = c(0, 0, 0, 0),
  R5 = c(0.69, -0.69, -0.69, 0.69)
)

#' Outcome parameters with a linear (baseline) or non-linear interaction
outcome_with_shape <- function(shape = c("linear", "nonlinear"), op = outcome_params) {
  shape <- match.arg(shape)
  if (shape == "nonlinear") op$region_age <- region_age_nonlinear
  op
}

# Selection (inclusion propensity, log scale) ----------------------------------
# Younger groups and regions R2, R3 oversampled; R5 undersampled. The region x
# age term is aligned with the outcome interaction: anti-A young in R2 are
# over-represented and pro-A young in R5 are under-represented, so a weighting
# scheme that only matches single-variable margins cannot remove the bias.
selection_params <- list(
  age    = c(0.90, 0.40, 0.00, -0.40),
  gender = c(0.00, 0.10),
  region = c(0.00, 0.90, 0.70, 0.00, -0.60),
  educ   = c(-0.40, 0.00, 0.40),
  region_age = rbind(
    R1 = c(0, 0, 0, 0),
    R2 = c(0.70, 0.30, -0.20, -0.60),
    R3 = c(0, 0, 0, 0),
    R4 = c(0, 0, 0, 0),
    R5 = c(-0.70, -0.30, 0.20, 0.50)
  )
)

#' Scale the region x age selection term (1 = baseline, 0 = margins-only bias)
#'
#' The population (including votes) is unaffected: log_incl is deterministic,
#' so the same seed yields the same population under every scaling.
scale_selection_ia <- function(sp = selection_params, k = 1) {
  sp$region_age <- k * sp$region_age
  sp
}

# Raking settings, shared by estimators (2) and (3) ----------------------------
rake_settings <- list(
  control  = list(maxit = 100, epsilon = 1e-6, verbose = FALSE),
  trim_cap = 5 # upper trim at 5 x mean weight, strict (iterative) trimming
)

# -----------------------------------------------------------------------------
# 1. Population
# -----------------------------------------------------------------------------

#' Full 120-cell frame with level indices (no counts)
make_cell_frame <- function() {
  expand_grid(
    region = factor(region_lv, levels = region_lv),
    age    = factor(age_lv, levels = age_lv),
    gender = factor(gender_lv, levels = gender_lv),
    educ   = factor(educ_lv, levels = educ_lv)
  ) |>
    mutate(
      cell       = str_c(region, age, gender, educ, sep = "|"),
      region_age = str_c(region, age, sep = ":"),
      across(c(region, age, gender, educ), as.integer, .names = "{.col}_i")
    )
}

#' Simulate the finite population
#'
#' Returns a list with individual-level data, joint cell counts with true
#' cell-level support, and true national / regional support.
simulate_population <- function(N = 50000,
                                pp = pop_params,
                                op = outcome_params,
                                sp = selection_params) {
  cells <- make_cell_frame() |>
    mutate(
      share = pp$region_share[region_i] *
        pp$age_by_region[cbind(region_i, age_i)] *
        pp$gender_share[gender_i] *
        pp$educ_by_age[cbind(age_i, educ_i)] *
        exp(rnorm(n(), 0, pp$cell_noise_sd)),
      share = share / sum(share),
      N_c   = as.vector(rmultinom(1, N, share)),
      eta   = op$intercept + op$age[age_i] + op$gender[gender_i] +
        op$region[region_i] + op$educ[educ_i] +
        op$region_age[cbind(region_i, age_i)],
      p_true = plogis(eta),
      log_incl = sp$age[age_i] + sp$gender[gender_i] + sp$region[region_i] +
        sp$educ[educ_i] + sp$region_age[cbind(region_i, age_i)]
    )

  persons <- cells |>
    select(cell, region, age, gender, educ, region_age, p_true, log_incl, N_c) |>
    uncount(N_c) |>
    mutate(
      id   = row_number(),
      vote = rbinom(n(), 1, p_true),
      incl = exp(log_incl)
    )

  # Realised (finite-population) counts and support
  cells <- cells |>
    left_join(
      persons |> summarise(y_c = sum(vote), .by = cell),
      by = "cell"
    ) |>
    mutate(y_c = replace_na(y_c, 0L))

  truth <- bind_rows(
    persons |> summarise(truth = mean(vote)) |> mutate(estimand = "National"),
    persons |>
      summarise(truth = mean(vote), .by = region) |>
      mutate(estimand = as.character(region)) |>
      select(-region)
  )

  list(persons = persons, cells = cells, truth = truth, N = N)
}

#' Population margins in survey::rake format
make_pop_margins <- function(cells) {
  list(
    age        = count(cells, age, wt = N_c, name = "Freq"),
    gender     = count(cells, gender, wt = N_c, name = "Freq"),
    region     = count(cells, region, wt = N_c, name = "Freq"),
    educ       = count(cells, educ, wt = N_c, name = "Freq"),
    region_age = count(cells, region, age, region_age, wt = N_c, name = "Freq")
  )
}

# -----------------------------------------------------------------------------
# 2. Sampling
# -----------------------------------------------------------------------------

#' Draw sample ids for all replications of one sample size
#'
#' Unequal-probability sampling without replacement (successive draws with
#' probability proportional to the inclusion propensity).
draw_sample_ids <- function(pop, n, R) {
  map(seq_len(R), \(r) sample.int(pop$N, n, replace = FALSE, prob = pop$persons$incl))
}

#' Sparsity diagnostics for one sample
sparsity_stats <- function(samp, cells) {
  n_c  <- tabulate(match(samp$cell, cells$cell), nbins = nrow(cells))
  n_ra <- table(factor(samp$region_age, levels = unique(cells$region_age)))
  tibble(
    empty_cells       = sum(n_c == 0),
    cells_1_2         = sum(n_c %in% 1:2),
    empty_region_age  = sum(n_ra == 0),
    empty_region      = sum(table(factor(samp$region, levels = region_lv)) == 0)
  )
}

# -----------------------------------------------------------------------------
# 3. Estimators
# -----------------------------------------------------------------------------

#' Point estimates (national + each region) from a vector of predictions or
#' weights. `w` are unit-level weights; `y` the outcome.
weighted_estimates <- function(y, w, region) {
  reg <- tibble(y, w, region) |>
    summarise(estimate = sum(w * y) / sum(w), .by = region) |>
    mutate(estimand = as.character(region)) |>
    select(estimand, estimate)
  bind_rows(
    tibble(estimand = "National", estimate = sum(w * y) / sum(w)),
    reg
  ) |>
    complete(estimand = c("National", region_lv))
}

kish_deff <- function(w) length(w) * sum(w^2) / sum(w)^2

# (1) Raw mean -------------------------------------------------------------------
est_raw <- function(samp) {
  list(
    estimates = weighted_estimates(samp$vote, rep(1, nrow(samp)), samp$region),
    diag      = tibble(deff = 1)
  )
}

#' Generic raking wrapper around survey::rake + survey::trimWeights
#'
#' Captures nonconvergence (warning) and hard failures (error) instead of
#' dropping the replication. `margin_vars` names the sample variables, and
#' `pop_margins` the matching population tables.
rake_core <- function(samp, margin_vars, pop_margins, N,
                      settings = rake_settings) {
  samp$w0 <- N / nrow(samp)
  des <- svydesign(ids = ~1, weights = ~w0, data = samp)
  converged <- TRUE

  raked <- tryCatch(
    withCallingHandlers(
      rake(des,
        sample.margins = map(margin_vars, \(v) reformulate(v)),
        population.margins = pop_margins,
        control = settings$control
      ),
      warning = function(w) {
        if (str_detect(conditionMessage(w), "did not converge")) {
          converged <<- FALSE
          invokeRestart("muffleWarning")
        }
      }
    ),
    error = function(e) e
  )

  if (inherits(raked, "error")) {
    return(list(
      estimates = weighted_estimates(samp$vote, rep(NA_real_, nrow(samp)), samp$region),
      diag = tibble(deff = NA_real_, converged = NA, failed = TRUE,
                    error = conditionMessage(raked))
    ))
  }

  w_raked <- weights(raked)
  trimmed <- trimWeights(raked, upper = settings$trim_cap * mean(w_raked), strict = TRUE)
  w <- weights(trimmed)

  list(
    estimates = weighted_estimates(samp$vote, w, samp$region),
    diag = tibble(
      deff = kish_deff(w), converged = converged, failed = FALSE,
      error = NA_character_, n_trimmed = sum(w_raked > settings$trim_cap * mean(w_raked))
    )
  )
}

# (2) Raking on the four single-variable margins --------------------------------
est_rake_margins <- function(samp, pop_margins, N) {
  pm <- pop_margins[c("age", "gender", "region", "educ")]
  out <- rake_core(samp, c("age", "gender", "region", "educ"), pm, N)
  out$diag <- out$diag |> mutate(empty_ra_collapsed = 0L)
  out
}

#' Collapse empty region x age cells into an adjacent age group within region
#'
#' Within each region, age groups are accumulated in order until the running
#' group contains at least one respondent; trailing empty groups are merged
#' into the preceding group. Returns a lookup region_age -> raking cell.
collapse_region_age <- function(samp) {
  counts <- samp |>
    count(region, age, .drop = FALSE) |>
    complete(region, age, fill = list(n = 0L)) |>
    arrange(region, age)

  counts |>
    group_by(region) |>
    group_modify(\(d, key) {
      g <- integer(nrow(d)); grp <- 1L; acc <- 0L
      for (a in seq_len(nrow(d))) {
        g[a] <- grp
        acc <- acc + d$n[a]
        if (acc >= 1 && a < nrow(d)) { grp <- grp + 1L; acc <- 0L }
      }
      # trailing empty group -> merge into previous one
      last <- g == max(g)
      if (sum(d$n[last]) == 0 && max(g) > 1) g[last] <- max(g) - 1L
      tibble(age = d$age, grp = g)
    }) |>
    ungroup() |>
    mutate(
      region_age = str_c(region, age, sep = ":"),
      ra_cell = str_c(region, "_g", grp)
    ) |>
    select(region_age, ra_cell)
}

# (3) Raking on region x age joint + gender + education --------------------------
est_rake_joint <- function(samp, pop_margins, N) {
  n_empty <- sum(table(factor(samp$region_age, levels = pop_margins$region_age$region_age)) == 0)

  # Fallback for empty joint cells: collapse within region (flagged)
  lookup <- collapse_region_age(samp)
  samp <- samp |> left_join(lookup, by = "region_age")
  pm_ra <- pop_margins$region_age |>
    left_join(lookup, by = "region_age") |>
    count(ra_cell, wt = Freq, name = "Freq")

  out <- rake_core(
    samp, c("ra_cell", "gender", "educ"),
    list(pm_ra, pop_margins$gender, pop_margins$educ), N
  )
  out$diag <- out$diag |> mutate(empty_ra_collapsed = n_empty)
  out
}

#' Aggregate a sample to observed cells (binomial counts)
aggregate_sample <- function(samp) {
  samp |>
    summarise(y = sum(vote), n = n(),
              .by = c(cell, region, age, gender, educ, region_age)) |>
    mutate(across(c(region, age, gender, educ), as.character),
           age_c = match(age, age_lv) - 2.5)
}

#' Prediction frame for all 120 cells, matching aggregate_sample() columns
mrp_newdata <- function(cells) {
  cells |>
    select(region, age, gender, educ, region_age) |>
    mutate(across(c(region, age, gender, educ), as.character),
           age_c = match(age, age_lv) - 2.5)
}

mrp_formulas <- list(
  mrp_int = cbind(y, n - y) ~ 1 + (1 | age) + (1 | gender) + (1 | region) + (1 | educ),
  mrp_ia  = cbind(y, n - y) ~ 1 + (1 | age) + (1 | gender) + (1 | region) + (1 | educ) +
    (1 | region_age),
  # (6) as (5), gender as fixed effect (two-level variance is unidentifiable)
  mrp_ia_fe = cbind(y, n - y) ~ gender + (1 | age) + (1 | region) + (1 | educ) +
    (1 | region_age),
  # (6b) as (6), plus a fixed linear trend over ordered age; (1 | age) then
  # only captures departures from linearity (decomposes (6) -> (7))
  mrp_ia_fe_lin = cbind(y, n - y) ~ gender + age_c + (1 | age) + (1 | region) +
    (1 | educ) + (1 | region_age),
  # (7) structured interaction: region-specific linear age slopes (uncorrelated
  # with region intercepts) around a fixed common slope, plus residual
  # region x age intercepts; (1 | age) captures non-linearity in age
  mrp_struct = cbind(y, n - y) ~ gender + age_c + (1 | age) + (1 | educ) +
    (1 + age_c || region) + (1 | region_age)
)

#' Poststratify cell predictions with the true joint counts
poststratify <- function(pred, cells) {
  tibble(pred, N_c = cells$N_c, region = as.character(cells$region)) |>
    (\(d) bind_rows(
      tibble(estimand = "National", estimate = sum(d$pred * d$N_c) / sum(d$N_c)),
      d |> summarise(estimate = sum(pred * N_c) / sum(N_c), .by = region) |>
        rename(estimand = region)
    ))()
}

# (4) / (5) MRP via lme4::glmer --------------------------------------------------
est_mrp <- function(samp_agg, cells, formula) {
  conv_warn <- FALSE
  fit <- tryCatch(
    withCallingHandlers(
      glmer(formula, data = samp_agg, family = binomial,
            control = glmerControl(optimizer = "bobyqa",
                                   optCtrl = list(maxfun = 2e4))),
      warning = function(w) {
        if (str_detect(conditionMessage(w), "converge|Hessian|gradient")) conv_warn <<- TRUE
        invokeRestart("muffleWarning")
      },
      message = function(m) invokeRestart("muffleMessage") # "boundary (singular)"
    ),
    error = function(e) e
  )

  if (inherits(fit, "error")) {
    return(list(
      estimates = tibble(estimand = c("National", region_lv), estimate = NA_real_),
      diag = tibble(singular = NA, conv_warning = NA, failed = TRUE,
                    error = conditionMessage(fit))
    ))
  }

  newdata <- mrp_newdata(cells)
  pred <- predict(fit, newdata = newdata, type = "response", allow.new.levels = TRUE)

  list(
    estimates = poststratify(pred, cells),
    diag = tibble(singular = isSingular(fit), conv_warning = conv_warn,
                  failed = FALSE, error = NA_character_)
  )
}

# Optional robustness check: rstanarm::stan_glmer ---------------------------------
est_mrp_stan <- function(samp_agg, cells, formula, seed,
                         chains = 4, iter = 2000, adapt_delta = 0.99) {
  low_ess <- FALSE
  fit <- withCallingHandlers(
    rstanarm::stan_glmer(
      formula, data = samp_agg, family = binomial(),
      chains = chains, iter = iter, cores = 1, seed = seed,
      adapt_delta = adapt_delta, refresh = 0
    ),
    warning = function(w) {
      if (str_detect(conditionMessage(w), "Effective Samples Size")) {
        low_ess <<- TRUE
        invokeRestart("muffleWarning")
      }
    }
  )
  newdata <- mrp_newdata(cells)
  # draws x cells; new region_age levels are drawn from the fitted distribution
  epred <- rstanarm::posterior_epred(fit, newdata = newdata)
  N_c <- cells$N_c
  reg <- as.character(cells$region)

  nat_draws <- as.vector(epred %*% N_c) / sum(N_c)
  reg_draws <- map(region_lv, \(r) {
    idx <- reg == r
    as.vector(epred[, idx, drop = FALSE] %*% N_c[idx]) / sum(N_c[idx])
  })

  sp <- rstan::get_sampler_params(fit$stanfit, inc_warmup = FALSE)
  list(
    estimates = tibble(
      estimand = c("National", region_lv),
      estimate = c(mean(nat_draws), map_dbl(reg_draws, mean))
    ),
    diag = tibble(
      divergent = sum(map_dbl(sp, \(x) sum(x[, "divergent__"]))),
      low_ess   = low_ess,
      max_rhat  = max(rstanarm:::summary.stanreg(fit)[, "Rhat"], na.rm = TRUE)
    )
  )
}

# -----------------------------------------------------------------------------
# 4. One replication: same sample -> all estimators
# -----------------------------------------------------------------------------

run_one_rep <- function(ids, pop, pop_margins) {
  samp  <- pop$persons[ids, ]
  agg   <- aggregate_sample(samp)
  cells <- pop$cells

  res <- list(
    raw         = est_raw(samp),
    rake_margin = est_rake_margins(samp, pop_margins, pop$N),
    rake_joint  = est_rake_joint(samp, pop_margins, pop$N),
    mrp_int     = est_mrp(agg, cells, mrp_formulas$mrp_int),
    mrp_ia      = est_mrp(agg, cells, mrp_formulas$mrp_ia),
    mrp_ia_fe   = est_mrp(agg, cells, mrp_formulas$mrp_ia_fe),
    mrp_ia_fe_lin = est_mrp(agg, cells, mrp_formulas$mrp_ia_fe_lin),
    mrp_struct  = est_mrp(agg, cells, mrp_formulas$mrp_struct)
  )

  list(
    estimates = imap(res, \(r, nm) mutate(r$estimates, estimator = nm)) |> list_rbind(),
    diag      = imap(res, \(r, nm) mutate(r$diag, estimator = nm)) |> list_rbind(),
    sparsity  = sparsity_stats(samp, cells)
  )
}

# -----------------------------------------------------------------------------
# 5. Evaluation
# -----------------------------------------------------------------------------

estimator_labels <- c(
  raw         = "(1) Raw mean",
  rake_margin = "(2) Raking: margins",
  rake_joint  = "(3) Raking: region x age",
  mrp_int     = "(4) MRP: intercepts",
  mrp_ia      = "(5) MRP: + region x age",
  mrp_ia_fe   = "(6) MRP: (5), gender fixed",
  mrp_ia_fe_lin = "(6b) MRP: (6) + linear age trend",
  mrp_struct  = "(7) MRP: (6b) + region age slopes"
)

# Categorical slots of the validated reference palette, fixed per estimator
# (colour follows the estimator, not its position in a given plot)
estimator_cols <- c(
  raw = "#2a78d6", rake_margin = "#eb6834", rake_joint = "#1baf7a",
  mrp_int = "#eda100", mrp_ia = "#e87ba4", mrp_ia_fe = "#008300",
  mrp_ia_fe_lin = "#e34948", mrp_struct = "#6250d6"
)
label_cols <- function() set_names(estimator_cols[names(estimator_labels)], estimator_labels)

#' Bias, empirical SD and RMSE per estimator x estimand x sample size
summarise_performance <- function(estimates, truth) {
  estimates |>
    left_join(truth, by = "estimand") |>
    summarise(
      truth    = first(truth),
      n_reps   = n(),
      n_na     = sum(is.na(estimate)),
      mean_est = mean(estimate, na.rm = TRUE),
      bias     = mean(estimate - truth, na.rm = TRUE),
      emp_sd   = sd(estimate, na.rm = TRUE),
      rmse     = sqrt(mean((estimate - truth)^2, na.rm = TRUE)),
      .by = c(n, estimator, estimand)
    )
}

#' Estimator-level diagnostics per sample size
summarise_diagnostics <- function(diag) {
  diag |>
    summarise(
      n_reps            = n(),
      mean_kish_deff    = mean(deff, na.rm = TRUE),
      pct_failed        = 100 * mean(coalesce(failed, FALSE)),
      pct_nonconverged  = 100 * mean(!coalesce(converged, TRUE)),
      pct_collapsed     = 100 * mean(coalesce(empty_ra_collapsed, 0L) > 0),
      mean_empty_ra     = mean(empty_ra_collapsed, na.rm = TRUE),
      pct_singular      = 100 * mean(singular, na.rm = TRUE),
      pct_conv_warning  = 100 * mean(conv_warning, na.rm = TRUE),
      .by = c(n, estimator)
    ) |>
    mutate(
      across(c(mean_kish_deff, pct_nonconverged, pct_collapsed, mean_empty_ra),
             \(x) if_else(estimator %in% c("rake_margin", "rake_joint"), x, NA_real_)),
      across(c(pct_singular, pct_conv_warning),
             \(x) if_else(str_starts(estimator, "mrp"), x, NA_real_))
    )
}
