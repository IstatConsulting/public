# =============================================================================
# MRP vs. raking simulation: driver
# -----------------------------------------------------------------------------
# Usage (from this folder):  Rscript run_simulation.R
# Optional env vars: SIM_R (replications, default 500), SIM_WORKERS,
#                    SIM_N_STAN (rstanarm subset per condition, 0 = skip),
#                    SIM_SEL_IA (scaling of the region x age selection term,
#                    default 1; results go to output/sel_ia_<value>/)
# =============================================================================

source("R/sim_functions.R")
suppressPackageStartupMessages({
  library(furrr)
})

# -----------------------------------------------------------------------------
# Settings
# -----------------------------------------------------------------------------
seed         <- 20261018
N_pop        <- 50000
sample_sizes <- c(200, 500, 2000)
R            <- as.integer(Sys.getenv("SIM_R", 500))
n_stan       <- as.integer(Sys.getenv("SIM_N_STAN", 20))
workers      <- as.integer(Sys.getenv("SIM_WORKERS", parallel::detectCores()))
sel_ia       <- as.numeric(Sys.getenv("SIM_SEL_IA", 1))
out_dir      <- file.path("output", sprintf("sel_ia_%.1f", sel_ia))
dir.create(out_dir, showWarnings = FALSE, recursive = TRUE)

# -----------------------------------------------------------------------------
# Population and samples (all randomness fixed here, before estimation)
# -----------------------------------------------------------------------------
set.seed(seed)
pop         <- simulate_population(N = N_pop, sp = scale_selection_ia(k = sel_ia))
pop_margins <- make_pop_margins(pop$cells)

samples <- map(sample_sizes, \(n) draw_sample_ids(pop, n, R)) |>
  set_names(sample_sizes)

pop$cells |>
  select(cell, region, age, gender, educ, N_c, y_c, p_true, log_incl) |>
  write_csv(file.path(out_dir, "population_cells.csv"))
write_csv(pop$truth, file.path(out_dir, "truth.csv"))

# -----------------------------------------------------------------------------
# Estimation: the same sample goes to all five estimators
# -----------------------------------------------------------------------------
plan(multisession, workers = workers)
sim_pkgs <- c("tidyverse", "survey", "lme4", "rstanarm")

reps <- imap(samples, \(id_list, n) {
  message("n = ", n, ": ", length(id_list), " replications")
  future_map(id_list, run_one_rep, pop = pop, pop_margins = pop_margins,
             .options = furrr_options(seed = TRUE, packages = sim_pkgs)) |>
    imap(\(x, r) map(x, \(tbl) mutate(tbl, n = as.integer(n), rep = r)))
}) |>
  list_flatten()

estimates <- map(reps, "estimates") |> list_rbind()
diag      <- map(reps, "diag") |> list_rbind()
sparsity  <- map(reps, "sparsity") |> list_rbind()

write_csv(estimates, file.path(out_dir, "replication_estimates.csv.gz"))
write_csv(diag, file.path(out_dir, "replication_diagnostics.csv"))
write_csv(sparsity, file.path(out_dir, "replication_sparsity.csv"))

# -----------------------------------------------------------------------------
# Evaluation
# -----------------------------------------------------------------------------
perf <- summarise_performance(estimates, pop$truth) |>
  mutate(estimator_label = estimator_labels[estimator]) |>
  arrange(estimand, n, estimator)

diag_sum <- summarise_diagnostics(diag) |>
  mutate(estimator_label = estimator_labels[estimator]) |>
  arrange(n, estimator)

sparsity_sum <- sparsity |>
  summarise(
    mean_empty_cells      = mean(empty_cells),
    mean_cells_1_2        = mean(cells_1_2),
    mean_empty_region_age = mean(empty_region_age),
    pct_any_empty_ra      = 100 * mean(empty_region_age > 0),
    pct_any_empty_region  = 100 * mean(empty_region > 0),
    .by = n
  )

# Tidy results table: performance joined with estimator diagnostics
results <- perf |>
  left_join(select(diag_sum, -n_reps, -estimator_label), by = c("n", "estimator"))

write_csv(results, file.path(out_dir, "results_table.csv"))
write_csv(diag_sum, file.path(out_dir, "results_diagnostics.csv"))
write_csv(sparsity_sum, file.path(out_dir, "results_sparsity.csv"))

# -----------------------------------------------------------------------------
# Plots
# -----------------------------------------------------------------------------
# Categorical slots 1-5 of the validated reference palette, fixed order
est_cols <- c("#2a78d6", "#eb6834", "#1baf7a", "#eda100", "#e87ba4", "#008300", "#6250d6") |>
  set_names(estimator_labels)

facet_lab <- c(National = "National", focus_regions)

plot_metric <- function(df, metric, ylab) {
  df |>
    filter(estimand %in% names(facet_lab)) |>
    mutate(
      estimand = factor(facet_lab[estimand], levels = facet_lab),
      estimator_label = factor(estimator_label, levels = estimator_labels),
      value = 100 * .data[[metric]]
    ) |>
    ggplot(aes(factor(n), value, colour = estimator_label, group = estimator_label)) +
    { if (metric == "bias") geom_hline(yintercept = 0, colour = "grey60", linewidth = 0.4) } +
    geom_line(linewidth = 0.7) +
    geom_point(size = 2.4) +
    facet_wrap(~estimand, scales = "free_y") +
    scale_colour_manual(values = est_cols, name = NULL) +
    labs(x = "Sample size (n)", y = ylab) +
    theme_minimal(base_size = 11) +
    theme(
      legend.position = "bottom",
      panel.grid.minor = element_blank(),
      panel.grid.major.x = element_blank(),
      panel.grid.major.y = element_line(colour = "grey90", linewidth = 0.3),
      strip.text = element_text(face = "bold", hjust = 0)
    ) +
    guides(colour = guide_legend(nrow = 2, byrow = TRUE))
}

p_rmse <- plot_metric(perf, "rmse", "RMSE (percentage points)") +
  labs(title = "RMSE by estimator and sample size",
       subtitle = sprintf("%d replications per condition; region x age selection term scaled by %.1f",
                          R, sel_ia))
ggsave(file.path(out_dir, "fig_rmse.png"), p_rmse, width = 10.5, height = 5.2, dpi = 200, bg = "white")

p_bias <- plot_metric(perf, "bias", "Bias (percentage points)") +
  labs(title = "Bias by estimator and sample size")
ggsave(file.path(out_dir, "fig_bias.png"), p_bias, width = 10.5, height = 5.2, dpi = 200, bg = "white")

# -----------------------------------------------------------------------------
# Optional robustness check: rstanarm on a subset of the same samples
# -----------------------------------------------------------------------------
if (n_stan > 0) {
  message("rstanarm robustness check: ", n_stan, " replications per condition")
  stan_jobs <- expand_grid(n = sample_sizes, rep = seq_len(min(n_stan, R)),
                           model = c("mrp_int", "mrp_ia"))
  stan_out <- future_pmap(stan_jobs, \(n, rep, model) {
    ids <- samples[[as.character(n)]][[rep]]
    agg <- aggregate_sample(pop$persons[ids, ])
    s <- est_mrp_stan(agg, pop$cells, mrp_formulas[[model]], seed = seed + rep)
    s$estimates |> mutate(n = n, rep = rep, estimator = str_c(model, "_stan")) |>
      bind_cols(s$diag)
  }, .options = furrr_options(seed = TRUE, packages = sim_pkgs)) |>
    list_rbind()
  write_csv(stan_out, file.path(out_dir, "replication_estimates_stan.csv"))

  # Compare glmer and stan on the identical subset of replications
  stan_check <- bind_rows(
    stan_out |> select(n, rep, estimator, estimand, estimate),
    estimates |>
      filter(estimator %in% c("mrp_int", "mrp_ia")) |>
      semi_join(distinct(stan_out, n, rep), by = c("n", "rep"))
  ) |>
    summarise_performance(pop$truth) |>
    filter(estimand %in% names(facet_lab)) |>
    left_join(
      stan_out |>
        distinct(n, rep, estimator, divergent, low_ess, max_rhat) |>
        summarise(pct_divergent_fits = 100 * mean(divergent > 0),
                  pct_low_ess = 100 * mean(low_ess),
                  max_rhat = max(max_rhat), .by = c(n, estimator)),
      by = c("n", "estimator")
    ) |>
    arrange(estimand, n, estimator)
  write_csv(stan_check, file.path(out_dir, "results_stan_check.csv"))
}

plan(sequential)

# -----------------------------------------------------------------------------
# Expected-pattern check (printed and saved)
# -----------------------------------------------------------------------------
rmse_wide <- perf |>
  select(n, estimand, estimator, rmse) |>
  pivot_wider(names_from = estimator, values_from = rmse)

check_lines <- perf |>
  filter(estimand %in% names(facet_lab)) |>
  summarise(
    best = estimator[which.min(rmse)],
    `(3) <= min MRP` = rmse[estimator == "rake_joint"] <=
      min(rmse[estimator %in% c("mrp_int", "mrp_ia")]),
    `(5) within 10% of best` = rmse[estimator == "mrp_ia"] <= 1.1 * min(rmse),
    `(7) within 10% of best` = rmse[estimator == "mrp_struct"] <= 1.1 * min(rmse),
    .by = c(n, estimand)
  )

sink(file.path(out_dir, "pattern_check.txt"))
cat("RMSE (pp) by condition\n")
print(rmse_wide |> mutate(across(raw:mrp_ia, \(x) round(100 * x, 2))), n = Inf, width = Inf)
cat("\nExpected-pattern checks\n")
print(check_lines, n = Inf, width = Inf)
cat("\nSparsity\n"); print(sparsity_sum, width = Inf)
cat("\nDiagnostics\n"); print(select(diag_sum, -estimator_label), n = Inf, width = Inf)
sink()
cat(readLines(file.path(out_dir, "pattern_check.txt")), sep = "\n")
