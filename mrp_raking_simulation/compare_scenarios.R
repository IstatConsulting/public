# =============================================================================
# MRP vs. raking simulation: comparison across selection scenarios
# -----------------------------------------------------------------------------
# Reads output/sel_ia_*/results_table.csv (written by run_simulation.R) and
# output/nonlin_sel_ia_*/ and produces a combined table, RMSE / bias figures
# faceted by estimand x strength of the region x age selection term (one per
# outcome-interaction shape), and a linear vs. non-linear comparison at k = 1.
# =============================================================================

source("R/sim_functions.R")

out_dir <- "output"

results <- list.files(out_dir, pattern = "^results_table\\.csv$", recursive = TRUE,
                      full.names = TRUE) |>
  set_names(\(p) basename(dirname(p))) |>
  map(\(p) read_csv(p, show_col_types = FALSE)) |>
  list_rbind(names_to = "scenario") |>
  mutate(
    shape  = if_else(str_starts(scenario, "nonlin_"), "nonlinear", "linear"),
    sel_ia = parse_number(str_remove(scenario, "^nonlin_"))
  )

write_csv(results, file.path(out_dir, "results_by_scenario.csv"))

# Compact RMSE grid (pp) for the focus estimands
rmse_grid <- results |>
  filter(estimand %in% c("National", "R2", "R5")) |>
  transmute(shape, sel_ia, estimand, n, estimator, rmse = round(100 * rmse, 2)) |>
  pivot_wider(names_from = estimator, values_from = rmse) |>
  arrange(shape, desc(sel_ia), estimand, n)
write_csv(rmse_grid, file.path(out_dir, "rmse_grid_by_scenario.csv"))

# Figures: raw mean omitted (its RMSE dominates the scale), see tables
est_cols <- c("#2a78d6", "#eb6834", "#1baf7a", "#eda100", "#e87ba4", "#008300", "#6250d6") |>
  set_names(estimator_labels)
facet_lab <- c(National = "National", focus_regions)

plot_scenarios <- function(df, metric, ylab, cols = "scenario") {
  df |>
    filter(estimand %in% names(facet_lab), estimator != "raw") |>
    mutate(
      estimand = factor(facet_lab[estimand], levels = facet_lab),
      scenario = factor(sprintf("Selection interaction x %.1f", sel_ia),
                        levels = sprintf("Selection interaction x %.1f", sort(unique(sel_ia), TRUE))),
      shape = factor(str_c(str_to_sentence(shape), " outcome interaction"),
                     levels = str_c(c("Linear", "Nonlinear"), " outcome interaction")),
      estimator_label = factor(estimator_labels[estimator], levels = estimator_labels),
      value = 100 * .data[[metric]]
    ) |>
    ggplot(aes(factor(n), value, colour = estimator_label, group = estimator_label)) +
    { if (metric == "bias") geom_hline(yintercept = 0, colour = "grey60", linewidth = 0.4) } +
    geom_line(linewidth = 0.6) +
    geom_point(size = 1.9) +
    facet_grid(reformulate(cols, "estimand"), scales = "free_y") +
    scale_colour_manual(values = est_cols, name = NULL, drop = TRUE) +
    labs(x = "Sample size (n)", y = ylab,
         caption = "Raw mean omitted (see rmse_grid_by_scenario.csv). (5) and (6) overlap almost exactly.") +
    theme_minimal(base_size = 11) +
    theme(
      legend.position = "bottom",
      panel.grid.minor = element_blank(),
      panel.grid.major.x = element_blank(),
      panel.grid.major.y = element_line(colour = "grey90", linewidth = 0.3),
      strip.text = element_text(face = "bold"),
      plot.caption = element_text(colour = "grey40")
    ) +
    guides(colour = guide_legend(nrow = 2, byrow = TRUE))
}

for (sh in unique(results$shape)) {
  d <- filter(results, shape == sh)
  ggsave(file.path(out_dir, sprintf("fig_rmse_by_scenario_%s.png", sh)),
         plot_scenarios(d, "rmse", "RMSE (percentage points)") +
           labs(title = sprintf("RMSE by estimator, sample size and selection interaction (%s outcome interaction)", sh)),
         width = 11, height = 8.5, dpi = 200, bg = "white")
  ggsave(file.path(out_dir, sprintf("fig_bias_by_scenario_%s.png", sh)),
         plot_scenarios(d, "bias", "Bias (percentage points)") +
           labs(title = sprintf("Bias by estimator, sample size and selection interaction (%s outcome interaction)", sh)),
         width = 11, height = 8.5, dpi = 200, bg = "white")
}

# Linear vs. non-linear outcome interaction at baseline selection (k = 1)
if (n_distinct(results$shape) > 1) {
  ggsave(file.path(out_dir, "fig_rmse_by_shape_k1.png"),
         plot_scenarios(filter(results, sel_ia == 1), "rmse", "RMSE (percentage points)",
                        cols = "shape") +
           labs(title = "RMSE by shape of the outcome interaction (baseline selection, k = 1)"),
         width = 9, height = 8.5, dpi = 200, bg = "white")
}

options(width = 200)
print(as.data.frame(rmse_grid))
