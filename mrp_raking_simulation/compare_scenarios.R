# =============================================================================
# MRP vs. raking simulation: comparison across selection scenarios
# -----------------------------------------------------------------------------
# Reads output/sel_ia_*/results_table.csv (written by run_simulation.R) and
# produces a combined table and RMSE / bias figures faceted by estimand x
# strength of the region x age selection term.
# =============================================================================

source("R/sim_functions.R")

out_dir <- "output"

results <- list.files(out_dir, pattern = "^results_table\\.csv$", recursive = TRUE,
                      full.names = TRUE) |>
  set_names(\(p) basename(dirname(p))) |>
  map(\(p) read_csv(p, show_col_types = FALSE)) |>
  list_rbind(names_to = "scenario") |>
  mutate(sel_ia = parse_number(scenario))

write_csv(results, file.path(out_dir, "results_by_scenario.csv"))

# Compact RMSE grid (pp) for the focus estimands
rmse_grid <- results |>
  filter(estimand %in% c("National", "R2", "R5")) |>
  transmute(sel_ia, estimand, n, estimator, rmse = round(100 * rmse, 2)) |>
  pivot_wider(names_from = estimator, values_from = rmse) |>
  arrange(desc(sel_ia), estimand, n)
write_csv(rmse_grid, file.path(out_dir, "rmse_grid_by_scenario.csv"))

# Figures: raw mean omitted (its RMSE dominates the scale), see tables
est_cols <- c("#2a78d6", "#eb6834", "#1baf7a", "#eda100", "#e87ba4", "#008300", "#6250d6") |>
  set_names(estimator_labels)
facet_lab <- c(National = "National", focus_regions)

plot_scenarios <- function(metric, ylab) {
  results |>
    filter(estimand %in% names(facet_lab), estimator != "raw") |>
    mutate(
      estimand = factor(facet_lab[estimand], levels = facet_lab),
      scenario = factor(sprintf("Selection interaction x %.1f", sel_ia),
                        levels = sprintf("Selection interaction x %.1f", sort(unique(sel_ia), TRUE))),
      estimator_label = factor(estimator_labels[estimator], levels = estimator_labels),
      value = 100 * .data[[metric]]
    ) |>
    ggplot(aes(factor(n), value, colour = estimator_label, group = estimator_label)) +
    { if (metric == "bias") geom_hline(yintercept = 0, colour = "grey60", linewidth = 0.4) } +
    geom_line(linewidth = 0.6) +
    geom_point(size = 1.9) +
    facet_grid(estimand ~ scenario, scales = "free_y") +
    scale_colour_manual(values = est_cols, name = NULL, drop = TRUE) +
    labs(x = "Sample size (n)", y = ylab,
         caption = "Raw mean omitted from the figure; see rmse_grid_by_scenario.csv.") +
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

ggsave(file.path(out_dir, "fig_rmse_by_scenario.png"),
       plot_scenarios("rmse", "RMSE (percentage points)") +
         labs(title = "RMSE by estimator, sample size and strength of the selection interaction"),
       width = 11, height = 8.5, dpi = 200, bg = "white")
ggsave(file.path(out_dir, "fig_bias_by_scenario.png"),
       plot_scenarios("bias", "Bias (percentage points)") +
         labs(title = "Bias by estimator, sample size and strength of the selection interaction"),
       width = 11, height = 8.5, dpi = 200, bg = "white")

options(width = 200)
print(as.data.frame(rmse_grid))
