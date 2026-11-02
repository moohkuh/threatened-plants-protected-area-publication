# 03b-robust-cover-change.R
# Plot-level sensitivity analysis for cover change (Table S12).
# Species x plot means are analysed directly with robust LMMs
# instead of averaging them to one value per species and side.
# Needs the cache written by 03-analyse-cover-change.R.


# 0) Setup ---------------------------------------------------------------

source("./R-scripts/resurvey/00-preamble-resurvey.R")

suppressPackageStartupMessages({
  library(robustlmm)
  library(openxlsx)
})

RDS_FILE <- file.path(OUT_DIR_RDS, "species_cover_means_objects.rds")
assert_files_exist(RDS_FILE)

cache <- readRDS(RDS_FILE)

# The cached data only contain species with at least five plots
# inside and five plots outside (same species as in the main analysis).
scenario_data <- list(
  "All PAs"     = cache$dat_all_long,
  "Strict"      = cache$dat_strict_long,
  "IUCN IV"     = cache$dat_iucn4_long,
  "IUCN V"      = cache$dat_iucn5_long,
  "Natura 2000" = cache$dat_n2k_long
)


# 1) Fit one robust model per protected-area category -------------------

robust_results <- data.frame()
weight_summary <- data.frame()
robust_models <- list()

for (scenario in names(scenario_data)) {

  # average repeated intervals within species x plot
  plot_means <- scenario_data[[scenario]] %>%
    dplyr::group_by(RS_PROJECT, RS_PROJECT_PLOT, wcvp_name, side) %>%
    dplyr::summarise(
      rel_change_py_plot = mean(rel_change_py, na.rm = TRUE),
      .groups = "drop"
    ) %>%
    dplyr::filter(is.finite(rel_change_py_plot)) %>%
    dplyr::mutate(side = factor(side, levels = c("outside", "inside")))

  # random intercepts for species, project and plot
  model <- robustlmm::rlmer(
    rel_change_py_plot ~ side +
      (1 | wcvp_name) +
      (1 | RS_PROJECT) +
      (1 | RS_PROJECT_PLOT),
    data = plot_means
  )

  robust_models[[scenario]] <- model

  # "outside" is the reference level, so sideinside = inside - outside
  beta <- lme4::fixef(model)["sideinside"]

  # covariance of the fixed effects = sigma^2 x unscaled covariance.
  # Calculated directly because vcov() of lme4 >= 2.0 fails for rlmerMod.
  vc <- stats::sigma(model)^2 * model@pp$unsc()
  colnames(vc) <- names(lme4::fixef(model))
  rownames(vc) <- names(lme4::fixef(model))
  se <- sqrt(vc["sideinside", "sideinside"])

  # Wald CI and normal-approximation p-value
  z_value <- beta / se
  p_value <- 2 * stats::pnorm(abs(z_value), lower.tail = FALSE)

  robust_results <- rbind(
    robust_results,
    data.frame(
      PA_category = scenario,
      n_species = dplyr::n_distinct(plot_means$wcvp_name),
      n_species_plot_means = nrow(plot_means),
      n_plots = dplyr::n_distinct(plot_means$RS_PROJECT_PLOT),
      n_projects = dplyr::n_distinct(plot_means$RS_PROJECT),
      estimate = unname(beta),
      SE = se,
      lower_CI = unname(beta) - stats::qnorm(0.975) * se,
      upper_CI = unname(beta) + stats::qnorm(0.975) * se,
      z_value = unname(z_value),
      p_value = unname(p_value)
    )
  )

  # robustness weights: 1 = no downweighting, small = little influence
  w <- robustlmm::getME(model, "w_e")

  weight_summary <- rbind(
    weight_summary,
    data.frame(
      PA_category = scenario,
      median_weight = stats::median(w),
      n_weight_below_0.5 = sum(w < 0.5),
      prop_weight_below_0.5 = mean(w < 0.5)
    )
  )

  cat("\n", scenario, "\n")
  print(lme4::fixef(model))
  print(lme4::VarCorr(model))
}

print(robust_results)
print(weight_summary)


# 2) Diagnostic plots (DHARMa does not support rlmerMod) ----------------

for (scenario in names(robust_models)) {

  model <- robust_models[[scenario]]
  w <- robustlmm::getME(model, "w_e")

  png(
    file.path(OUT_DIAGNOSTIC_DIR, paste0("robust_cover_", gsub(" ", "_", scenario), ".png")),
    width = 1600, height = 1600, res = 200
  )
  par(mfrow = c(2, 2))

  plot(stats::fitted(model), stats::residuals(model),
       pch = 16, cex = 0.35,
       xlab = "Fitted values", ylab = "Residuals",
       main = paste(scenario, "- residuals vs fitted"))
  abline(h = 0, lty = 2)

  qqnorm(stats::residuals(model), pch = 16, cex = 0.35,
         main = paste(scenario, "- normal Q-Q"))
  qqline(stats::residuals(model), lty = 2)

  hist(w, breaks = 30, xlim = c(0, 1),
       xlab = "Robustness weight", main = paste(scenario, "- weights"))

  plot(stats::model.response(stats::model.frame(model)), w,
       pch = 16, cex = 0.35,
       xlab = "Mean annual relative cover change", ylab = "Robustness weight",
       main = paste(scenario, "- weights vs response"))
  abline(h = 1, lty = 2)

  dev.off()
}


# 3) Save results --------------------------------------------------------

readr::write_csv(
  robust_results,
  file.path(OUT_RESULTS_DIR, "robust_plotlevel_cover_change_results.csv")
)

readr::write_csv(
  weight_summary,
  file.path(OUT_DIAGNOSTIC_DIR, "robust_cover_weights.csv")
)

# formatted table as in Table S12
p_text <- formatC(robust_results$p_value, format = "f", digits = 3)
p_text[robust_results$p_value < 0.001] <- "<0.001"

table_s12 <- data.frame(
  "PA category" = robust_results$PA_category,
  "n species" = robust_results$n_species,
  "n species x plot means" = robust_results$n_species_plot_means,
  "n plots" = robust_results$n_plots,
  "Estimate" = formatC(robust_results$estimate, format = "f", digits = 4),
  "SE" = formatC(robust_results$SE, format = "f", digits = 4),
  "95% CI" = paste0(
    "[", formatC(robust_results$lower_CI, format = "f", digits = 4),
    ", ", formatC(robust_results$upper_CI, format = "f", digits = 4), "]"
  ),
  "p" = p_text,
  check.names = FALSE
)

wb <- openxlsx::createWorkbook()
openxlsx::addWorksheet(wb, "Table_S12")
openxlsx::writeData(wb, "Table_S12", table_s12)
openxlsx::addWorksheet(wb, "Raw_results")
openxlsx::writeData(wb, "Raw_results", robust_results)
openxlsx::saveWorkbook(
  wb,
  file.path(OUT_RESULTS_DIR, "robust_plotlevel_cover_change_results.xlsx"),
  overwrite = TRUE
)
