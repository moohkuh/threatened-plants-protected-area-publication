# 02-Fit-plotlevel-event-models.R
#
# Plot-level colonization and final-extinction analysis
#
# Main analysis:
#   Excludes Red List category V.
#
# Sensitivity analysis:
#   Repeats the analysis including Red List category V.
#
# Event definitions:
#   Colonization = colonization not followed by a later extinction.
#   Final extinction = extinction not followed by later recolonization.
#
# Colonization model:
#   Negative-binomial GLMM
#   Effect measure = incidence rate ratio (IRR)
#
# Final-extinction model:
#   Binomial GLMM
#   Effect measure = odds ratio (OR)
#
# Final-extinction risk set:
#   Unique threatened species that were present at the beginning
#   of at least one observed survey interval.
#
# Each species contributes at most once to the numerator and
# once to the denominator of the final-extinction analysis.

# 0) Setup ----------------------------------------------------------------

source("./R-scripts/resurvey/00-preamble-resurvey.R")

suppressPackageStartupMessages({
  library(dplyr)
  library(purrr)
  library(ggplot2)
  library(glmmTMB)
  library(DHARMa)
  library(patchwork)
})

# ---------------------- Output directories ---------------------------------

OUT_EVENT_RESULTS_DIR <- file.path(OUT_RESULTS_DIR, "events")
OUT_EVENT_DIAGNOSTIC_DIR <- file.path(OUT_DIAGNOSTIC_DIR, "events")
OUT_EVENT_FIG_DIR <- OUT_FIG_RS_DIR
OUT_EVENT_MODELS_DIR <- file.path(OUT_DIR_MODELS, "events")

for (directory in c(OUT_EVENT_RESULTS_DIR, OUT_EVENT_DIAGNOSTIC_DIR, OUT_EVENT_FIG_DIR, OUT_EVENT_MODELS_DIR)) {
  dir.create(directory, recursive = TRUE, showWarnings = FALSE)
}

assert_files_exist(c(OUT_RS_CLEAN))


# 1) Plotting settings -----------------------------------------

scenario_levels <- c("All PAs", "Strict", "IUCN IV", "IUCN V", "Natura 2000")

scenario_cols <- c(
  "All PAs" = "#5D7388FF",
  "Strict" = "#8AAE8AFF",
  "IUCN IV" = "#4F4D1DFF",
  "IUCN V" = "#A08F27FF",
  "Natura 2000" = "#E5A729FF"
)


# 2) Load clean ReSurvey data ----------------------------------

RS_clean_flagged <- readr::read_csv(OUT_RS_CLEAN, show_col_types = FALSE)

assert_columns(
  RS_clean_flagged,
  c(
    "RS_PROJECT", "RS_PLOT", "RS_PROJECT_PLOT", "wcvp_name", "threat_2018", "from", "to",
    "resurvey_years", "SURF_AREA", "LONGITUDE", "LATITUDE", "absolute.change.colonizer",
    "absolute.change.extinct", "any_PA", "any_strict", "any_n2k", "any_iucn4", "any_iucn5", "is_nat_rl"
  ),
  "clean ReSurvey data"
)

# Keep nationally Red-Listed species only.
rs_base <- RS_clean_flagged %>%
  dplyr::filter(is_nat_rl == TRUE) %>%
  dplyr::mutate(threat_2018 = stringr::str_trim(as.character(threat_2018)))


# 3) Collapse duplicate rows caused by overlapping PA polygons -------------------

# Build exactly one row per plot x species x survey interval.
# PA membership is summarized using binary flags.

obs_raw <- rs_base %>%
  dplyr::group_by(RS_PROJECT, RS_PLOT, RS_PROJECT_PLOT, wcvp_name, threat_2018, from, to) %>%
  dplyr::summarise(
    years = dplyr::first(resurvey_years),
    absolute.change.colonizer = dplyr::first(absolute.change.colonizer),
    absolute.change.extinct = dplyr::first(absolute.change.extinct),
    SURF_AREA = dplyr::first(SURF_AREA),
    any_PA = dplyr::first(any_PA),
    any_strict = dplyr::first(any_strict),
    any_n2k = dplyr::first(any_n2k),
    any_iucn4 = dplyr::first(any_iucn4),
    any_iucn5 = dplyr::first(any_iucn5),
    .groups = "drop"
  ) %>%
  dplyr::filter(!is.na(years), years > 0)


# 4) Define retained colonization and extinction events -------------------

# Raw colonizations and extinctions come from the ReSurvey
# interval-change table.
#
# For plots with multiple survey intervals:
#
# - a colonization is retained only if it is not followed by
#   a later extinction;
#
# - an extinction is retained only if it is not followed by
#   a later recolonization.

obs_raw <- obs_raw %>%
  dplyr::mutate(
    is_col = !is.na(absolute.change.colonizer),
    is_ext = !is.na(absolute.change.extinct),
    # A species is present at interval start unless the
    # interval itself represents a colonization.
    present_at_start = !is_col
  ) %>%
  dplyr::arrange(RS_PROJECT_PLOT, wcvp_name, from, to) %>%
  dplyr::group_by(RS_PROJECT_PLOT, wcvp_name) %>%
  dplyr::mutate(
    # Position of the last extinction event.
    last_ext_row = ifelse(any(is_ext), max(which(is_ext)), NA_integer_),
    # TRUE for observations occurring before the last extinction.
    any_future_ext = dplyr::row_number() < last_ext_row,
    # Retain colonizations that persist until the end
    # of the observed time series.
    is_col_strict = dplyr::case_when(!is_col ~ FALSE, any_future_ext ~ FALSE, TRUE ~ TRUE),
    # Position of the last colonization event.
    last_col_row = ifelse(any(is_col), max(which(is_col)), NA_integer_),
    # TRUE for observations occurring before the last colonization.
    any_future_col = dplyr::row_number() < last_col_row,
    # Retain extinctions that are not followed by recolonization.
    is_ext_strict = dplyr::case_when(!is_ext ~ FALSE, any_future_col ~ FALSE, TRUE ~ TRUE)
  ) %>%
  dplyr::ungroup() %>%
  dplyr::select(-last_ext_row, -last_col_row)

cat("\nRetained colonizations:", sum(obs_raw$is_col_strict, na.rm = TRUE), "\n")
cat("Retained final extinctions:", sum(obs_raw$is_ext_strict, na.rm = TRUE), "\n")


# 5) Use frozen protected-area classification ---------------------------

# Protected-area status is already attached in 01-build-resurvey-data.R
# from the frozen plot-level classification derived from exact coordinates.


# 6) Create one PA classification per plot -----------------------------

# PA membership is a plot-level property.
# It is therefore collapsed once and later joined to the
# colonization and extinction datasets.

plot_info <- obs_raw %>%
  dplyr::group_by(RS_PROJECT, RS_PLOT, RS_PROJECT_PLOT) %>%
  dplyr::summarise(
    SURF_AREA = dplyr::first(SURF_AREA),
    any_PA = any(any_PA, na.rm = TRUE),
    any_strict = any(any_strict, na.rm = TRUE),
    any_n2k = any(any_n2k, na.rm = TRUE),
    any_iucn4 = any(any_iucn4, na.rm = TRUE),
    any_iucn5 = any(any_iucn5, na.rm = TRUE),
    .groups = "drop"
  )

stopifnot(!anyDuplicated(plot_info$RS_PROJECT_PLOT))


# 7) Extract one coordinate pair per plot -------------------------

coords_per_plot <- RS_clean_flagged %>%
  dplyr::filter(is.finite(LONGITUDE), is.finite(LATITUDE)) %>%
  dplyr::group_by(RS_PROJECT_PLOT) %>%
  dplyr::summarise(
    LONGITUDE = median(LONGITUDE, na.rm = TRUE),
    LATITUDE = median(LATITUDE, na.rm = TRUE),
    .groups = "drop"
  )


# 8) Helper for category-specific inside-vs-outside comparisons ---------------
# For category-specific comparisons:
#
# inside  = plot lies inside the focal PA category
# outside = plot lies outside the complete PA network
#
# Plots inside another PA category only are excluded.

make_scenario <- function(data, inside_flag) {
  data %>%
    dplyr::mutate(
      side = dplyr::case_when(.data[[inside_flag]] ~ "inside", !any_PA ~ "outside", TRUE ~ NA_character_),
      side = factor(side, levels = c("outside", "inside"))
    ) %>%
    dplyr::filter(!is.na(side))
}


# 9) Main analysis function -------------------------------------

run_event_analysis <- function(include_v, run_dharma = FALSE) {

  suffix <- if (include_v) { "includingV" } else { "withoutV" }

  # ----------- Select species set --------------------------

  obs_analysis <- obs_raw

  if (!include_v) {
    obs_analysis <- obs_analysis %>% dplyr::filter(threat_2018 != "V")
  }

  # 10) COLONIZATION ------------------------------------------------------

  # ----------- Aggregate events to plot x survey interval ----------------

  plot_interval <- obs_analysis %>%
    dplyr::group_by(RS_PROJECT, RS_PLOT, RS_PROJECT_PLOT, from, to) %>%
    dplyr::summarise(
      years = dplyr::first(years),
      n_species = dplyr::n_distinct(wcvp_name),
      n_col = sum(is_col_strict, na.rm = TRUE),
      .groups = "drop"
    )

  # --------------- Aggregate colonization events to one row per plot ---------------

  plot_level_col <- plot_interval %>%
    dplyr::group_by(RS_PROJECT, RS_PLOT, RS_PROJECT_PLOT) %>%
    dplyr::summarise(
      n_intervals = dplyr::n(),
      total_years = sum(years, na.rm = TRUE),
      total_col = sum(n_col, na.rm = TRUE),
      mean_n_species = mean(n_species, na.rm = TRUE),
      .groups = "drop"
    ) %>%
    dplyr::left_join(plot_info, by = c("RS_PROJECT", "RS_PLOT", "RS_PROJECT_PLOT")) %>%
    dplyr::filter(!is.na(SURF_AREA), SURF_AREA > 0, total_years > 0) %>%
    dplyr::mutate(
      log_area = log(SURF_AREA),
      log_area_z = as.numeric(scale(log_area, center = TRUE, scale = TRUE)),
      total_years_c = as.numeric(scale(total_years, center = TRUE, scale = FALSE))
    ) %>%
    dplyr::left_join(coords_per_plot, by = "RS_PROJECT_PLOT")

  # -------------------- Colonization scenarios ------------------

  pl_all <- plot_level_col %>%
    dplyr::mutate(
      side = dplyr::if_else(any_PA, "inside", "outside"),
      side = factor(side, levels = c("outside", "inside"))
    )

  pl_strict <- make_scenario(plot_level_col, "any_strict")
  pl_iucn4 <- make_scenario(plot_level_col, "any_iucn4")
  pl_iucn5 <- make_scenario(plot_level_col, "any_iucn5")
  pl_n2k <- make_scenario(plot_level_col, "any_n2k")

  # --------- Negative-binomial colonization models ---------------------------

  fit_col_model <- function(data) {
    glmmTMB::glmmTMB(
      total_col ~ side + log_area_z + total_years_c + (1 | RS_PROJECT),
      data = data,
      family = glmmTMB::nbinom2()
    )
  }

  m_all_col <- fit_col_model(pl_all)
  m_strict_col <- fit_col_model(pl_strict)
  m_iucn4_col <- fit_col_model(pl_iucn4)
  m_iucn5_col <- fit_col_model(pl_iucn5)
  m_n2k_col <- fit_col_model(pl_n2k)

  # 11) FINAL EXTINCTION -----------------------------------------------------------

  # ------- Reduce each plot x species time series to one observation --------------

  final_loss_species <- obs_analysis %>%
    dplyr::group_by(RS_PROJECT, RS_PLOT, RS_PROJECT_PLOT, wcvp_name) %>%
    dplyr::summarise(
      # Species was present at the start of at least one interval.
      ever_at_risk = any(present_at_start, na.rm = TRUE),
      # Species underwent a retained final extinction.
      final_extinction = any(is_ext_strict, na.rm = TRUE),
      .groups = "drop"
    )

  bad_final_extinction <- final_loss_species %>%
    dplyr::filter(final_extinction, !ever_at_risk)

  # Keep only species that had at least one observed opportunity
  # to go locally extinct.
  final_loss_species <- final_loss_species %>%
    dplyr::filter(ever_at_risk)

  # ------------- Aggregate final-extinction risk set to one row per plot ------------------------

  plot_level_ext <- final_loss_species %>%
    dplyr::group_by(RS_PROJECT, RS_PLOT, RS_PROJECT_PLOT) %>%
    dplyr::summarise(
      n_at_risk_final = dplyr::n(),
      n_final_ext = sum(final_extinction, na.rm = TRUE),
      .groups = "drop"
    )

  # -------------- Summarize survey metadata independently ------------------

  plot_metadata_ext <- obs_analysis %>%
    dplyr::distinct(RS_PROJECT, RS_PLOT, RS_PROJECT_PLOT, from, to, years) %>%
    dplyr::group_by(RS_PROJECT, RS_PLOT, RS_PROJECT_PLOT) %>%
    dplyr::summarise(
      n_intervals = dplyr::n(),
      total_years = sum(years, na.rm = TRUE),
      .groups = "drop"
    )

  # --------------- Add plot information and covariates -------------------------------

  plot_level_ext <- plot_level_ext %>%
    dplyr::left_join(plot_info, by = c("RS_PROJECT", "RS_PLOT", "RS_PROJECT_PLOT")) %>%
    dplyr::left_join(plot_metadata_ext, by = c("RS_PROJECT", "RS_PLOT", "RS_PROJECT_PLOT")) %>%
    dplyr::filter(n_at_risk_final > 0, !is.na(SURF_AREA), SURF_AREA > 0, total_years > 0) %>%
    dplyr::mutate(
      log_area = log(SURF_AREA),
      log_area_z = as.numeric(scale(log_area, center = TRUE, scale = TRUE)),
      total_years_c = as.numeric(scale(total_years, center = TRUE, scale = FALSE))
    ) %>%
    dplyr::left_join(coords_per_plot, by = "RS_PROJECT_PLOT")

  # Final-extinction scenarios --------------------------------------------------

  ext_all <- plot_level_ext %>%
    dplyr::mutate(
      side = dplyr::if_else(any_PA, "inside", "outside"),
      side = factor(side, levels = c("outside", "inside"))
    )

  ext_strict <- make_scenario(plot_level_ext, "any_strict")
  ext_iucn4 <- make_scenario(plot_level_ext, "any_iucn4")
  ext_iucn5 <- make_scenario(plot_level_ext, "any_iucn5")
  ext_n2k <- make_scenario(plot_level_ext, "any_n2k")

  # Binomial final-extinction models ---------------------------

  fit_ext_model <- function(data) {
    glmmTMB::glmmTMB(
      cbind(n_final_ext, n_at_risk_final - n_final_ext) ~ side + log_area_z + total_years_c + (1 | RS_PROJECT),
      data = data,
      family = binomial()
    )
  }

  m_all_ext <- fit_ext_model(ext_all)
  m_strict_ext <- fit_ext_model(ext_strict)
  m_iucn4_ext <- fit_ext_model(ext_iucn4)
  m_iucn5_ext <- fit_ext_model(ext_iucn5)
  m_n2k_ext <- fit_ext_model(ext_n2k)

  # 12) Store models ---------------------------------------------

  models_final <- list(
    all_col = m_all_col, all_ext = m_all_ext,
    strict_col = m_strict_col, strict_ext = m_strict_ext,
    iucn4_col = m_iucn4_col, iucn4_ext = m_iucn4_ext,
    iucn5_col = m_iucn5_col, iucn5_ext = m_iucn5_ext,
    n2k_col = m_n2k_col, n2k_ext = m_n2k_ext
  )

  # 13) Model validity ----------------------------------------

  model_validity <- purrr::imap_dfr(
    models_final,
    ~ tibble::tibble(
      model = .y,
      convergence_code = .x$fit$convergence,
      pdHess = isTRUE(.x$sdr$pdHess),
      finite_se = tryCatch(
        all(is.finite(summary(.x)$coefficients$cond[, "Std. Error"])),
        error = function(e) FALSE
      )
    )
  ) %>%
    dplyr::mutate(valid = convergence_code == 0 & pdHess & finite_se)

  print(model_validity)

  readr::write_csv(model_validity, file.path(OUT_EVENT_RESULTS_DIR, paste0("model_validity_", suffix, ".csv")))

  # 14) Extract model effects --------------------------------

  get_effect <- function(model, scenario_label, outcome_label, effect_type) {

    coefs <- summary(model)$coefficients$cond

    beta <- coefs["sideinside", "Estimate"]
    se <- coefs["sideinside", "Std. Error"]
    p <- coefs["sideinside", "Pr(>|z|)"]

    tibble::tibble(
      scenario = scenario_label,
      outcome = outcome_label,
      effect_type = effect_type,
      estimate = exp(beta),
      lcl = exp(beta - 1.96 * se),
      ucl = exp(beta + 1.96 * se),
      p_val = p
    )
  }

  effect_tbl <- dplyr::bind_rows(
    get_effect(m_all_col, "All PAs", "Colonization", "IRR"),
    get_effect(m_all_ext, "All PAs", "Final extinction", "OR"),
    get_effect(m_strict_col, "Strict", "Colonization", "IRR"),
    get_effect(m_strict_ext, "Strict", "Final extinction", "OR"),
    get_effect(m_iucn4_col, "IUCN IV", "Colonization", "IRR"),
    get_effect(m_iucn4_ext, "IUCN IV", "Final extinction", "OR"),
    get_effect(m_iucn5_col, "IUCN V", "Colonization", "IRR"),
    get_effect(m_iucn5_ext, "IUCN V", "Final extinction", "OR"),
    get_effect(m_n2k_col, "Natura 2000", "Colonization", "IRR"),
    get_effect(m_n2k_ext, "Natura 2000", "Final extinction", "OR")
  ) %>%
    dplyr::mutate(
      sig_flag = dplyr::if_else(p_val < 0.05, "significant", "not significant"),
      p_fmt = format.pval(p_val, digits = 2, eps = 0.001),
      scenario = factor(scenario, levels = scenario_levels)
    )

  cat("\nModel effects:\n")

  print(effect_tbl %>% dplyr::mutate(dplyr::across(c(estimate, lcl, ucl, p_val), ~ round(.x, 4))))

  readr::write_csv(effect_tbl, file.path(OUT_EVENT_RESULTS_DIR, paste0("event_effects_", suffix, ".csv")))

  # 15) Supplement table: colonization ---------------------------

  supplement_col <- purrr::imap_dfr(
    list(
      "All PAs" = pl_all, "Strict" = pl_strict, "IUCN IV" = pl_iucn4,
      "IUCN V" = pl_iucn5, "Natura 2000" = pl_n2k
    ),
    ~ .x %>%
      dplyr::group_by(side) %>%
      dplyr::summarise(
        n_plots = dplyr::n(),
        total_intervals = sum(n_intervals, na.rm = TRUE),
        median_intervals = median(n_intervals, na.rm = TRUE),
        sum_col = sum(total_col, na.rm = TRUE),
        mean_col = round(mean(total_col, na.rm = TRUE), 2),
        median_years = round(median(total_years, na.rm = TRUE), 1),
        .groups = "drop"
      ) %>%
      dplyr::mutate(scenario = .y)
  )

  readr::write_csv(
    supplement_col,
    file.path(OUT_EVENT_RESULTS_DIR, paste0("supplement_colonization_", suffix, ".csv"))
  )

  # 16) Supplement table: final extinction ------------------------

  supplement_ext <- purrr::imap_dfr(
    list(
      "All PAs" = ext_all, "Strict" = ext_strict, "IUCN IV" = ext_iucn4,
      "IUCN V" = ext_iucn5, "Natura 2000" = ext_n2k
    ),
    ~ .x %>%
      dplyr::group_by(side) %>%
      dplyr::summarise(
        n_plots = dplyr::n(),
        total_intervals = sum(n_intervals, na.rm = TRUE),
        median_intervals = median(n_intervals, na.rm = TRUE),
        sum_final_ext = sum(n_final_ext, na.rm = TRUE),
        mean_final_ext = round(mean(n_final_ext, na.rm = TRUE), 2),
        sum_at_risk_final = sum(n_at_risk_final, na.rm = TRUE),
        mean_at_risk_final = round(mean(n_at_risk_final, na.rm = TRUE), 2),
        median_years = round(median(total_years, na.rm = TRUE), 1),
        .groups = "drop"
      ) %>%
      dplyr::mutate(scenario = .y)
  )

  readr::write_csv(
    supplement_ext,
    file.path(OUT_EVENT_RESULTS_DIR, paste0("supplement_final_extinction_", suffix, ".csv"))
  )

  # 17) DHARMa diagnostics ------------------------------------

  dharma_final <- NULL

  if (run_dharma) {

    diag_dir <- file.path(OUT_EVENT_DIAGNOSTIC_DIR, suffix)

    dir.create(diag_dir, recursive = TRUE, showWarnings = FALSE)

    dharma_final <- purrr::imap_dfr(models_final, ~ {

      sim <- DHARMa::simulateResiduals(.x, n = 1000, seed = 42)

      grDevices::png(
        filename = file.path(diag_dir, paste0(.y, "_DHARMa.png")),
        width = 1400, height = 700, res = 130
      )

      plot(sim, main = .y)

      grDevices::dev.off()

      tibble::tibble(
        model = .y,
        uniformity_p = tryCatch(DHARMa::testUniformity(sim, plot = FALSE)$p.value, error = function(e) NA_real_),
        dispersion_p = tryCatch(DHARMa::testDispersion(sim, plot = FALSE)$p.value, error = function(e) NA_real_),
        zero_p = tryCatch(DHARMa::testZeroInflation(sim, plot = FALSE)$p.value, error = function(e) NA_real_),
        outlier_p = tryCatch(
          suppressWarnings(DHARMa::testOutliers(sim, plot = FALSE)$p.value),
          error = function(e) NA_real_
        )
      )
    })

    print(dharma_final)

    readr::write_csv(dharma_final, file.path(diag_dir, "DHARMa_summary.csv"))
  }

  # 18) Save model objects and plot-level datasets ---------------

  purrr::iwalk(
    models_final,
    ~ saveRDS(.x, file.path(OUT_EVENT_MODELS_DIR, paste0("plotlevel_", suffix, "_", .y, ".rds")))
  )

  readr::write_csv(
    plot_level_col,
    file.path(OUT_EVENT_RESULTS_DIR, paste0("plot_level_colonization_", suffix, ".csv"))
  )

  readr::write_csv(
    plot_level_ext,
    file.path(OUT_EVENT_RESULTS_DIR, paste0("plot_level_final_extinction_", suffix, ".csv"))
  )

  # 19) Return analysis objects ----------------------

  list(
    plot_interval = plot_interval,
    plot_level_col = plot_level_col,

    final_loss_species = final_loss_species,
    plot_level_ext = plot_level_ext,

    pl_all = pl_all, pl_strict = pl_strict, pl_iucn4 = pl_iucn4, pl_iucn5 = pl_iucn5, pl_n2k = pl_n2k,
    ext_all = ext_all, ext_strict = ext_strict, ext_iucn4 = ext_iucn4, ext_iucn5 = ext_iucn5, ext_n2k = ext_n2k,

    models = models_final,
    effects = effect_tbl,
    diagnostics = dharma_final,

    supplement_col = supplement_col,
    supplement_ext = supplement_ext
  )
}


# 20) Run main analysis and including-V sensitivity ---------------------

main_without_v <- run_event_analysis(include_v = FALSE, run_dharma = FALSE)

sensitivity_with_v <- run_event_analysis(include_v = TRUE, run_dharma = FALSE)


# 21) Forest plot ------------------------------------------------------

make_event_forest <- function(analysis_object, suffix) {

  # Extract sample sizes from the fitted models

  sample_sizes <- tibble::tibble(
    scenario = rep(c("All PAs", "Strict", "IUCN IV", "IUCN V", "Natura 2000"), each = 2),
    outcome = rep(c("Colonization", "Final extinction"), times = 5),
    n = c(
      nobs(analysis_object$models$all_col), nobs(analysis_object$models$all_ext),
      nobs(analysis_object$models$strict_col), nobs(analysis_object$models$strict_ext),
      nobs(analysis_object$models$iucn4_col), nobs(analysis_object$models$iucn4_ext),
      nobs(analysis_object$models$iucn5_col), nobs(analysis_object$models$iucn5_ext),
      nobs(analysis_object$models$n2k_col), nobs(analysis_object$models$n2k_ext)
    )
  )

  # Prepare model results for plotting ----------------

  plot_data <- analysis_object$effects %>%
    dplyr::left_join(sample_sizes, by = c("scenario", "outcome")) %>%
    dplyr::mutate(
      scenario_short = dplyr::recode(as.character(scenario), "All PAs" = "All"),
      scenario_plot = paste0(scenario_short, " (n = ", n, ")")
    )

  # Preserve the original x-axis limits ------------------------

  if (suffix == "withoutV") {
    x_limits <- c(0.4, 2.3)
  } else {
    x_limits <- c(0.7, 2.52)
  }

  # Build one forest-plot panel ---------------------------------

  make_event_panel <- function(data, outcome_name) {

    dat_panel <- data %>%
      dplyr::filter(outcome == outcome_name) %>%
      dplyr::mutate(
        scenario_short = factor(
          scenario_short,
          levels = c("Natura 2000", "IUCN V", "IUCN IV", "Strict", "All")
        )
      ) %>%
      dplyr::arrange(scenario_short) %>%
      dplyr::mutate(scenario_plot = factor(scenario_plot, levels = scenario_plot[order(scenario_short)]))

    ggplot2::ggplot(dat_panel, ggplot2::aes(y = scenario_plot, x = estimate, color = scenario)) +
      ggplot2::geom_vline(xintercept = 1, color = "black", linetype = "dashed", linewidth = 0.5) +
      ggplot2::geom_errorbar(
        ggplot2::aes(xmin = lcl, xmax = ucl, linetype = sig_flag),
        orientation = "y", height = 0.18, linewidth = 0.7
      ) +
      ggplot2::geom_point(shape = 21, fill = "white", size = 2.2, stroke = 1.3) +
      ggplot2::scale_color_manual(values = scenario_cols) +
      ggplot2::scale_linetype_manual(
        values = c("significant" = "solid", "not significant" = "dashed"),
        guide = "none"
      ) +
      ggplot2::scale_x_log10(
        breaks = c(0.5, 0.75, 1, 1.5, 2, 2.5),
        labels = scales::label_number(accuracy = 0.01)
      ) +
      ggplot2::coord_cartesian(xlim = x_limits) +
      ggplot2::labs(x = NULL, y = NULL) +
      ggplot2::theme_classic(base_family = base_family) +
      ggplot2::theme(
        legend.position = "none",
        axis.text.x = ggplot2::element_text(size = 34, family = base_family),
        axis.text.y = ggplot2::element_text(size = 32, family = base_family),
        plot.margin = ggplot2::margin(10, 20, 10, 10)
      )
  }

  # Create both panels --------------------------------------

  p_colonization <- make_event_panel(data = plot_data, outcome_name = "Colonization")

  p_extinction <- make_event_panel(data = plot_data, outcome_name = "Final extinction")

  # Combine panels exactly as in the current figure -----------

  p_fig <- (p_colonization / p_extinction) +
    patchwork::plot_layout(heights = c(1, 1)) +
    patchwork::plot_annotation(
      tag_levels = list(c("a) Colonization", "b) Extinction")),
      caption = "Effect ratio inside vs. outside"
    ) &
    ggplot2::theme(
      plot.tag = ggplot2::element_text(size = 36, family = base_family, hjust = 0, vjust = 1),
      plot.tag.position = c(0.02, 0.99),
      plot.caption = ggplot2::element_text(
        size = 35, family = base_family, hjust = 0.5,
        margin = ggplot2::margin(t = 10, l = 20)
      )
    )

  print(p_fig)

  # Save figure using the original dimensions ---------------

  ggplot2::ggsave(
    filename = file.path(OUT_EVENT_FIG_DIR, paste0("Figure_plotlevel_events_", suffix, "_forest.png")),
    plot = p_fig,
    width = 10, height = 9, units = "in", dpi = 280
  )

  p_fig
}


# 22) Create main and sensitivity figures ------------------

p_event_forest_without_v <- make_event_forest(analysis_object = main_without_v, suffix = "withoutV")

p_event_forest_with_v <- make_event_forest(analysis_object = sensitivity_with_v, suffix = "includingV")


# 23) Save session information -----------------------------

capture.output(sessionInfo(), file = file.path(OUT_EVENT_RESULTS_DIR, "sessionInfo_ReSurvey_event_models.txt"))
