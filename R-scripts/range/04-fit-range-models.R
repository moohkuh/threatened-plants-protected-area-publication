# 04-fit-range-models.R
# All range-scale models.
#
# PURPOSE
#   Fit the same range-persistence model structure to all four
#   analyses reported in Table S6:
#     1) T1 -> T3, main analysis, excluding Red List category V (n = 658;
#        this is the analysis reported in the main text and Figure 1)
#     2) T2 -> T3, temporal sensitivity, excluding Red List category V
#        (consistent species set with the main analysis; decision E-1,
#        29.10.2026)
#     3) T1 -> T3, sensitivity analysis including Red List category V
#        (n = 860; reported in Figure S4)
#     4) T2 -> T3, temporal sensitivity, including Red List category V
#        (n = 855; requested for Table S6 so both T1-T3 and T2-T3 are
#        reported with AND without category V, 29.10.2026)
#
#   NOTE (29.10.2026): analysis_id "T1_T3_main" previously referred to the
#   WITH-V species set and "T1_T3_without_V" to the actual main analysis
#   used in the manuscript - the opposite of what the names suggested.
#   The ids below have been swapped/relabelled so that "T1_T3_main" is
#   genuinely the main analysis (without V) and "T1_T3_with_V_sensitivity"
#   is the sensitivity analysis including V (see C-01 in
#   02_Pruefbericht_2026-10-29.md).
#
# CORE DEFINITION
#   Baseline cells are cells with OP_start > 0.
#   SOP_start and SOP_end are both calculated over those baseline cells.
#   Protected range share is the OP_start-weighted mean of PRESENT-DAY
#   PA coverage over those same baseline cells.
#
# MODEL
#   log(SOP_end / SOP_start) ~ protected range [%] + log(SOP_start)
#
# Separate linear models are fitted for:
#   All PAs, Strict, IUCN IV, IUCN V, Natura 2000 Habitats sites.
#
# No objects are taken from .GlobalEnv. This script is standalone after
# 03-build-range-dataset.R.
# QC summaries are kept in the R session and printed, but not saved as files.

suppressPackageStartupMessages({
  library(dplyr)
  library(readr)
  library(purrr)
  library(tibble)
})

OUT_DIR <- "./data/derived/range_repro"
MODEL_DIR <- file.path(OUT_DIR, "models")
dir.create(MODEL_DIR, recursive = TRUE, showWarnings = FALSE)

PATH_DATA <- file.path(OUT_DIR, "range_analysis_wide.rds")
stopifnot(file.exists(PATH_DATA))

dat <- readRDS(PATH_DATA)

required_cols <- c(
  "wcvp_name", "MTB_Q", "threat_2018", "OP_T1", "OP_T2", "OP_T3", "urban_prop",
  "cov_anyPA", "cov_strict", "cov_IV", "cov_V", "cov_n2000_habitats_only"
)

missing_cols <- setdiff(required_cols, names(dat))
if (length(missing_cols) > 0) {
  stop("Missing required columns: ", paste(missing_cols, collapse = ", "))
}

if (anyDuplicated(dat[c("wcvp_name", "MTB_Q")])) {
  stop("Input data are not unique by species x MTBQ.")
}

# New spatial pipeline guarantees 0-1 coverage; do not silently rescale.
coverage_cols <- c("cov_anyPA", "cov_strict", "cov_IV", "cov_V", "cov_n2000_habitats_only")

coverage_matrix <- as.matrix(dat[, coverage_cols])
if (any(!is.finite(coverage_matrix))) {
  stop("PA coverage contains non-finite values.")
}
if (any(coverage_matrix < 0 | coverage_matrix > 1)) {
  stop("PA coverage must be on a 0-1 scale.")
}

# Only non-urban cells enter range models.
dat_nonurban <- dat %>%
  dplyr::filter(is.finite(urban_prop), urban_prop <= 0.5)


# 1) Specifications --------------------------------------------------------

PA_SPECS <- tibble::tribble(
  ~pa_category,   ~coverage_col,
  "All PAs",      "cov_anyPA",
  "Strict",       "cov_strict",
  "IUCN IV",      "cov_IV",
  "IUCN V",       "cov_V",
  "Natura 2000",  "cov_n2000_habitats_only"
)

SPECIES_SETS <- list(
  with_V    = c("0", "1", "2", "3", "G", "R", "V"),
  without_V = c("0", "1", "2", "3", "G", "R")
)

ANALYSIS_SPECS <- tibble::tribble(
  ~analysis_id,                   ~start, ~end, ~species_set,
  "T1_T3_main",                   "T1",   "T3", "without_V",
  "T2_T3_temporal",               "T2",   "T3", "without_V",
  "T1_T3_with_V_sensitivity",     "T1",   "T3", "with_V",
  "T2_T3_with_V_sensitivity",     "T2",   "T3", "with_V"
)


# 2) Build one species-level interval dataset ------------------------------

build_interval_summary <- function(data, start, end, species_set) {

  op_start <- paste0("OP_", start)
  op_end   <- paste0("OP_", end)

  allowed_categories <- SPECIES_SETS[[species_set]]

  d <- data %>%
    dplyr::filter(threat_2018 %in% allowed_categories)

  species_universe <- d %>%
    dplyr::distinct(wcvp_name, threat_2018)

  # Baseline cells are defined ONLY by OP_start > 0.
  baseline_all <- d %>%
    dplyr::filter(is.finite(.data[[op_start]]), .data[[op_start]] > 0)

  no_baseline <- species_universe %>%
    dplyr::anti_join(baseline_all %>% dplyr::distinct(wcvp_name), by = "wcvp_name") %>%
    dplyr::mutate(
      analysis_start = start, analysis_end = end,
      species_set = species_set, exclusion_reason = "no_positive_baseline_OP"
    )

  # If any baseline cell is missing endpoint OP, exclude that species
  # explicitly rather than silently dropping individual cells.
  missing_endpoint <- baseline_all %>%
    dplyr::filter(!is.finite(.data[[op_end]])) %>%
    dplyr::distinct(wcvp_name, threat_2018) %>%
    dplyr::mutate(
      analysis_start = start, analysis_end = end,
      species_set = species_set, exclusion_reason = "missing_endpoint_OP_in_baseline_cell"
    )

  bad_endpoint_species <- missing_endpoint$wcvp_name

  baseline <- baseline_all %>%
    dplyr::filter(is.finite(.data[[op_end]]), !wcvp_name %in% bad_endpoint_species)

  # All PA coverage values must be known in baseline cells.
  if (any(!is.finite(as.matrix(baseline[, coverage_cols])))) {
    stop("Non-finite PA coverage in baseline cells for ", start, " -> ", end)
  }

  species_summary <- baseline %>%
    dplyr::group_by(wcvp_name, threat_2018) %>%
    dplyr::summarise(
      SOP_start = sum(.data[[op_start]]),
      SOP_end   = sum(.data[[op_end]]),
      n_baseline_cells = dplyr::n(),
      dplyr::across(
        dplyr::all_of(coverage_cols),
        ~ sum(.data[[op_start]] * .x) / sum(.data[[op_start]]),
        .names = "range_w_{.col}"
      ),
      .groups = "drop"
    )

  zero_endpoint <- species_summary %>%
    dplyr::filter(SOP_end <= 0 | !is.finite(SOP_end)) %>%
    dplyr::transmute(
      wcvp_name, threat_2018,
      analysis_start = start, analysis_end = end,
      species_set = species_set, exclusion_reason = "SOP_end_zero_or_nonfinite"
    )

  species_summary <- species_summary %>%
    dplyr::filter(is.finite(SOP_start), SOP_start > 0, is.finite(SOP_end), SOP_end > 0) %>%
    dplyr::mutate(
      trend = log(SOP_end / SOP_start),
      log_SOP_start = log(SOP_start),
      interval = paste0(start, "-", end),
      species_set = species_set
    )

  exclusions <- dplyr::bind_rows(no_baseline, missing_endpoint, zero_endpoint) %>%
    dplyr::distinct(wcvp_name, analysis_start, analysis_end, species_set, exclusion_reason, .keep_all = TRUE)

  list(summary = species_summary, exclusions = exclusions, n_species_universe = nrow(species_universe))
}


# 3) Fit the five PA models to one interval summary ------------------------

fit_one_pa <- function(species_summary, pa_category, coverage_col) {
  # This function receives a table with one row per species.
  # coverage_col specifies which PA coverage variable is used in the model.
  # pa_category provides the corresponding PA category name for the output.
  # This allows the same function to be used for all five PA categories.

  weighted_col <- paste0("range_w_", coverage_col) # creates the seperate weighted coverage columns

  model_data <- species_summary %>% # we create seperate model dataframes for each pa category
    dplyr::filter(
      is.finite(.data[[weighted_col]]),
      # .data[[weighted_col]] tells dplyr to use the column whose name is stored
      # in weighted_col from the current data frame.
      is.finite(trend),
      is.finite(log_SOP_start)
    ) %>%
    dplyr::transmute( # transmute just keeps columns (variables) that i need for the model and convert protected range
      # share from proportions (0-1) to percentages (0-100)
      wcvp_name, threat_2018, interval, species_set,
      SOP_start, SOP_end, trend, log_SOP_start,
      range_pct = 100 * .data[[weighted_col]]
    )
  # model fitting
  fit <- stats::lm(trend ~ range_pct + log_SOP_start, data = model_data)

  sm <- summary(fit) # save the summarys

  expected_terms <- c("range_pct", "log_SOP_start")

  # Check that both predictors were actually estimated
  coef_vec <- stats::coef(fit)

  # Standard errors directly from variance-covariance matrix
  se_vec <- sqrt(diag(stats::vcov(fit)))

  # Residual degrees of freedom for t tests and confidence intervals
  df_resid <- stats::df.residual(fit)
  crit <- stats::qt(0.975, df = df_resid) # gets the critical t from the t-distribution belonging to the residual df to get 95% CIs

  term_lookup <- c(
    "range_pct" = "range protected",
    "log_SOP_start" = paste0("log(SOP_", sub("-.*$", "", unique(model_data$interval)), ")")
  )

  result_terms <- purrr::map_dfr(expected_terms, function(term) {

    beta <- unname(coef_vec[term])
    se   <- unname(se_vec[term])

    t_value <- beta / se
    p_value <- 2 * stats::pt(abs(t_value), df = df_resid, lower.tail = FALSE)

    tibble::tibble(
      pa_category = pa_category,
      term = unname(term_lookup[term]),
      n_species = nrow(model_data),
      beta = beta,
      se = se,
      conf_low = beta - crit * se, # formula to calculate CIs
      conf_high = beta + crit * se,
      t_value = t_value,
      p_value = p_value,
      r_squared = sm$r.squared,
      adj_r_squared = sm$adj.r.squared
    )
  })

  # With two predictors, VIF for range_pct is 1/(1-r^2) from the
  # pairwise predictor correlation.
  pred_r <- stats::cor(model_data$range_pct, model_data$log_SOP_start, use = "complete.obs", method = "pearson")
  vif_range <- 1 / (1 - pred_r^2)

  cor_tbl <- tibble::tibble(
    pa_category = pa_category,
    pearson_r = pred_r,
    spearman_rho = stats::cor(model_data$range_pct, model_data$log_SOP_start, use = "complete.obs", method = "spearman"),
    vif_range = vif_range
  )

  # Prediction line at median baseline range size.
  pred_grid <- tibble::tibble(
    range_pct = seq(min(model_data$range_pct, na.rm = TRUE), max(model_data$range_pct, na.rm = TRUE), length.out = 200),
    log_SOP_start = stats::median(model_data$log_SOP_start, na.rm = TRUE)
  )

  pr <- stats::predict(fit, newdata = pred_grid, se.fit = TRUE)

  pred_grid <- pred_grid %>%
    dplyr::mutate(
      predicted = as.numeric(pr$fit),
      se_fit = as.numeric(pr$se.fit),
      conf_low = predicted - 1.96 * se_fit,
      conf_high = predicted + 1.96 * se_fit,
      pa_category = pa_category
    )

  list(fit = fit, data = model_data, results = result_terms, correlation = cor_tbl, prediction = pred_grid)
}


# 4) Run all three analyses with exactly the same functions -----------------

all_objects <- list()
all_results <- list()
all_species_summaries <- list()
all_exclusions <- list()
all_correlations <- list()
all_predictions <- list()
analysis_qc <- list()

for (i in seq_len(nrow(ANALYSIS_SPECS))) {

  spec <- ANALYSIS_SPECS[i, ]

  interval_obj <- build_interval_summary(
    data = dat_nonurban, start = spec$start, end = spec$end, species_set = spec$species_set
  )

  species_summary <- interval_obj$summary

  message(spec$analysis_id, ": ", nrow(species_summary), " species retained from ", interval_obj$n_species_universe)

  # Fit the same five PA categories for all retained analyses
  fits <- purrr::map2(
    PA_SPECS$pa_category, PA_SPECS$coverage_col,
    ~ fit_one_pa(species_summary = species_summary, pa_category = .x, coverage_col = .y)
  )

  names(fits) <- PA_SPECS$pa_category

  results_i <- purrr::map_dfr(fits, "results") %>%
    dplyr::mutate(
      analysis_id = spec$analysis_id, interval = paste0(spec$start, "-", spec$end),
      species_set = spec$species_set, .before = 1
    )

  correlations_i <- purrr::map_dfr(fits, "correlation") %>%
    dplyr::mutate(
      analysis_id = spec$analysis_id, interval = paste0(spec$start, "-", spec$end),
      species_set = spec$species_set, .before = 1
    )

  predictions_i <- purrr::map_dfr(fits, "prediction") %>%
    dplyr::mutate(
      analysis_id = spec$analysis_id, interval = paste0(spec$start, "-", spec$end),
      species_set = spec$species_set, .before = 1
    )

  all_objects[[spec$analysis_id]] <- fits
  all_results[[spec$analysis_id]] <- results_i
  all_species_summaries[[spec$analysis_id]] <- species_summary %>%
    dplyr::mutate(analysis_id = spec$analysis_id, .before = 1)
  all_exclusions[[spec$analysis_id]] <- interval_obj$exclusions %>%
    dplyr::mutate(analysis_id = spec$analysis_id, .before = 1)
  all_correlations[[spec$analysis_id]] <- correlations_i
  all_predictions[[spec$analysis_id]] <- predictions_i

  analysis_qc[[spec$analysis_id]] <- tibble::tibble(
    analysis_id = spec$analysis_id,
    interval = paste0(spec$start, "-", spec$end),
    species_set = spec$species_set,
    n_species_universe = interval_obj$n_species_universe,
    n_species_retained = nrow(species_summary),
    n_species_excluded = dplyr::n_distinct(interval_obj$exclusions$wcvp_name),
    n_negative_trends = sum(species_summary$trend < 0),
    pct_negative_trends = 100 * mean(species_summary$trend < 0),
    median_trend = stats::median(species_summary$trend),
    mean_trend = mean(species_summary$trend),
    sd_trend = stats::sd(species_summary$trend)
  )
}

results_all <- dplyr::bind_rows(all_results)
species_summaries_all <- dplyr::bind_rows(all_species_summaries)
exclusions_all <- dplyr::bind_rows(all_exclusions)
correlations_all <- dplyr::bind_rows(all_correlations)
predictions_all <- dplyr::bind_rows(all_predictions)
analysis_qc <- dplyr::bind_rows(analysis_qc)


# 5) Save model outputs used for the publication -----------------------------

saveRDS(all_objects, file.path(MODEL_DIR, "range_lm_objects.rds"))
saveRDS(all_species_summaries, file.path(MODEL_DIR, "range_species_summaries.rds"))
saveRDS(all_predictions, file.path(MODEL_DIR, "range_prediction_grids.rds"))

readr::write_csv(results_all, file.path(OUT_DIR, "range_model_results_all_terms.csv"))

readr::write_csv(
  results_all %>% dplyr::filter(term == "range protected"),
  file.path(OUT_DIR, "range_model_results_PA_term.csv")
)

readr::write_csv(species_summaries_all, file.path(OUT_DIR, "range_species_summaries_all_analyses.csv"))


# 6) Share of declining species in the main analysis --------------------------

# Number of declining species out of total number of threatened species
range_species <- read_csv("./data/derived/range_repro/range_species_summaries_all_analyses.csv", show_col_types = FALSE)

# Count declining species in the analysis excluding category V
decline_main <- range_species %>%
  filter(analysis_id == "T1_T3_main") %>%
  summarise(
    n_species = n(),
    n_declining = sum(trend < 0),
    percent_declining = round(100 * n_declining / n_species, 1)
  )

print(decline_main)


# 7) Save R/package environment for reproducibility ----------------------------

capture.output(sessionInfo(), file = file.path(OUT_DIR, "sessionInfo_range_analysis.txt"))
