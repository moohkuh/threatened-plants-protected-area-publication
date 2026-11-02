# 04-analyse-pa-size.R
# Supplementary protected-area size analyses
#   A) annual relative cover change
#   B) colonization counts and final-extinction risk
#
# Key decisions implemented here:
# - PA-size analyses are restricted to plots inside each focal PA category.
# - If a plot intersects >1 PA of the same category, it is assigned to
#   the largest intersecting PA for the size analysis.
# - Cover change is annualized per survey interval and THEN averaged
#   across intervals within species x plot before modelling.
# - Colonization = retained arrival not followed by later extinction.
# - Final extinction = retained loss not followed by later recolonization.
# - Colonization is modelled as a count response.
# - Final extinction is modelled as a binomial response relative to the
#   number of unique threatened species at risk.
# - Survey duration is included as a covariate.
# - "All PAs" is not analysed for PA size because a unique PA-size
#   definition is ambiguous when designation categories overlap.


# 0. SETUP -----------------------------------------------------------------

source("./R-scripts/resurvey/00-preamble-resurvey.R")

assert_files_exist(OUT_RS_CLEAN)

suppressPackageStartupMessages({
  library(dplyr)
  library(tidyr)
  library(purrr)
  library(stringr)
  library(readr)
  library(ggplot2)
  library(lme4)
  library(lmerTest)
  library(glmmTMB)
  library(DHARMa)
  library(openxlsx)
  library(tibble)
})

set.seed(42)

OUT_PA_SIZE_RESULTS_DIR <- file.path(OUT_RESULTS_DIR, "pa_size")
OUT_FIG_DIR <- OUT_FIG_RS_DIR
OUT_PA_SIZE_DIAGNOSTIC_DIR <- file.path(OUT_DIAGNOSTIC_DIR, "pa_size")
OUT_MODEL_DIR <- file.path(OUT_DIR_MODELS, "pa_size")

dir.create(OUT_PA_SIZE_RESULTS_DIR, recursive = TRUE, showWarnings = FALSE)
dir.create(OUT_FIG_DIR, recursive = TRUE, showWarnings = FALSE)
dir.create(OUT_PA_SIZE_DIAGNOSTIC_DIR, recursive = TRUE, showWarnings = FALSE)
dir.create(OUT_MODEL_DIR, recursive = TRUE, showWarnings = FALSE)

# Main colonization/extinction analysis excludes near-threatened category V.
# Set TRUE only for an additional sensitivity analysis.
INCLUDE_V_EVENTS <- FALSE
event_suffix <- ifelse(INCLUDE_V_EVENTS, "includingV", "withoutV")


# 1. LOAD CLEAN RESURVEY DATA -----------------------------------------------

rs_clean_flagged <- readr::read_csv(OUT_RS_CLEAN, show_col_types = FALSE)

# Keep nationally red-listed species.
rs_base <- rs_clean_flagged %>%
  dplyr::filter(is_nat_rl == TRUE) %>%
  dplyr::mutate(threat_2018 = stringr::str_trim(as.character(threat_2018)))

# IMPORTANT FOR PA-SIZE ANALYSES:
# Theexact Bavarian coordinates are confidential and the provider override
# contains only binary PA status (not WDPAID + PA area), that override alone
# is insufficient for a PA-size analysis.These plots were excluded
# from the PA-size analysis.


# 2. DEFINE FOCAL PA CATEGORIES ---------------------------------------------

# A row can belong to one WDPA designation; plots can nevertheless occur in
# multiple focal PA categories because different WDPA polygons overlap.

rs_pa <- rs_base %>%
  dplyr::mutate(
    pa_cat = dplyr::case_when(
      # Strict = IUCN I-III + biosphere reserve core zones.
      IUCN_label %in% strict_labels |
        grepl("core", paste0(DESIG, " ", NAME), ignore.case = TRUE) ~ "Strict",
      # Natura 2000 Habitats Directive.
      IUCN_label %in% n2k_labels ~ "Natura 2000",
      IUCN_label == "IV" ~ "IUCN IV",
      IUCN_label == "V"  ~ "IUCN V",
      TRUE ~ NA_character_
    )
  )


# 3. ASSIGN ONE PA SIZE PER PLOT x CATEGORY ---------------------------------

# For PA-size analyses, each plot needs one PA identity and one PA size per
# focal category. If a plot intersects several PAs of the same category,
# assign the largest intersecting PA.

pa_assignment <- rs_pa %>%
  dplyr::filter(!is.na(pa_cat), !is.na(WDPAID), is.finite(area_km2), area_km2 > 0) %>%
  dplyr::distinct(RS_PROJECT, RS_PROJECT_PLOT, pa_cat, WDPAID, area_km2) %>%
  dplyr::group_by(RS_PROJECT, RS_PROJECT_PLOT, pa_cat) %>%
  dplyr::slice_max(order_by = area_km2, n = 1, with_ties = FALSE) %>%
  dplyr::ungroup() %>%
  dplyr::mutate(pa_cat = factor(pa_cat, levels = c("Strict", "IUCN IV", "IUCN V", "Natura 2000")))

print(
  pa_assignment %>%
    dplyr::group_by(pa_cat) %>%
    dplyr::summarise(
      n_plots = dplyr::n_distinct(RS_PROJECT_PLOT),
      n_PAs   = dplyr::n_distinct(WDPAID),
      median_area_km2 = median(area_km2, na.rm = TRUE),
      .groups = "drop"
    )
)


# 4. COLLAPSE WDPA DUPLICATES WHILE RETAINING SURVEY INTERVALS -------------

obs_interval <- rs_base %>%
  dplyr::group_by(RS_PROJECT, RS_PLOT, RS_PROJECT_PLOT, wcvp_name, from, to) %>%
  dplyr::summarise(
    years                     = dplyr::first(resurvey_years),
    relative.change           = dplyr::first(relative.change),
    absolute.change           = dplyr::first(absolute.change),
    absolute.change.colonizer = dplyr::first(absolute.change.colonizer),
    absolute.change.extinct   = dplyr::first(absolute.change.extinct),
    SURF_AREA                 = dplyr::first(SURF_AREA),
    EUNIS_L1                  = dplyr::first(EUNIS_L1),
    threat_2018               = dplyr::first(threat_2018),
    .groups = "drop"
  ) %>%
  dplyr::filter(!is.na(years), years > 0)


# PART A. PA SIZE AND ANNUAL RELATIVE COVER CHANGE --------------------------


# 5. ANNUALIZE COVER CHANGE PER INTERVAL ------------------------------------

cover_interval <- obs_interval %>%
  dplyr::filter(!is.na(relative.change)) %>%
  dplyr::mutate(rel_change_py = relative.change / years)


# 6. JOIN PA ASSIGNMENT, THEN AVERAGE INTERVALS WITHIN SPECIES x PLOT ------

# A plot can contribute to more than one PA category if designation polygons
# overlap. For each category, however, it has only one selected WDPAID/area
# from pa_assignment.

cover_inside <- cover_interval %>%
  dplyr::inner_join(pa_assignment, by = c("RS_PROJECT", "RS_PROJECT_PLOT")) %>%
  dplyr::group_by(RS_PROJECT, RS_PROJECT_PLOT, wcvp_name, pa_cat, WDPAID, area_km2) %>%
  dplyr::summarise(
    mean_rel_change_py = mean(rel_change_py, na.rm = TRUE),
    n_intervals = dplyr::n(),
    .groups = "drop"
  ) %>%
  dplyr::filter(is.finite(mean_rel_change_py))


# 7. CENTER LOG10 PA SIZE WITHIN EACH CATEGORY ------------------------------

cover_inside <- cover_inside %>%
  dplyr::group_by(pa_cat) %>%
  dplyr::mutate(
    log_pa_area = log10(area_km2),
    log_pa_area_c = log_pa_area - mean(log_pa_area, na.rm = TRUE)
  ) %>%
  dplyr::ungroup()


# 8. FIT COVER-CHANGE PA-SIZE MODELS ----------------------------------------

fit_cover_size_model <- function(df) {
  lmerTest::lmer(
    mean_rel_change_py ~ log_pa_area_c + (1 | wcvp_name) + (1 | RS_PROJECT_PLOT) + (1 | WDPAID),
    data = df,
    REML = TRUE
  )
}

cover_split <- split(cover_inside, cover_inside$pa_cat)
cover_models <- purrr::imap(cover_split, ~ fit_cover_size_model(.x))

# Save models.
purrr::iwalk(
  cover_models,
  ~ saveRDS(.x, file.path(OUT_MODEL_DIR, paste0("cover_pa_size_", gsub("[^A-Za-z0-9]+", "_", .y), ".rds")))
)


# 9. EXTRACT COVER MODEL RESULTS --------------------------------------------

extract_cover_size <- function(model, df, scenario_label) {

  cs <- coef(summary(model))

  est <- cs["log_pa_area_c", "Estimate"]
  se  <- cs["log_pa_area_c", "Std. Error"]
  p   <- cs["log_pa_area_c", "Pr(>|t|)"]

  tibble::tibble(
    scenario = scenario_label,
    outcome  = "Cover change",
    beta     = est,
    lcl      = est - 1.96 * se,
    ucl      = est + 1.96 * se,
    p_val    = p,
    n_species = dplyr::n_distinct(df$wcvp_name),
    n_plots   = dplyr::n_distinct(df$RS_PROJECT_PLOT),
    n_PAs     = dplyr::n_distinct(df$WDPAID),
    singular  = lme4::isSingular(model, tol = 1e-4)
  )
}

cover_size_tbl <- purrr::imap_dfr(
  cover_models,
  ~ extract_cover_size(model = .x, df = cover_split[[.y]], scenario_label = .y)
) %>%
  dplyr::mutate(scenario = factor(scenario, levels = c("Strict", "IUCN IV", "IUCN V", "Natura 2000"))) %>%
  dplyr::arrange(scenario)

print(cover_size_tbl)

readr::write_csv(cover_size_tbl, file.path(OUT_PA_SIZE_RESULTS_DIR, "supplement_pa_size_cover_change.csv"))

wb_cover <- openxlsx::createWorkbook()
openxlsx::addWorksheet(wb_cover, "PA_size_cover_change")
openxlsx::writeData(wb_cover, "PA_size_cover_change", cover_size_tbl)
openxlsx::saveWorkbook(
  wb_cover, file.path(OUT_PA_SIZE_RESULTS_DIR, "supplement_pa_size_cover_change.xlsx"), overwrite = TRUE
)


# 10. COVER-MODEL DIAGNOSTICS -----------------------------------------------

cover_dharma <- purrr::imap_dfr(cover_models, ~ {
  sim <- DHARMa::simulateResiduals(.x, n = 1000, seed = 42)

  tibble::tibble(
    scenario = .y,
    dispersion_p = DHARMa::testDispersion(sim, plot = FALSE)$p.value,
    outlier_p = suppressWarnings(DHARMa::testOutliers(sim, plot = FALSE)$p.value),
    singular = lme4::isSingular(.x, tol = 1e-4)
  )
})

print(cover_dharma)

readr::write_csv(cover_dharma, file.path(OUT_PA_SIZE_DIAGNOSTIC_DIR, "cover_pa_size_DHARMa_summary.csv"))

purrr::iwalk(cover_models, ~ {
  sim <- DHARMa::simulateResiduals(.x, n = 1000, seed = 42)

  png(
    file.path(OUT_PA_SIZE_DIAGNOSTIC_DIR, paste0("cover_pa_size_", gsub("[^A-Za-z0-9]+", "_", .y), "_DHARMa.png")),
    width = 1400, height = 700, res = 130
  )

  plot(sim, main = paste("Cover change PA size -", .y))
  dev.off()
})


# PART B. PA SIZE AND COLONIZATION / EXTINCTION -----------------------------


# 11. SELECT SPECIES SET FOR EVENT ANALYSIS ---------------------------------

events_obs <- obs_interval

if (!INCLUDE_V_EVENTS) {
  events_obs <- events_obs %>%
    dplyr::filter(threat_2018 != "V")
}


# 12. DEFINE COLONIZATION AND FINAL EXTINCTION ------------------------------

# Use exactly the same retained-event definition as in the
# main plot-level colonization/extinction analysis.

events_obs <- events_obs %>%
  dplyr::mutate(
    # Raw colonization transition
    is_col = !is.na(absolute.change.colonizer),
    # Raw extinction transition
    is_ext = !is.na(absolute.change.extinct),
    # Species is assumed to be present at interval start
    # unless the interval represents a colonization
    present_at_start = !is_col
  ) %>%
  # Sort each plot-species time series chronologically
  dplyr::arrange(RS_PROJECT_PLOT, wcvp_name, from) %>%
  dplyr::group_by(RS_PROJECT_PLOT, wcvp_name) %>%
  dplyr::mutate(
    # Position of the last extinction event
    last_ext_row = ifelse(any(is_ext), max(which(is_ext)), NA_integer_),
    # TRUE if another extinction occurs later
    any_future_ext = dplyr::row_number() < last_ext_row,
    # Retain colonizations that are not followed
    # by a later extinction
    is_col_strict = dplyr::case_when(!is_col ~ FALSE, any_future_ext ~ FALSE, TRUE ~ TRUE),
    # Position of the last colonization event
    last_col_row = ifelse(any(is_col), max(which(is_col)), NA_integer_),
    # TRUE if another colonization occurs later
    any_future_col = dplyr::row_number() < last_col_row,
    # Retain extinctions that are not followed
    # by later recolonization
    is_ext_strict = dplyr::case_when(!is_ext ~ FALSE, any_future_col ~ FALSE, TRUE ~ TRUE)
  ) %>%
  dplyr::ungroup() %>%
  dplyr::select(-last_ext_row, -last_col_row)


# COLONIZATION ---------------------------------------------------------------


# 13. AGGREGATE COLONIZATION TO PLOT LEVEL ----------------------------------

plot_interval_col <- events_obs %>%
  dplyr::group_by(RS_PROJECT, RS_PLOT, RS_PROJECT_PLOT, from, to) %>%
  dplyr::summarise(
    years = dplyr::first(years), # Survey duration
    n_col = sum(is_col_strict, na.rm = TRUE), # Number of retained colonization events
    .groups = "drop"
  )

plot_level_col <- plot_interval_col %>%
  dplyr::group_by(RS_PROJECT, RS_PLOT, RS_PROJECT_PLOT) %>%
  dplyr::summarise(
    n_intervals = dplyr::n(), # Number of survey intervals
    total_years = sum(years, na.rm = TRUE), # Total observed survey duration
    total_col = sum(n_col, na.rm = TRUE), # Total retained colonizations
    .groups = "drop"
  ) %>%
  # Add plot area once per plot
  dplyr::left_join(
    obs_interval %>% dplyr::distinct(RS_PROJECT, RS_PLOT, RS_PROJECT_PLOT, SURF_AREA),
    by = c("RS_PROJECT", "RS_PLOT", "RS_PROJECT_PLOT")
  ) %>%
  dplyr::filter(!is.na(SURF_AREA), SURF_AREA > 0, total_years > 0) %>%
  dplyr::mutate(
    log_plot_area = log(SURF_AREA), # Log-transform plot area
    # Standardize plot area
    log_plot_area_z = as.numeric(scale(log_plot_area, center = TRUE, scale = TRUE)),
    # Center total survey duration
    total_years_c = as.numeric(scale(total_years, center = TRUE, scale = FALSE))
  )


# 14. JOIN PA SIZE FOR COLONIZATION -----------------------------------------

col_inside <- plot_level_col %>%
  dplyr::inner_join(pa_assignment, by = c("RS_PROJECT", "RS_PROJECT_PLOT")) %>%
  dplyr::group_by(pa_cat) %>%
  dplyr::mutate(
    # A one-unit increase represents a 10-fold increase in PA area
    log_pa_area = log10(area_km2),
    # Center PA size within each PA category
    log_pa_area_c = log_pa_area - mean(log_pa_area, na.rm = TRUE)
  ) %>%
  dplyr::ungroup()


# EXTINCTION ------------------------------------------------------------------


# 15. CREATE UNIQUE-SPECIES EXTINCTION RISK SET -----------------------------

final_loss_species <- events_obs %>%
  dplyr::group_by(RS_PROJECT, RS_PLOT, RS_PROJECT_PLOT, wcvp_name) %>%
  dplyr::summarise(
    # Species was present at the start of at least one
    # observed interval and therefore had an opportunity
    # to go locally extinct
    ever_at_risk = any(present_at_start, na.rm = TRUE),
    # Species experienced a retained final extinction
    final_extinction = any(is_ext_strict, na.rm = TRUE),
    .groups = "drop"
  )

# Check that every final extinction belongs to the risk set
bad_final_extinction <- final_loss_species %>%
  dplyr::filter(final_extinction, !ever_at_risk)

if (nrow(bad_final_extinction) > 0) {
  stop("At least one final extinction occurred for a species ", "that was never classified as at risk.")
}

# Keep only species that were at risk at least once
final_loss_species <- final_loss_species %>%
  dplyr::filter(ever_at_risk)


# 16. AGGREGATE EXTINCTION RISK SET TO PLOT LEVEL ---------------------------

plot_level_ext <- final_loss_species %>%
  dplyr::group_by(RS_PROJECT, RS_PLOT, RS_PROJECT_PLOT) %>%
  dplyr::summarise(
    n_at_risk_final = dplyr::n(), # Number of unique threatened species ever at risk
    n_final_ext = sum(final_extinction, na.rm = TRUE), # Number of unique threatened species finally lost
    .groups = "drop"
  )

# Create survey metadata for extinction models
plot_metadata_ext <- events_obs %>%
  dplyr::distinct(RS_PROJECT, RS_PLOT, RS_PROJECT_PLOT, from, to, years, SURF_AREA) %>%
  dplyr::group_by(RS_PROJECT, RS_PLOT, RS_PROJECT_PLOT) %>%
  dplyr::summarise(
    n_intervals = dplyr::n(), # Number of observed intervals
    total_years = sum(years, na.rm = TRUE), # Total observed survey duration
    SURF_AREA = dplyr::first(SURF_AREA), # Plot area
    .groups = "drop"
  )

plot_level_ext <- plot_level_ext %>%
  dplyr::left_join(plot_metadata_ext, by = c("RS_PROJECT", "RS_PLOT", "RS_PROJECT_PLOT")) %>%
  dplyr::filter(n_at_risk_final > 0, !is.na(SURF_AREA), SURF_AREA > 0, total_years > 0) %>%
  dplyr::mutate(
    log_plot_area = log(SURF_AREA), # Log-transform plot area
    # Standardize plot area
    log_plot_area_z = as.numeric(scale(log_plot_area, center = TRUE, scale = TRUE)),
    # Center total survey duration
    total_years_c = as.numeric(scale(total_years, center = TRUE, scale = FALSE))
  )


# 17. JOIN PA SIZE FOR EXTINCTION -------------------------------------------

ext_inside <- plot_level_ext %>%
  dplyr::inner_join(pa_assignment, by = c("RS_PROJECT", "RS_PROJECT_PLOT")) %>%
  dplyr::group_by(pa_cat) %>%
  dplyr::mutate(
    # A one-unit increase represents a 10-fold increase in PA area
    log_pa_area = log10(area_km2),
    # Center PA size within each PA category
    log_pa_area_c = log_pa_area - mean(log_pa_area, na.rm = TRUE)
  ) %>%
  dplyr::ungroup()


# MODELS ----------------------------------------------------------------------


# 18. SPLIT DATA BY PA CATEGORY ----------------------------------------------

col_split <- split(col_inside, col_inside$pa_cat)
ext_split <- split(ext_inside, ext_inside$pa_cat)


# 19. FIT PA-SIZE MODELS ------------------------------------------------------

# Colonization model accounting for plots nested within the same protected area
fit_col_size_model <- function(data) {
  glmmTMB::glmmTMB(
    total_col ~ log_pa_area_c + total_years_c + log_plot_area_z + (1 | RS_PROJECT) + (1 | WDPAID),
    data = data,
    family = glmmTMB::nbinom2()
  )
}

# Extinction model accounting for plots nested within the same protected area
fit_ext_size_model <- function(data) {
  glmmTMB::glmmTMB(
    cbind(n_final_ext, n_at_risk_final - n_final_ext) ~
      log_pa_area_c + total_years_c + log_plot_area_z + (1 | RS_PROJECT) + (1 | WDPAID),
    data = data,
    family = binomial()
  )
}

col_models <- purrr::imap(col_split, ~ fit_col_size_model(.x))
ext_models <- purrr::imap(ext_split, ~ fit_ext_size_model(.x))


# 20. SAVE MODELS --------------------------------------------------------------

purrr::iwalk(col_models, ~ saveRDS(
  .x,
  file.path(OUT_MODEL_DIR, paste0("colonization_pa_size_", gsub("[^A-Za-z0-9]+", "_", .y), "_", event_suffix, ".rds"))
))

purrr::iwalk(ext_models, ~ saveRDS(
  .x,
  file.path(OUT_MODEL_DIR, paste0("extinction_pa_size_", gsub("[^A-Za-z0-9]+", "_", .y), "_", event_suffix, ".rds"))
))


# EXTRACT MODEL RESULTS --------------------------------------------------------


# 21. CHECK MODEL VALIDITY ----------------------------------------------------

model_is_valid <- function(model) {

  convergence_ok <- model$fit$convergence == 0
  pd_hess_ok <- isTRUE(model$sdr$pdHess)

  se_ok <- tryCatch(
    {
      se <- summary(model)$coefficients$cond[, "Std. Error"]
      all(is.finite(se))
    },
    error = function(e) FALSE
  )

  convergence_ok && pd_hess_ok && se_ok
}


# 22. EXTRACT PA-SIZE EFFECT ---------------------------------------------------

extract_size_effect <- function(model, data, scenario_label, outcome_label, effect_type) {

  coefs <- summary(model)$coefficients$cond

  beta <- coefs["log_pa_area_c", "Estimate"]
  se   <- coefs["log_pa_area_c", "Std. Error"]
  p    <- coefs["log_pa_area_c", "Pr(>|z|)"]

  tibble::tibble(
    scenario = scenario_label,
    outcome = outcome_label,
    effect_type = effect_type,
    estimate = exp(beta),
    lcl = exp(beta - 1.96 * se),
    ucl = exp(beta + 1.96 * se),
    p_val = p,
    n_plots = dplyr::n_distinct(data$RS_PROJECT_PLOT),
    n_PAs = dplyr::n_distinct(data$WDPAID),
    valid = model_is_valid(model)
  )
}


# 23. CREATE RESULT TABLE ------------------------------------------------------

col_size_tbl <- purrr::imap_dfr(
  col_models,
  ~ extract_size_effect(model = .x, data = col_split[[.y]], scenario_label = .y, outcome_label = "Colonization", effect_type = "IRR")
)

ext_size_tbl <- purrr::imap_dfr(
  ext_models,
  ~ extract_size_effect(model = .x, data = ext_split[[.y]], scenario_label = .y, outcome_label = "Extinction", effect_type = "OR")
)

event_size_tbl <- dplyr::bind_rows(col_size_tbl, ext_size_tbl) %>%
  dplyr::mutate(
    scenario = factor(scenario, levels = c("Strict", "IUCN IV", "IUCN V", "Natura 2000")),
    outcome = factor(outcome, levels = c("Colonization", "Extinction"))
  ) %>%
  dplyr::arrange(scenario, outcome)

print(event_size_tbl)

readr::write_csv(
  event_size_tbl,
  file.path(OUT_PA_SIZE_RESULTS_DIR, paste0("supplement_pa_size_col_ext_", event_suffix, ".csv"))
)

# Save the result table as Excel
wb_events <- openxlsx::createWorkbook()
openxlsx::addWorksheet(wb_events, "PA_size_col_ext")
openxlsx::writeData(wb_events, "PA_size_col_ext", event_size_tbl)
openxlsx::saveWorkbook(
  wb_events,
  file.path(OUT_PA_SIZE_RESULTS_DIR, paste0("supplement_pa_size_col_ext_", event_suffix, ".xlsx")),
  overwrite = TRUE
)


# DIAGNOSTICS -------------------------------------------------------------------


# 24. COMBINE MODELS FOR DIAGNOSTICS ------------------------------------------

event_models <- list(Colonization = col_models, Extinction = ext_models)

event_model_flat <- purrr::imap(event_models, ~ setNames(.x, paste(names(.x), .y, sep = "__"))) %>%
  unlist(recursive = FALSE)


# 25. DHARMA DIAGNOSTICS -------------------------------------------------------

event_dharma <- purrr::imap_dfr(event_model_flat, ~ {
  sim <- DHARMa::simulateResiduals(.x, n = 1000, seed = 42)

  tibble::tibble(
    model = .y,
    uniformity_p = tryCatch(DHARMa::testUniformity(sim, plot = FALSE)$p.value, error = function(e) NA_real_),
    dispersion_p = tryCatch(DHARMa::testDispersion(sim, plot = FALSE)$p.value, error = function(e) NA_real_),
    zero_p = tryCatch(DHARMa::testZeroInflation(sim, plot = FALSE)$p.value, error = function(e) NA_real_),
    outlier_p = tryCatch(suppressWarnings(DHARMa::testOutliers(sim, plot = FALSE)$p.value), error = function(e) NA_real_),
    pdHess = isTRUE(.x$sdr$pdHess),
    valid = model_is_valid(.x)
  )
})

print(event_dharma)

readr::write_csv(
  event_dharma,
  file.path(OUT_PA_SIZE_DIAGNOSTIC_DIR, paste0("event_pa_size_DHARMa_summary_", event_suffix, ".csv"))
)

# Save diagnostic plots
purrr::iwalk(event_model_flat, ~ {
  sim <- DHARMa::simulateResiduals(.x, n = 1000, seed = 42)

  grDevices::png(
    file.path(
      OUT_PA_SIZE_DIAGNOSTIC_DIR,
      paste0("event_pa_size_", gsub("[^A-Za-z0-9]+", "_", .y), "_", event_suffix, "_DHARMa.png")
    ),
    width = 1400, height = 700, res = 130
  )

  plot(sim, main = .y)
  grDevices::dev.off()
})


# 26. SESSION INFO / FINAL CHECKS ---------------------------------------------

cat("\nCover model\n")
print(cover_size_tbl %>% dplyr::select(scenario, n_species, n_plots, n_PAs, singular))

cat("\nEvent model\n")
print(event_size_tbl %>% dplyr::select(scenario, outcome, n_plots, n_PAs, valid))

cat("\nOutput files written to:", normalizePath(OUT_PA_SIZE_RESULTS_DIR), "\n")

capture.output(sessionInfo(), file = file.path(OUT_PA_SIZE_RESULTS_DIR, "sessionInfo_ReSurvey_PA_size.txt"))
