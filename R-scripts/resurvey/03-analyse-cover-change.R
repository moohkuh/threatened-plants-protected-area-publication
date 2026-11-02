# 03-analyse-cover-change.R
# Species-level cover change inside vs. outside protected areas (Fig. 2, Tables S8-S11)
#
# Main analysis: Red List categories 0-3, G, R and V.
# Relative cover change is annualized per survey interval, averaged within
# species x plot and then within species x side (inside/outside).
# Only species with >= 5 plots inside and >= 5 plots outside are kept.
# Sensitivity: paired Wilcoxon tests and LMMs without category V.


# 0) Setup ---------------------------------------------------------------

source("./R-scripts/resurvey/00-preamble-resurvey.R")

suppressPackageStartupMessages({
  library(openxlsx)
})

assert_files_exist(c(OUT_RS_CLEAN))

MIN_PLOTS <- 5

RDS_FILE <- file.path(OUT_DIR_RDS, "species_cover_means_objects.rds")
TABLE_FILE <- file.path(OUT_RESULTS_DIR, "species_cover_means_model_tables.xlsx")


# 1) Red-listed species only ---------------------------------------------

rs_clean <- readr::read_csv(OUT_RS_CLEAN, show_col_types = FALSE) %>%
  dplyr::filter(is_nat_rl == TRUE)


# 2) One row per plot x species x interval -------------------------------

# a plot can lie in several PA polygons, which duplicates rows after the
# spatial join; all values are identical within a group, so first() is fine
obs_interval <- rs_clean %>%
  dplyr::group_by(RS_PROJECT, RS_PLOT, RS_PROJECT_PLOT, wcvp_name, threat_2018, from, to) %>%
  dplyr::summarise(
    years = dplyr::first(resurvey_years),
    relative.change = dplyr::first(relative.change),
    absolute.change = dplyr::first(absolute.change),
    absolute.change.colonizer = dplyr::first(absolute.change.colonizer),
    absolute.change.extinct = dplyr::first(absolute.change.extinct),
    LONGITUDE = dplyr::first(LONGITUDE),
    LATITUDE = dplyr::first(LATITUDE),
    any_PA = dplyr::first(any_PA),
    any_strict = dplyr::first(any_strict),
    any_n2k = dplyr::first(any_n2k),
    any_iucn4 = dplyr::first(any_iucn4),
    any_iucn5 = dplyr::first(any_iucn5),
    .groups = "drop"
  ) %>%
  dplyr::mutate(rel_change_py = relative.change / years) %>%
  dplyr::filter(!is.na(years), years > 0, is.finite(rel_change_py))


# 3) Function: data for one inside/outside comparison --------------------

build_cover_scenario <- function(data, scenario, inside_flag = NULL, exclude_v = FALSE) {

  scenario_data <- data

  if (exclude_v) {
    scenario_data <- scenario_data %>% dplyr::filter(threat_2018 != "V")
  }

  # All PAs: inside = any PA. Other categories: inside = this category,
  # outside = no PA at all (make_pa_scenario() in the preamble)
  if (scenario == "All PAs") {
    scenario_data <- scenario_data %>%
      dplyr::mutate(
        side = dplyr::if_else(any_PA, "inside", "outside"),
        side = factor(side, levels = c("outside", "inside"))
      )
  } else {
    scenario_data <- make_pa_scenario(scenario_data, inside_flag)
  }

  # keep species with at least MIN_PLOTS plots on both sides
  species_keep <- scenario_data %>%
    dplyr::group_by(wcvp_name, side) %>%
    dplyr::summarise(n_plots = dplyr::n_distinct(RS_PROJECT_PLOT), .groups = "drop") %>%
    tidyr::pivot_wider(names_from = side, values_from = n_plots, values_fill = 0) %>%
    dplyr::filter(inside >= MIN_PLOTS, outside >= MIN_PLOTS) %>%
    dplyr::pull(wcvp_name)

  scenario_data <- scenario_data %>% dplyr::filter(wcvp_name %in% species_keep)

  # step 1: mean over intervals within species x plot
  plot_means <- scenario_data %>%
    dplyr::group_by(wcvp_name, side, RS_PROJECT_PLOT) %>%
    dplyr::summarise(rel_change_py_plot = mean(rel_change_py, na.rm = TRUE), .groups = "drop")

  # step 2: mean over plots within species x side
  species_means <- plot_means %>%
    dplyr::group_by(wcvp_name, side) %>%
    dplyr::summarise(mean_rel = mean(rel_change_py_plot, na.rm = TRUE), .groups = "drop")

  # every species needs exactly one inside and one outside value
  pairing_check <- species_means %>% dplyr::count(wcvp_name)

  if (any(pairing_check$n != 2)) {
    stop("Paired inside/outside structure failed for scenario: ", scenario)
  }

  list(long = scenario_data, plot_means = plot_means, species_means = species_means)
}


scenario_specs <- tibble::tribble(
  ~scenario,     ~inside_flag,
  "All PAs",     NA_character_,
  "Strict",      "any_strict",
  "IUCN IV",     "any_iucn4",
  "IUCN V",      "any_iucn5",
  "Natura 2000", "any_n2k"
)


# 4) Build the five comparisons ------------------------------------------

cover_main <- purrr::pmap(
  scenario_specs,
  function(scenario, inside_flag) {
    build_cover_scenario(data = obs_interval, scenario = scenario,
                         inside_flag = inside_flag, exclude_v = FALSE)
  }
)

names(cover_main) <- scenario_specs$scenario


# 5) Paired LMMs (species as random intercept) ---------------------------

fit_cover_lmm <- function(species_means, scenario) {

  model <- lmerTest::lmer(mean_rel ~ side + (1 | wcvp_name), data = species_means, REML = TRUE)

  # estimated means per side (Table S8)
  emm <- as.data.frame(confint(emmeans::emmeans(model, ~ side))) %>%
    dplyr::transmute(
      scenario = scenario,
      side = factor(side, levels = c("outside", "inside")),
      mean = emmean,
      lower = lower.CL,
      upper = upper.CL
    )

  # inside - outside contrast with CI and p-value
  contrast <- as.data.frame(
    summary(
      emmeans::contrast(emmeans::emmeans(model, ~ side),
                        method = list("inside - outside" = c(-1, 1))),
      infer = c(TRUE, TRUE)
    )
  ) %>%
    dplyr::transmute(
      scenario = scenario,
      estimate = estimate,
      se = SE,
      lower = lower.CL,
      upper = upper.CL,
      p_value = p.value
    )

  list(model = model, emm = emm, contrast = contrast)
}


cover_fits <- purrr::imap(cover_main, ~ fit_cover_lmm(species_means = .x$species_means, scenario = .y))

emm_tbl <- purrr::map_dfr(cover_fits, "emm")

effects_tbl <- purrr::map_dfr(cover_fits, "contrast") %>%
  dplyr::mutate(
    scenario = factor(scenario, levels = c("Natura 2000", "IUCN V", "IUCN IV", "Strict", "All PAs"))
  )


# 6) Sample sizes and direction of change (Tables S9, S11) ---------------

sample_size_tbl <- purrr::imap_dfr(
  cover_main,
  function(scenario_object, scenario_name) {

    plot_counts <- scenario_object$long %>%
      dplyr::group_by(side) %>%
      dplyr::summarise(n_plots = dplyr::n_distinct(RS_PROJECT_PLOT), .groups = "drop") %>%
      tidyr::pivot_wider(names_from = side, values_from = n_plots,
                         values_fill = 0, names_prefix = "n_plots_")

    tibble::tibble(
      scenario = scenario_name,
      n_species = dplyr::n_distinct(scenario_object$species_means$wcvp_name)
    ) %>%
      dplyr::bind_cols(plot_counts)
  }
)

median_tbl <- purrr::imap_dfr(
  cover_main,
  function(scenario_object, scenario_name) {
    scenario_object$species_means %>%
      dplyr::group_by(side) %>%
      dplyr::summarise(median = median(mean_rel, na.rm = TRUE), .groups = "drop") %>%
      dplyr::mutate(scenario = scenario_name)
  }
)

species_diff_tbl <- purrr::imap_dfr(
  cover_main,
  function(scenario_object, scenario_name) {
    scenario_object$species_means %>%
      tidyr::pivot_wider(names_from = side, values_from = mean_rel) %>%
      dplyr::mutate(diff_inside_outside = inside - outside) %>%
      dplyr::summarise(
        scenario = scenario_name,
        proportion_species_inside_gt_outside = mean(diff_inside_outside > 0, na.rm = TRUE),
        proportion_species_negative_inside = mean(inside < 0, na.rm = TRUE),
        proportion_species_negative_outside = mean(outside < 0, na.rm = TRUE),
        median_diff = median(diff_inside_outside, na.rm = TRUE),
        IQR_diff = IQR(diff_inside_outside, na.rm = TRUE)
      )
  }
)


# 7) Paired Wilcoxon tests (Table S10) -----------------------------------

make_wilcoxon <- function(species_means, scenario) {

  paired_data <- species_means %>%
    tidyr::pivot_wider(names_from = side, values_from = mean_rel)

  test <- stats::wilcox.test(paired_data$inside, paired_data$outside,
                             paired = TRUE, exact = FALSE, conf.int = TRUE)

  tibble::tibble(
    scenario = scenario,
    n_species = nrow(paired_data),
    V_statistic = as.numeric(test$statistic),
    pseudo_median = as.numeric(test$estimate),
    ci_lower = as.numeric(test$conf.int[1]),
    ci_upper = as.numeric(test$conf.int[2]),
    p_value = test$p.value
  )
}

wilcox_tbl <- purrr::imap_dfr(cover_main, ~ make_wilcoxon(.x$species_means, .y)) %>%
  dplyr::mutate(scenario = factor(scenario, levels = scenario_levels))

readr::write_csv(wilcox_tbl, file.path(OUT_RESULTS_DIR, "cover_change_wilcoxon_signed_rank.csv"))


# 8) Sensitivity without category V (Figure S5) --------------------------

cover_no_v <- purrr::pmap(
  scenario_specs,
  function(scenario, inside_flag) {

    # Strict is left out: too few species for a stable model without V
    if (scenario == "Strict") {
      return(NULL)
    }

    build_cover_scenario(data = obs_interval, scenario = scenario,
                         inside_flag = inside_flag, exclude_v = TRUE)
  }
)

names(cover_no_v) <- scenario_specs$scenario

cover_no_v <- cover_no_v[!purrr::map_lgl(cover_no_v, is.null)]

cover_no_v_fits <- purrr::imap(cover_no_v, ~ fit_cover_lmm(species_means = .x$species_means, scenario = .y))

effects_no_v <- purrr::imap_dfr(
  cover_no_v_fits,
  function(fit_object, scenario_name) {
    fit_object$contrast %>%
      dplyr::mutate(n_species = dplyr::n_distinct(cover_no_v[[scenario_name]]$species_means$wcvp_name))
  }
) %>%
  dplyr::mutate(scenario = factor(scenario, levels = c("Natura 2000", "IUCN V", "IUCN IV", "All PAs")))

readr::write_csv(effects_no_v, file.path(OUT_RESULTS_DIR, "cover_change_sensitivity_without_V.csv"))


# 9) Model diagnostics (DHARMa) ------------------------------------------

diagnostic_summary <- purrr::imap_dfr(
  cover_fits,
  function(fit_object, scenario_name) {

    simulation <- DHARMa::simulateResiduals(fit_object$model, n = 1000, seed = 42)

    tibble::tibble(
      scenario = scenario_name,
      dispersion_p = DHARMa::testDispersion(simulation, plot = FALSE)$p.value,
      outlier_p = suppressWarnings(DHARMa::testOutliers(simulation, plot = FALSE)$p.value),
      singular = lme4::isSingular(fit_object$model, tol = 1e-4)
    )
  }
)

readr::write_csv(diagnostic_summary, file.path(OUT_DIAGNOSTIC_DIR, "cover_change_DHARMa_summary.csv"))


# 10) Save tables --------------------------------------------------------

workbook <- openxlsx::createWorkbook()

for (sheet_name in c("EMMs", "Contrasts_LMM", "Wilcoxon", "Medians",
                     "Species_diff", "Sample_sizes", "Sensitivity_without_V")) {
  openxlsx::addWorksheet(workbook, sheet_name)
}

openxlsx::writeData(workbook, "EMMs", emm_tbl)
openxlsx::writeData(workbook, "Contrasts_LMM", effects_tbl)
openxlsx::writeData(workbook, "Wilcoxon", wilcox_tbl)
openxlsx::writeData(workbook, "Medians", median_tbl)
openxlsx::writeData(workbook, "Species_diff", species_diff_tbl)
openxlsx::writeData(workbook, "Sample_sizes", sample_size_tbl)
openxlsx::writeData(workbook, "Sensitivity_without_V", effects_no_v)

openxlsx::saveWorkbook(workbook, TABLE_FILE, overwrite = TRUE)


# 11) Cache for Figure 2, Figure S5 and 03b ------------------------------

saveRDS(
  list(
    obs_interval = obs_interval,
    dat_all_long = cover_main[["All PAs"]]$long,
    dat_strict_long = cover_main[["Strict"]]$long,
    dat_iucn4_long = cover_main[["IUCN IV"]]$long,
    dat_iucn5_long = cover_main[["IUCN V"]]$long,
    dat_n2k_long = cover_main[["Natura 2000"]]$long,
    sp_all = cover_main[["All PAs"]]$species_means,
    sp_strict = cover_main[["Strict"]]$species_means,
    sp_iucn4 = cover_main[["IUCN IV"]]$species_means,
    sp_iucn5 = cover_main[["IUCN V"]]$species_means,
    sp_n2k = cover_main[["Natura 2000"]]$species_means,
    m_all = cover_fits[["All PAs"]]$model,
    m_strict = cover_fits[["Strict"]]$model,
    m_iucn4 = cover_fits[["IUCN IV"]]$model,
    m_iucn5 = cover_fits[["IUCN V"]]$model,
    m_n2k = cover_fits[["Natura 2000"]]$model,
    emm_tbl = emm_tbl,
    effects_tbl = effects_tbl,
    wilcox_tbl = wilcox_tbl,
    median_tbl = median_tbl,
    species_diff_tbl = species_diff_tbl,
    sample_size_tbl = sample_size_tbl,
    effects_no_v = effects_no_v
  ),
  RDS_FILE
)


# 12) Main results -------------------------------------------------------

print(effects_tbl)
print(wilcox_tbl)
print(effects_no_v)

capture.output(sessionInfo(), file = file.path(OUT_RESULTS_DIR, "sessionInfo_ReSurvey_cover_change.txt"))
