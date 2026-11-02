# 01-build-resurvey-data.R
# Build the clean ReSurvey Germany dataset used by all analyses.
#
# Main steps:
#   1. Load ReSurvey Germany change data.
#   2. Attach public plot metadata and frozen PA classification.
#   3. Harmonize taxa to WCVP names and attach national Red List status.
#   4. Attach frozen plot-to-PA memberships for PA-size analyses.
#   5. Save one clean table for downstream analyses.


# 0) Setup ---------------------------------------------------------------

source("./R-scripts/resurvey/00-preamble-resurvey.R")

assert_files_exist(c(
  PATH_RS_RAW, PATH_RS_PLOT_CLASSIFICATION, PATH_RS_PA_MEMBERSHIP,
  PATH_WCVP_MATCH_RS, PATH_RL_NAT_MATCH
))


# 1) Load raw ReSurvey Germany data --------------------------------------

rs_raw <- readr::read_csv(PATH_RS_RAW, show_col_types = FALSE)

if (!"RS_PROJECT_PLOT" %in% names(rs_raw)) {
  rs_raw <- rs_raw %>%
    dplyr::mutate(RS_PROJECT_PLOT = paste(RS_PROJECT, RS_PLOT, sep = "_"))
}


# 2) Attach public plot metadata and frozen PA classification -------------

plot_classification <- readr::read_csv(PATH_RS_PLOT_CLASSIFICATION, show_col_types = FALSE)

assert_columns(plot_classification, c(
  "RS_PROJECT_PLOT", "LONGITUDE", "LATITUDE", "SURF_AREA", "EUNIS_plot", "EUNIS_L1", "EUNIS_L2",
  "any_PA", "any_strict", "any_iucn4", "any_iucn5", "any_n2k"
))

rs_with_plot_data <- rs_raw %>%
  dplyr::left_join(plot_classification, by = "RS_PROJECT_PLOT") %>%
  dplyr::mutate(
    resurvey_years = dplyr::if_else(!is.na(from) & !is.na(to), as.numeric(to - from), NA_real_)
  )


# 3) Attach harmonized taxonomy ------------------------------------------

# Load the ReSurvey -> WCVP lookup created by the shared taxonomy pipeline.
rs_taxonomy_lookup <- readr::read_csv(PATH_WCVP_MATCH_RS, show_col_types = FALSE) %>%
  dplyr::mutate(species = stringr::str_squish(species), wcvp_name = stringr::str_squish(wcvp_name))

rs_taxonomy <- rs_with_plot_data %>%
  dplyr::mutate(species = stringr::str_squish(species)) %>%
  dplyr::left_join(rs_taxonomy_lookup, by = "species")


# 4) Attach national Red List status -------------------------------------

# Load the already cleaned national Red List lookup created by the
# shared taxonomy pipeline.
nat_rl_lookup <- readr::read_csv(PATH_RL_NAT_MATCH, show_col_types = FALSE) %>%
  dplyr::mutate(wcvp_name = stringr::str_squish(wcvp_name))

rs_taxonomy <- rs_taxonomy %>%
  dplyr::left_join(nat_rl_lookup, by = "wcvp_name") %>%
  dplyr::mutate(is_nat_rl = dplyr::coalesce(is_nat_rl, FALSE))


# 5) Attach plot-to-PA memberships ---------------------------------------

# Attach the frozen plot-to-PA memberships used for the PA-size analysis.
plot_pa_membership <- readr::read_csv(PATH_RS_PA_MEMBERSHIP, show_col_types = FALSE)

assert_columns(plot_pa_membership, c(
  "RS_PROJECT_PLOT", "WDPAID", "NAME", "DESIG", "IUCN_CAT", "IUCN_label", "area_km2"
))

RS_clean_flagged <- rs_taxonomy %>%
  dplyr::left_join(plot_pa_membership, by = "RS_PROJECT_PLOT")


# 6) Basic QC -------------------------------------------------------------

qc_summary <- RS_clean_flagged %>%
  dplyr::summarise(
    n_rows = dplyr::n(),
    n_plots = dplyr::n_distinct(RS_PROJECT_PLOT),
    n_taxa = dplyr::n_distinct(wcvp_name),
    n_nat_rl_rows = sum(is_nat_rl, na.rm = TRUE),
    n_rows_in_pa = sum(any_PA, na.rm = TRUE),
    n_rows_missing_eunis = sum(is.na(EUNIS_plot))
  )

print(qc_summary)


# 7) Table S7: overlapping PA types per plot ------------------------------

# Number of focal PA categories (Strict, IUCN IV, IUCN V, Natura 2000)
# intersecting each plot, from the frozen plot classification.
plot_n_pa_types <- plot_classification %>%
  dplyr::mutate(
    n_pa_types = as.integer(any_strict) + as.integer(any_iucn4) +
      as.integer(any_iucn5) + as.integer(any_n2k)
  ) %>%
  dplyr::select(RS_PROJECT_PLOT, n_pa_types)

table_s7 <- RS_clean_flagged %>%
  dplyr::left_join(plot_n_pa_types, by = "RS_PROJECT_PLOT") %>%
  dplyr::group_by(n_pa_types) %>%
  dplyr::summarise(
    n_plots = dplyr::n_distinct(RS_PROJECT_PLOT),
    n_species_plot_interval_records = dplyr::n(),
    .groups = "drop"
  ) %>%
  dplyr::arrange(n_pa_types)

print(table_s7)
readr::write_csv(table_s7, file.path(OUT_RESULTS_DIR, "TableS7_pa_types_per_plot.csv"))


# 8) Save clean dataset ---------------------------------------------------

readr::write_csv(RS_clean_flagged, OUT_RS_CLEAN)

capture.output(sessionInfo(), file = file.path(OUT_RS_DIR, "sessionInfo_build_resurvey_data.txt"))
