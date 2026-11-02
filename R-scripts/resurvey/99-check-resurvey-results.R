# 99-check-resurvey-results.R
# Compare key ReSurvey outputs with manuscript/SI values.
#
# This script does not refit models. It reads result tables created by
# the pipeline and checks whether key rounded estimates reproduce the
# values reported in the current supplementary information.

source("./R-scripts/resurvey/00-preamble-resurvey.R")


# 1) Expected protection-status effects ---------------------------------

expected_without_v <- tibble::tribble(
  ~scenario,     ~outcome,            ~estimate_expected,
  "All PAs",     "Colonization",       1.124,
  "All PAs",     "Final extinction",   0.864,
  "Strict",      "Colonization",       1.517,
  "Strict",      "Final extinction",   0.997,
  "IUCN IV",     "Colonization",       1.238,
  "IUCN IV",     "Final extinction",   0.651,
  "IUCN V",      "Colonization",       1.024,
  "IUCN V",      "Final extinction",   1.323,
  "Natura 2000", "Colonization",       1.335,
  "Natura 2000", "Final extinction",   0.634
)

expected_including_v <- tibble::tribble(
  ~scenario,     ~outcome,            ~estimate_expected,
  "All PAs",     "Colonization",       1.353,
  "All PAs",     "Final extinction",   1.063,
  "Strict",      "Colonization",       1.251,
  "Strict",      "Final extinction",   1.554,
  "IUCN IV",     "Colonization",       1.668,
  "IUCN IV",     "Final extinction",   0.814,
  "IUCN V",      "Colonization",       1.212,
  "IUCN V",      "Final extinction",   1.149,
  "Natura 2000", "Colonization",       1.653,
  "Natura 2000", "Final extinction",   0.937
)

# Compares one event-model result file (withoutV / includingV) against
# the matching expected table above.
check_event_table <- function(suffix, expected_table) {

  result_file <- file.path(OUT_RESULTS_DIR, "events", paste0("event_effects_", suffix, ".csv"))
  assert_files_exist(result_file)

  observed <- readr::read_csv(result_file, show_col_types = FALSE) %>%
    dplyr::mutate(scenario = as.character(scenario), outcome = as.character(outcome)) %>%
    dplyr::select(scenario, outcome, estimate)

  expected_table %>%
    dplyr::left_join(observed, by = c("scenario", "outcome")) %>%
    dplyr::mutate(
      estimate_observed = estimate,
      difference = estimate_observed - estimate_expected,
      # reported values have three decimals, so the difference
      # must be smaller than half of the last digit
      matches_to_3_decimals = abs(difference) <= 0.0005
    ) %>%
    dplyr::select(scenario, outcome, estimate_expected, estimate_observed, difference, matches_to_3_decimals)
}

event_check_without_v <- check_event_table("withoutV", expected_without_v)
event_check_including_v <- check_event_table("includingV", expected_including_v)


# 2) Expected PA-size cover effects --------------------------------------

expected_pa_size_cover <- tibble::tribble(
  ~scenario,     ~beta_expected,
  "Strict",      -0.089,
  "IUCN IV",      0.002,
  "IUCN V",      -0.006,
  "Natura 2000",  0.004
)

pa_size_cover_file <- file.path(OUT_RESULTS_DIR, "pa_size", "supplement_pa_size_cover_change.csv")
assert_files_exist(pa_size_cover_file)

pa_size_cover_check <- readr::read_csv(pa_size_cover_file, show_col_types = FALSE) %>%
  dplyr::mutate(scenario = as.character(scenario)) %>%
  dplyr::left_join(expected_pa_size_cover, by = "scenario") %>%
  dplyr::mutate(
    beta_observed = beta,
    matches_to_3_decimals = abs(beta_observed - beta_expected) <= 0.0005
  ) %>%
  dplyr::select(scenario, beta_expected, beta_observed, matches_to_3_decimals)


# 3) Print checks ---------------------------------------------------------

cat("\nMain event analysis excluding category V:\n")
print(event_check_without_v)

cat("\nEvent sensitivity including category V:\n")
print(event_check_including_v)

cat("\nPA-size cover analysis:\n")
print(pa_size_cover_check)

all_event_matches <- all(event_check_without_v$matches_to_3_decimals, event_check_including_v$matches_to_3_decimals)
all_pa_size_cover_matches <- all(pa_size_cover_check$matches_to_3_decimals)

if (all_event_matches && all_pa_size_cover_matches) {
  cat("\nAll checked estimates reproduce the reported values ", "to three decimal places.\n")
} else {
  warning(
    "At least one checked estimate differs from the reported value. ",
    "Inspect the comparison tables above before proceeding."
  )
}
