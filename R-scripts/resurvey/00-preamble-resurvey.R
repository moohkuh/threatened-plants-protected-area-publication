# 00-preamble-resurvey.R
# Shared packages, paths, constants and helper functions
# for the ReSurvey Germany publication pipeline.

suppressPackageStartupMessages({
  library(dplyr)
  library(tidyr)
  library(stringr)
  library(readr)
  library(sf)
  library(lme4)
  library(lmerTest)
  library(emmeans)
  library(ggplot2)
  library(scales)
  library(purrr)
  library(tibble)
  library(glmmTMB)
  library(DHARMa)
})

sf::sf_use_s2(FALSE)
set.seed(42)


# Input paths ------------------------------------------------------

# These paths intentionally retain the frozen input files used for the
# manuscript. After the numerical results have been reproduced, the
# input files can be moved into a cleaner data/raw structure by changing
# only the paths in this preamble.

PATH_RS_VEGETATION <- "./data/raw/resurvey/ReSurveyGermany.csv"
PATH_RS_PLOT_CLASSIFICATION <- "./data/resurvey/resurvey_plot_classification.csv"
PATH_RS_PA_MEMBERSHIP <- "./data/resurvey/resurvey_plot_pa_membership.csv"

# Derived ReSurvey change table
PATH_RS_RAW <- "./data/derived/resurvey/RS_raw.csv"

# Taxonomic lookup tables created by the shared taxonomy pipeline
PATH_WCVP_MATCH_RS <- "./data/derived/taxonomy/resurvey_wcvp_lookup.csv"
PATH_RL_NAT_MATCH <- "./data/derived/taxonomy/national_red_list_lookup.csv"

PATH_PA_GPKG <- "./data/derived/protected_areas/PA_polygons_slim.gpkg"
LAYER_PA <- "pa_polygons_slim"


# Output paths --------------------------------------------------

OUT_RS_DIR <- "./data/derived/resurvey_repro"
OUT_RS_CLEAN <- file.path(OUT_RS_DIR, "resurvey_data_clean.csv")
OUT_DIR_RDS <- file.path(OUT_RS_DIR, "cache")
OUT_DIR_MODELS <- file.path(OUT_RS_DIR, "models")

OUT_RESULTS_DIR <- "./results/resurvey"
OUT_FIG_RS_DIR <- "./figures/resurvey"
OUT_DIAGNOSTIC_DIR <- "./diagnostics/resurvey"

for (directory in c(OUT_RS_DIR, OUT_DIR_RDS, OUT_DIR_MODELS, OUT_RESULTS_DIR, OUT_FIG_RS_DIR, OUT_DIAGNOSTIC_DIR)) {
  dir.create(directory, recursive = TRUE, showWarnings = FALSE)
}


# Protected-area definitions -------------------------------

# Protected-area labels used in the shared PA layer
strict_labels <- c("Ia/Ib", "II", "BR core", "III / National Natural Monument")
n2k_labels <- "N2000 Habitats"

scenario_levels <- c("All PAs", "Strict", "IUCN IV", "IUCN V", "Natura 2000")
outcome_levels <- c("Colonization", "Final extinction")


# Figure settings ----------------------------------------

# Register Roboto Condensed for consistent figure rendering.
library(sysfonts)
library(showtext)

font_add_google("Roboto Condensed", "roboto_condensed")
showtext_auto()
base_family <- "roboto_condensed"

col_outside <- "#6E6E6E"
col_strict  <- "#8AAE8AFF"
col_n2k     <- "#E5A729FF"
col_iucn4   <- "#4F4D1DFF"
col_iucn5   <- "#A08F27FF"
col_allpa   <- "#5D7388FF"
scenario_cols <- c(
  "All PAs" = col_allpa, "Strict" = col_strict, "IUCN IV" = col_iucn4,
  "IUCN V" = col_iucn5, "Natura 2000" = col_n2k
)


# Validation helpers --------------------------

assert_files_exist <- function(paths) {
  missing_paths <- paths[!file.exists(paths)]
  if (length(missing_paths) > 0) {
    stop("Required input file(s) are missing:\n", paste(missing_paths, collapse = "\n"))
  }
  invisible(TRUE)
}

# Check whether all required columns are available.
assert_columns <- function(data, columns, object_name = "data") {
  missing_columns <- setdiff(columns, names(data))
  if (length(missing_columns) > 0) {
    stop(object_name, " is missing required column(s): ", paste(missing_columns, collapse = ", "))
  }
  invisible(TRUE)
}


# Inside-versus-outside scenario helper ------------------------

# For category-specific comparisons:
#   inside  = plot lies in the focal PA category
#   outside = plot lies outside the complete PA network
# Plots that lie only in another PA category are excluded.

make_pa_scenario <- function(data, inside_flag) {
  data %>%
    dplyr::mutate(
      side = dplyr::case_when(
        .data[[inside_flag]] ~ "inside",
        !any_PA ~ "outside",
        TRUE ~ NA_character_
      ),
      side = factor(side, levels = c("outside", "inside"))
    ) %>%
    dplyr::filter(!is.na(side))
}

cat("ReSurvey preamble loaded.\n")
