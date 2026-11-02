# 01-taxonomic-harmonisation.R
# Prepare frozen taxonomic lookup tables for the workflow.
#
# PURPOSE
#   1) Read manually checked taxonomic matches for:
#        - the German national Red List
#        - sMon occurrence-probability data
#        - ReSurvey Germany data
#   2) Standardize names and verify that each original taxon name maps
#      to only one WCVP name.
#   3) Attach national Red List information to the sMon and ReSurvey
#      taxonomic lookups.
#   4) Save clean lookup tables for the downstream analysis pipelines.
#
# IMPORTANT
#   Taxonomic matching to WCVP is NOT repeated here.
#   Fuzzy matches and manual corrections were reviewed during data
#   preparation and are treated as frozen inputs.
#
# RUN
#   Run this script from the root of the RStudio project.


# 0) Packages -------------------------------------------------------------

suppressPackageStartupMessages({
  library(dplyr)
  library(readr)
  library(stringr)
})


# 1) Paths ----------------------------------------------------------------

# Manually checked WCVP matches (included in the repository).
# The ReSurvey table only contains the name columns of the original
# matching table, without plot data or coordinates.
PATH_NAT_RL_CHECKED <- "./data/taxonomy/red_list_wcvp_matching.csv"
PATH_SMON_CHECKED   <- "./data/taxonomy/smon_wcvp_matching.csv"
PATH_RS_CHECKED     <- "./data/taxonomy/resurvey_wcvp_matching.csv"

# Clean derived taxonomy tables used by the analysis pipelines.
OUT_DIR <- "./data/derived/taxonomy"

OUT_NAT_RL_LOOKUP <- file.path(OUT_DIR, "national_red_list_lookup.csv")
OUT_SMON_LOOKUP   <- file.path(OUT_DIR, "smon_wcvp_lookup.csv")
OUT_SMON_MASTER   <- file.path(OUT_DIR, "smon_species_master.csv")
OUT_RS_LOOKUP     <- file.path(OUT_DIR, "resurvey_wcvp_lookup.csv")
OUT_RS_MASTER     <- file.path(OUT_DIR, "resurvey_species_master.csv")

dir.create(OUT_DIR, recursive = TRUE, showWarnings = FALSE)


# 2) Helper functions -----------------------------------------------------

# Stop early if a required input file is missing.
check_files <- function(paths) {
  missing_files <- paths[!file.exists(paths)]
  if (length(missing_files) > 0) {
    stop("Missing required taxonomy input file(s):\n", paste0("  - ", missing_files, collapse = "\n"))
  }
}

# Stop if required columns are missing from a table.
check_columns <- function(data, required, object_name) {
  missing_columns <- setdiff(required, names(data))
  if (length(missing_columns) > 0) {
    stop(object_name, " is missing required column(s): ", paste(missing_columns, collapse = ", "))
  }
}

# Standardize whitespace in taxon names without changing taxonomy.
clean_name <- function(x) {
  stringr::str_squish(as.character(x))
}

check_files(c(PATH_NAT_RL_CHECKED, PATH_SMON_CHECKED, PATH_RS_CHECKED))


# 3) National Red List lookup --------------------------------------------

# This is the frozen, manually checked national Red List -> WCVP table.
# Read the manually checked Red List taxonomy using its original encoding
nat_rl_raw <- readr::read_csv(
  PATH_NAT_RL_CHECKED,
  locale = readr::locale(encoding = "Windows-1252"),
  show_col_types = FALSE
)

check_columns(nat_rl_raw, c("wcvp_name", "threat_2018", "threat_1998", "native"), "National Red List taxonomy table")

nat_rl_raw <- nat_rl_raw %>%
  dplyr::mutate(
    wcvp_name = clean_name(wcvp_name),
    threat_2018 = as.character(threat_2018),
    threat_1998 = as.character(threat_1998)
  ) %>%
  dplyr::filter(!is.na(wcvp_name), wcvp_name != "")

# Check whether one WCVP name has conflicting Red List information.
nat_rl_conflicts <- nat_rl_raw %>%
  dplyr::group_by(wcvp_name) %>%
  dplyr::summarise(
    n_threat_2018 = dplyr::n_distinct(threat_2018, na.rm = TRUE),
    n_threat_1998 = dplyr::n_distinct(threat_1998, na.rm = TRUE),
    n_native = dplyr::n_distinct(native, na.rm = TRUE),
    .groups = "drop"
  ) %>%
  dplyr::filter(n_threat_2018 > 1 | n_threat_1998 > 1 | n_native > 1)

if (nrow(nat_rl_conflicts) > 0) {
  print(nat_rl_conflicts)
  stop("Conflicting national Red List information was found for the same WCVP name.")
}

# Keep the nationally red-listed categories used in the publication.
red_list_categories <- c("0", "1", "2", "3", "G", "R", "V")

nat_rl_lookup <- nat_rl_raw %>%
  dplyr::filter(threat_2018 %in% red_list_categories) %>%
  dplyr::transmute(wcvp_name, is_nat_rl = TRUE, threat_2018, threat_1998, native) %>%
  dplyr::distinct(wcvp_name, .keep_all = TRUE)

readr::write_csv(nat_rl_lookup, OUT_NAT_RL_LOOKUP)


# 4) sMon -> WCVP lookup --------------------------------------------------

# This file contains the frozen, manually checked sMon taxonomic matches.
smon_raw <- readr::read_csv(
  PATH_SMON_CHECKED,
  locale = readr::locale(encoding = "Windows-1252"),
  show_col_types = FALSE
)

check_columns(smon_raw, c("TaxonName", "wcvp_name"), "sMon taxonomy table")

smon_lookup <- smon_raw %>%
  dplyr::mutate(TaxonName = clean_name(TaxonName), wcvp_name = clean_name(wcvp_name)) %>%
  dplyr::filter(!is.na(TaxonName), TaxonName != "", !is.na(wcvp_name), wcvp_name != "") %>%
  dplyr::select(
    TaxonName,
    dplyr::any_of(c("TaxonName_cleaned", "wcvp_name", "wcvp_status", "wcvp_accepted_id"))
  ) %>%
  dplyr::distinct()

# Check that each original sMon taxon name maps to only one WCVP name.
smon_conflicts <- smon_lookup %>%
  dplyr::distinct(TaxonName, wcvp_name) %>%
  dplyr::count(TaxonName, name = "n_wcvp_names") %>%
  dplyr::filter(n_wcvp_names > 1)

if (nrow(smon_conflicts) > 0) {
  print(smon_conflicts)
  stop("Some sMon TaxonName values map to more than one WCVP name.")
}

smon_lookup <- smon_lookup %>%
  dplyr::distinct(TaxonName, .keep_all = TRUE)

readr::write_csv(smon_lookup, OUT_SMON_LOOKUP)

# Attach national Red List information to every harmonized sMon taxon.
smon_species_master <- smon_lookup %>%
  dplyr::left_join(nat_rl_lookup, by = "wcvp_name") %>%
  dplyr::mutate(is_nat_rl = dplyr::coalesce(is_nat_rl, FALSE))

readr::write_csv(smon_species_master, OUT_SMON_MASTER)


# 5) ReSurvey -> WCVP lookup ---------------------------------------------

# This file contains the WCVP matching table used for ReSurvey Germany.
rs_raw <- readr::read_csv(PATH_RS_CHECKED, show_col_types = FALSE)

check_columns(rs_raw, c("species", "wcvp_name", "wcvp_status", "wcvp_accepted_id"), "ReSurvey taxonomy table")

# Reproduce the match-priority rule used in the current ReSurvey pipeline.
# Accepted names and synonyms are preferred, followed by accepted names.
rs_lookup <- rs_raw %>%
  dplyr::mutate(species = clean_name(species), wcvp_name = clean_name(wcvp_name)) %>%
  dplyr::arrange(
    dplyr::desc(wcvp_status %in% c("Accepted", "Synonym")),
    dplyr::desc(wcvp_status == "Accepted")
  ) %>%
  dplyr::distinct(species, .keep_all = TRUE) %>%
  dplyr::select(species, wcvp_name, wcvp_status, wcvp_accepted_id)

# Check that original ReSurvey names are unique after match selection.
if (anyDuplicated(rs_lookup$species)) {
  stop("Duplicated ReSurvey species names remain after taxonomic harmonization.")
}

readr::write_csv(rs_lookup, OUT_RS_LOOKUP)

# Attach national Red List information to every harmonized ReSurvey taxon.
rs_species_master <- rs_lookup %>%
  dplyr::left_join(nat_rl_lookup, by = "wcvp_name") %>%
  dplyr::mutate(is_nat_rl = dplyr::coalesce(is_nat_rl, FALSE))

readr::write_csv(rs_species_master, OUT_RS_MASTER)


# 6) In-session QC summary -----------------------------------------------

taxonomy_qc <- tibble::tibble(
  dataset = c("National Red List", "sMon", "ReSurvey"),
  n_original_names = c(
    nrow(nat_rl_lookup),
    dplyr::n_distinct(smon_lookup$TaxonName),
    dplyr::n_distinct(rs_lookup$species)
  ),
  n_wcvp_names = c(
    dplyr::n_distinct(nat_rl_lookup$wcvp_name),
    dplyr::n_distinct(smon_lookup$wcvp_name),
    dplyr::n_distinct(rs_lookup$wcvp_name, na.rm = TRUE)
  ),
  n_national_red_list = c(
    nrow(nat_rl_lookup),
    sum(smon_species_master$is_nat_rl, na.rm = TRUE),
    sum(rs_species_master$is_nat_rl, na.rm = TRUE)
  )
)

print(taxonomy_qc)
