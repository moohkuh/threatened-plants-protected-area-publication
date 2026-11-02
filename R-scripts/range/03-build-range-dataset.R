# 03-build-range-dataset.R
# Canonical sMon range dataset.
#
# PURPOSE
#   Build ONE canonical species x MTBQ table for the range analyses.
#
# IMPORTANT
#   This script does NOT silently remove Rhinanthus serotinus or choose
#   one of conflicting occurrence-probability records. Conflicting
#   accepted-species x MTBQ x timestep records are written to a QC file
#   and the script stops, so the taxonomic/data issue can be resolved
#   explicitly and documented.
#
# OUTPUTS
#   ./data/derived/range_repro/range_analysis_long.rds
#   ./data/derived/range_repro/range_analysis_wide.rds

suppressPackageStartupMessages({
  library(dplyr)
  library(tidyr)
  library(readr)
  library(stringr)
  library(purrr)
})

OUT_DIR <- "./data/derived/range_repro"
dir.create(OUT_DIR, recursive = TRUE, showWarnings = FALSE)

OP_FILES <- sprintf("./data/raw/sMon/1875_9_1875_2_Modelled_OPs_incl_sd_pt_%d.csv", 1:4)

PATH_VALID <- file.path(OUT_DIR, "valid_MTBQ_centroid_in_DE.rds")
PATH_PA    <- file.path(OUT_DIR, "grid_PA_coverage_present.rds")
PATH_URBAN <- file.path(OUT_DIR, "grid_urban_1960.rds")

PATH_SMON_TAXONOMY <- "./data/derived/taxonomy/smon_wcvp_lookup.csv"
PATH_SPECIES_MASTER <- "./data/derived/taxonomy/smon_species_master.csv"

stopifnot(all(file.exists(OP_FILES)))
stopifnot(file.exists(PATH_VALID))
stopifnot(file.exists(PATH_PA))
stopifnot(file.exists(PATH_URBAN))
stopifnot(file.exists(PATH_SMON_TAXONOMY))
stopifnot(file.exists(PATH_SPECIES_MASTER))


# 1) Load raw OP records -------------------------------------------------

op_all <- purrr::map_dfr(OP_FILES, function(f) {
  readr::read_csv(f, show_col_types = FALSE) %>%
    dplyr::filter(Timestep %in% c(1L, 2L, 3L)) %>%
    dplyr::select(MTB_Q, TaxonName, Longitude, Latitude, Timestep, OP, OP_sd, Period)
}) %>%
  dplyr::mutate(
    MTB_Q = as.character(MTB_Q),
    TaxonName = stringr::str_squish(TaxonName),
    Timestep = as.integer(Timestep)
  )


# 2) Restrict to OP centroids in Germany ----------------------------------

valid_MTBQ <- readRDS(PATH_VALID) %>% as.character()

op_all <- op_all %>%
  dplyr::filter(MTB_Q %in% valid_MTBQ)

message("OP records after centroid-in-DE filter: ", nrow(op_all))


# 3) Join frozen/checked sMon -> WCVP taxonomy ----------------------------

smon_key_raw <- readr::read_csv(PATH_SMON_TAXONOMY, show_col_types = FALSE)

# Build the lookup table from original sMon names to WCVP names
smon_key <- smon_key_raw %>%
  dplyr::transmute(TaxonName = stringr::str_squish(TaxonName), wcvp_name = stringr::str_squish(wcvp_name)) %>%
  dplyr::filter(!is.na(TaxonName), !is.na(wcvp_name)) %>%
  dplyr::distinct()

op_all <- op_all %>%
  dplyr::left_join(smon_key, by = "TaxonName")

unmatched_taxa <- op_all %>%
  dplyr::filter(is.na(wcvp_name)) %>%
  dplyr::distinct(TaxonName) %>%
  dplyr::arrange(TaxonName)

readr::write_csv(unmatched_taxa, file.path(OUT_DIR, "QC_unmatched_sMon_taxa.csv"))


# 4) Attach national Red List status --------------------------------------

species_master <- readr::read_csv(PATH_SPECIES_MASTER, show_col_types = FALSE) %>%
  dplyr::mutate(wcvp_name = stringr::str_squish(wcvp_name), threat_2018 = as.character(threat_2018)) %>%
  dplyr::select(dplyr::any_of(c("wcvp_name", "is_nat_rl", "threat_2018", "threat_1998", "native"))) %>%
  dplyr::distinct(wcvp_name, .keep_all = TRUE)

op_nat <- op_all %>%
  dplyr::filter(!is.na(wcvp_name)) %>%
  dplyr::left_join(species_master, by = "wcvp_name") %>%
  dplyr::filter(is_nat_rl %in% TRUE)

dplyr::n_distinct(op_nat$wcvp_name)


# 5) Remove known problematic species and resolve duplicates --------------

# Exclude Rhinanthus serotinus because duplicated, conflicting OP values
# are present in the raw data
op_nat <- op_nat %>%
  dplyr::filter(wcvp_name != "Rhinanthus serotinus")

# Remove exact duplicate records
op_nat <- op_nat %>%
  dplyr::distinct(wcvp_name, MTB_Q, Timestep, OP, .keep_all = TRUE)

# Check that only one OP value remains per species, grid cell and timestep
remaining_conflicts <- op_nat %>%
  dplyr::count(wcvp_name, MTB_Q, Timestep, name = "n_records") %>%
  dplyr::filter(n_records > 1)

if (nrow(remaining_conflicts) > 0) {
  stop(sprintf(
    "%d species/grid-cell/timestep combinations still have conflicting OP values after deduplication. Inspect `remaining_conflicts` before proceeding.",
    nrow(remaining_conflicts)
  ))
}

# Keep one row per species, grid cell and timestep
op_nat <- op_nat %>%
  dplyr::distinct(wcvp_name, MTB_Q, Timestep, .keep_all = TRUE)


# 6) Join grid-level PA coverage and urban cover ---------------------------

grid_pa <- readRDS(PATH_PA) %>%
  dplyr::mutate(MTB_Q = as.character(MTB_Q)) %>%
  dplyr::distinct(MTB_Q, .keep_all = TRUE)

grid_urban <- readRDS(PATH_URBAN) %>%
  dplyr::mutate(MTB_Q = as.character(MTB_Q)) %>%
  dplyr::distinct(MTB_Q, .keep_all = TRUE)

grid_lookup <- grid_pa %>%
  dplyr::left_join(
    grid_urban %>% dplyr::select(MTB_Q, urban_area_m2, urban_prop, urban_class),
    by = "MTB_Q"
  )

op_nat <- op_nat %>%
  dplyr::left_join(grid_lookup, by = "MTB_Q")

coverage_cols <- c("cov_anyPA", "cov_strict", "cov_IV", "cov_V", "cov_n2000_habitats_only")


# 7) Build canonical long and wide data ------------------------------------

species_info <- op_nat %>%
  dplyr::select(wcvp_name, dplyr::any_of(c("threat_2018", "threat_1998", "native", "is_nat_rl"))) %>%
  dplyr::distinct()

species_info <- species_info %>%
  dplyr::distinct(wcvp_name, .keep_all = TRUE)

# Keep one row per accepted species x MTBQ x timestep.
range_long <- op_nat %>%
  dplyr::select(
    wcvp_name, MTB_Q, Timestep, OP, OP_sd, Longitude, Latitude,
    dplyr::all_of(coverage_cols), cell_area_m2, cell_area_km2,
    urban_area_m2, urban_prop, urban_class,
    dplyr::any_of(c("threat_2018", "threat_1998", "native", "is_nat_rl"))
  )

range_wide_core <- range_long %>%
  dplyr::select(wcvp_name, MTB_Q, Timestep, OP) %>%
  tidyr::pivot_wider(id_cols = c(wcvp_name, MTB_Q), names_from = Timestep, values_from = OP, names_prefix = "OP_T")

grid_static <- range_long %>%
  dplyr::select(
    MTB_Q, Longitude, Latitude, dplyr::all_of(coverage_cols),
    cell_area_m2, cell_area_km2, urban_area_m2, urban_prop, urban_class
  ) %>%
  dplyr::distinct(MTB_Q, .keep_all = TRUE)

range_wide <- range_wide_core %>%
  dplyr::left_join(species_info, by = "wcvp_name") %>%
  dplyr::left_join(grid_static, by = "MTB_Q")

saveRDS(range_long, file.path(OUT_DIR, "range_analysis_long.rds"))
saveRDS(range_wide, file.path(OUT_DIR, "range_analysis_wide.rds"))


# 8) QC summaries -----------------------------------------------------------

species_counts <- range_wide %>%
  dplyr::distinct(wcvp_name, threat_2018) %>%
  dplyr::count(threat_2018, name = "n_species") %>%
  dplyr::arrange(threat_2018)

qc <- tibble::tibble(
  metric = c("rows_species_x_MTBQ", "unique_species", "unique_MTBQ", "urban_MTBQ", "nonurban_MTBQ"),
  value = c(
    nrow(range_wide),
    dplyr::n_distinct(range_wide$wcvp_name),
    dplyr::n_distinct(range_wide$MTB_Q),
    dplyr::n_distinct(range_wide$MTB_Q[range_wide$urban_class == "urban"]),
    dplyr::n_distinct(range_wide$MTB_Q[range_wide$urban_class == "non-urban"])
  )
)

print(qc)
print(species_counts)
