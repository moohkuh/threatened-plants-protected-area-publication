# 01-build_grid_and_pa_coverage.R
# Clean range-scale spatial preparation.
#
# PURPOSE
#   1) Reconstruct the MTBQ analysis grid from the grid_joined layer.
#   2) Use the ORIGINAL sMon OP coordinates to decide whether an MTBQ
#      centroid lies within Germany.
#   3) Build present-day protected-area coverage for the retained full
#      5 x 5 km cells.
#   4) Produce the five focal coverage variables used in the paper.
#
# DESIGN
#   A cell is retained when the sMon OP centroid for that MTBQ lies in
#   Germany. Border cells are NOT discarded merely because part of their
#   polygon extends across the national boundary.
#
#   PA coverage is calculated relative to the FULL MTBQ cell area
#
# OUTPUTS are written to ./data/derived/range_repro

suppressPackageStartupMessages({
  library(sf)
  library(dplyr)
  library(stringr)
  library(readr)
  library(purrr)
  library(tidyr)
})

sf::sf_use_s2(FALSE)


# 0) Paths -----------------------------------------------------------------

OUT_DIR <- "./data/derived/range_repro"
dir.create(OUT_DIR, recursive = TRUE, showWarnings = FALSE)

OP_FILES <- sprintf(
  "./data/raw/sMon/1875_9_1875_2_Modelled_OPs_incl_sd_pt_%d.csv", #%d is placeholder for an integer from 1 to 4, sprintf puts the numbers 1-4 instead of %f
  1:4
)

PATH_GRID_JOINED <- "./data/spatial/MTBQ/MTBQ_25832.shp"
PATH_DE   <- "./data/spatial/DE-boundaries/VG2500_STA.shp"
PATH_SMON_MASTER <- "./data/derived/taxonomy/smon_species_master.csv"
WDPA_FILES <- c(
  "./data/spatial/Protected_areas/WDPA_WDOECM_Jan2025_Public_DEU_shp/WDPA_WDOECM_Jan2025_Public_DEU_shp_0/WDPA_WDOECM_Jan2025_Public_DEU_shp-polygons.shp",
  "./data/spatial/Protected_areas/WDPA_WDOECM_Jan2025_Public_DEU_shp/WDPA_WDOECM_Jan2025_Public_DEU_shp_1/WDPA_WDOECM_Jan2025_Public_DEU_shp-polygons.shp",
  "./data/spatial/Protected_areas/WDPA_WDOECM_Jan2025_Public_DEU_shp/WDPA_WDOECM_Jan2025_Public_DEU_shp_2/WDPA_WDOECM_Jan2025_Public_DEU_shp-polygons.shp"
)

# Check that all required input files are available before starting
required_files <- c(OP_FILES, PATH_GRID_JOINED, PATH_DE, WDPA_FILES)

if (!all(file.exists(required_files))) {
  missing_files <- required_files[!file.exists(required_files)]
  stop("The following required input files are missing:\n", paste(missing_files, collapse = "\n"))
}


# 1) Reconstruct one original OP centroid per MTBQ -------------------------
# The four raw sMon files repeat MTBQ coordinates for many species.
# Small numerical differences between repeated coordinates are possible.
# We therefore retain all distinct coordinate candidates per MTBQ

centroid_candidates <- purrr::map_dfr(OP_FILES, function(f) { # map_dfr() reads each file and combines the results into a single data frame - it works similar as if we just used bind_rows
  readr::read_csv(
    f, # this is the file path to read
    show_col_types = FALSE, # this just suppresses the printing of column types when reading the CSV
    col_select = c(MTB_Q, Longitude, Latitude) # this selects only the columns MTB_Q, Longitude, and Latitude from the CSV
  )
}) %>% # now that we have our very long dataframe with MTBQs, Lat and Longs, we will check if there are NAs and if the coordinates are finite numbers, then we will round the coordinates to 5 decimal places and keep only distinct row
  dplyr::filter(
    !is.na(MTB_Q), #, we check if there are NAs
    is.finite(Longitude), #if the coordinates are finite numbers
    is.finite(Latitude)
  ) %>%
  dplyr::transmute( # this transmutes the columns to the desired types
    MTB_Q = as.character(MTB_Q),
    Longitude = as.numeric(Longitude),
    Latitude  = as.numeric(Latitude)
  ) %>%
  dplyr::mutate( # we round the coordinates to 5 decimal places
    Longitude = round(Longitude, 5),
    Latitude  = round(Latitude, 5)
  ) %>%
  dplyr::distinct(MTB_Q, Longitude, Latitude) # and keep only distinct row
# This is because the distinct coordinates per grid cell have several hits otherwise due to inaccuracy in the least significant digits of the coordinates

# Check that rounding leaves exactly one coordinate pair per MTBQ
centroid_check <- centroid_candidates %>%
  dplyr::count(MTB_Q, name = "n_coordinate_pairs")


# 2) Germany polygon + validation of repeated centroid coordinates ---------

de_poly <- sf::st_read(PATH_DE, quiet = TRUE) %>% # to check which unique MTBQs are located within the German border, we load our Germany polygon
  dplyr::mutate(BEGINN = as.Date(BEGINN)) %>%  # the BEGINN column is a date, so we convert it to a Date object
  dplyr::arrange(dplyr::desc(BEGINN)) %>% # we arrange the rows in descending order of the BEGINN column, so that the most recent polygons are first
  dplyr::slice(1) %>% # we take the first row, which is the most recent polygon
  sf::st_transform(25832) %>% # we transform the coordinate reference system to EPSG:25832, which is a projected coordinate system for Germany
  sf::st_make_valid() %>%# we make the geometry valid, in case there are any invalid geometries
  dplyr::summarise(geometry = sf::st_union(geometry)) # we unionize all geometries (in case there are multiple)

candidate_sf <- sf::st_as_sf(
  centroid_candidates,
  coords = c("Longitude", "Latitude"),
  crs = 4326,
  remove = FALSE
) %>%
  sf::st_transform(25832) # we transform our grid centroids to the same CRS as the germany polygon

candidate_tbl <- candidate_sf %>%
  sf::st_drop_geometry() %>%
  dplyr::mutate(
    in_DE = lengths(sf::st_intersects(candidate_sf, de_poly)) > 0 # this checks if each point intersects with the Germany polygon, and returns TRUE if it does, FALSE otherwise
  )

# Save the centroid table used for the Germany-membership decision
readr::write_csv(candidate_tbl, file.path(OUT_DIR, "MTBQ_centroid_DE_status.csv"))

valid_MTBQ <- candidate_tbl %>%
  dplyr::filter(in_DE) %>%
  dplyr::pull(MTB_Q) %>%
  unique()

saveRDS(valid_MTBQ, file.path(OUT_DIR, "valid_MTBQ_centroid_in_DE.rds"))


# 3) Read the historical MTBQ polygon layer ---------------------------------
# and retain cells by the ORIGINAL OP-centroid Germany criterion
#
#   We use the polygon geometry only for PA-area calculations.
#   Whether a cell belongs to the German analysis universe is determined
#   by the ORIGINAL OP centroid from Sections 1-2 above, NOT by
#   st_point_on_surface() of the polygon.

smon_master <- readr::read_csv(PATH_SMON_MASTER, show_col_types = FALSE) %>%
  dplyr::mutate(TaxonName = stringr::str_squish(as.character(TaxonName)))

if (!all(c("TaxonName", "is_nat_rl") %in% names(smon_master))) {
  stop("smon_species_master.csv must contain TaxonName and is_nat_rl.")
}

smon_rl_key <- smon_master %>%
  dplyr::filter(is_nat_rl) %>%
  dplyr::select(TaxonName) %>%
  dplyr::distinct()

# Recreate the set of MTBQ centroids used for the historical grid_joined layer.
analysis_centroids <- purrr::map_dfr(OP_FILES, function(f) {
  readr::read_csv(f, show_col_types = FALSE, col_select = c(MTB_Q, TaxonName, Longitude, Latitude, Timestep))
}) %>%
  dplyr::filter(Timestep %in% c(1L, 2L, 3L)) %>%
  dplyr::mutate(
    MTB_Q = as.character(MTB_Q),
    TaxonName = stringr::str_squish(as.character(TaxonName)),
    Longitude = as.numeric(Longitude),
    Latitude = as.numeric(Latitude)
  ) %>%
  dplyr::inner_join(smon_rl_key, by = "TaxonName") %>%
  dplyr::filter(!is.na(MTB_Q), is.finite(Longitude), is.finite(Latitude)) %>%
  dplyr::mutate(Longitude = round(Longitude, 5), Latitude = round(Latitude, 5)) %>%
  dplyr::select(MTB_Q, Longitude, Latitude) %>%
  dplyr::distinct()

# The historical construction assumes one centroid per MTBQ.
analysis_centroid_check <- analysis_centroids %>%
  dplyr::count(MTB_Q, name = "n_coordinate_pairs") %>%
  dplyr::filter(n_coordinate_pairs != 1)

if (nrow(analysis_centroid_check) > 0) {
  print(analysis_centroid_check)
  stop("The red-listed sMon data contain more than one centroid for some MTBQ cells.")
}

# Convert the analysis centroids to spatial points for matching them
# to the historical MTBQ polygons.
analysis_centroids_sf <- analysis_centroids %>%
  dplyr::mutate(id = dplyr::row_number()) %>%
  sf::st_as_sf(coords = c("Longitude", "Latitude"), crs = 4326, remove = FALSE) %>%
  sf::st_transform(25832)

# Load the original MTBQ polygon grid.
# Only the geometry is needed; all original attribute columns are removed.
grid_base <- sf::st_read(PATH_GRID_JOINED, quiet = TRUE, options = "ENCODING=LATIN1") %>%
  sf::st_transform(25832) %>%
  sf::st_make_valid()

# Remove all original attribute columns and keep only the geometry.
grid_base <- grid_base %>%
  dplyr::select(geometry)

# Spatially attach the sMon MTBQ identifier to the polygon
# containing the corresponding occurrence-probability centroid.
grid_raw <- sf::st_join(
  grid_base,
  analysis_centroids_sf %>% dplyr::select(id, MTB_Q),
  join = sf::st_contains,
  left = FALSE
) %>%
  dplyr::mutate(MTB_Q = as.character(MTB_Q)) %>%
  dplyr::select(id, MTB_Q, geometry)

# Check whether any MTBQ identifier occurs more than once.
duplicated_MTBQ <- grid_raw %>%
  sf::st_drop_geometry() %>%
  dplyr::count(MTB_Q, name = "n_polygons") %>%
  dplyr::filter(n_polygons > 1)

if (nrow(duplicated_MTBQ) > 0) {
  print(duplicated_MTBQ)
  stop("Some MTBQ identifiers were assigned to more than one grid polygon.")
}

# Check whether all analysis centroids received a grid polygon.
unmatched_MTBQ <- setdiff(analysis_centroids$MTB_Q, grid_raw$MTB_Q)

if (length(unmatched_MTBQ) > 0) {
  stop("Some MTBQ centroids could not be matched to a grid polygon: ", paste(unmatched_MTBQ, collapse = ", "))
}

# Retain only cells whose original sMon centroid lies within Germany.
grid_cells <- grid_raw %>%
  dplyr::filter(MTB_Q %in% valid_MTBQ) %>%
  dplyr::mutate(cell_area_m2 = as.numeric(sf::st_area(geometry)), cell_area_km2 = cell_area_m2 / 1e6) %>%
  dplyr::select(id, MTB_Q, cell_area_m2, cell_area_km2, geometry)

# Save a clean grid.
sf::st_write(grid_cells, file.path(OUT_DIR, "range_spatial.gpkg"), layer = "grid_centroid_in_DE", delete_dsn = TRUE, quiet = TRUE)


# 4) Load and clean WDPA polygons -------------------------------------------

pa_parts <- purrr::map(WDPA_FILES, ~ sf::st_read(.x, quiet = TRUE))  # reads and combines the WDPA shapefiles to a list

# remove marine and one known WH polygon.
pa_parts[[1]] <- pa_parts[[1]] %>% dplyr::filter(MARINE != "2")
pa_parts[[2]] <- pa_parts[[2]] %>% dplyr::filter(MARINE != "2", WDPAID != 903141)
pa_parts[[3]] <- pa_parts[[3]] %>% dplyr::filter(MARINE != "2")

PAs <- dplyr::bind_rows(pa_parts) %>%
  sf::st_transform(25832) %>%
  sf::st_make_valid() %>%
  dplyr::mutate(
    STATUS_YR_int = suppressWarnings(as.integer(STATUS_YR)),
    desig_text = stringr::str_to_lower(paste(
      dplyr::coalesce(as.character(DESIG), ""),
      dplyr::coalesce(as.character(DESIG_ENG), "")
    )),
    is_national_natural_monument = stringr::str_detect(
      desig_text, "national.*naturmonument|national natural monument"
    ),
    is_biosphere_core = stringr::str_detect(
      desig_text, "biosphärenreservat\\s*-\\s*kernzone|biosphere reserve.*core"
    ),
    is_habitats = DESIG %in% c(
      "Special Areas of Conservation (Habitats Directive)",
      "Site of Community Importance (Habitats Directive)"
    ),
    is_birds = DESIG == "Special Protection Area (Birds Directive)"
  )

# Category assignment used downstream.
# National Natural Monuments are deliberately grouped with strict protection.
# Polygon-only WDPA input means point natural monuments are not included.
PAs <- PAs %>%
  dplyr::mutate(
    IUCN_CAT_numeric = dplyr::case_when(
      is_national_natural_monument ~ 4L,
      IUCN_CAT %in% c("Ia", "Ib", "Ia/Ib") ~ 1L,
      IUCN_CAT == "II" ~ 2L,
      is_biosphere_core ~ 3L,
      IUCN_CAT == "III" & stringr::str_detect(desig_text, "landschafts?schutz") ~ 6L,
      IUCN_CAT == "III" ~ 4L,
      IUCN_CAT == "IV" ~ 5L,
      IUCN_CAT == "V" ~ 6L,
      is_habitats ~ 7L,
      is_birds ~ 8L,
      DESIG == "Ramsar Site, Wetland of International Importance" ~ 9L,
      DESIG %in% c("Biosphärenreservat - Pflegezone", "Biosphärenreservat - Entwicklungszone") ~ 10L,
      IUCN_CAT %in% c("Not Applicable", "Not Assigned", "Not Reported") ~ 11L,
      TRUE ~ 12L
    ),
    IUCN_label = dplyr::case_when(
      IUCN_CAT_numeric == 1L  ~ "Ia/Ib",
      IUCN_CAT_numeric == 2L  ~ "II",
      IUCN_CAT_numeric == 3L  ~ "BR core",
      IUCN_CAT_numeric == 4L  ~ "III / National Natural Monument",
      IUCN_CAT_numeric == 5L  ~ "IV",
      IUCN_CAT_numeric == 6L  ~ "V",
      IUCN_CAT_numeric == 7L  ~ "N2000 Habitats",
      IUCN_CAT_numeric == 8L  ~ "N2000 Birds",
      IUCN_CAT_numeric == 9L  ~ "Ramsar",
      IUCN_CAT_numeric == 10L ~ "BR buffer/transition",
      IUCN_CAT_numeric == 11L ~ "Other",
      TRUE ~ "Mixed/Undefined"
    )
  ) %>%
  dplyr::select(-desig_text)

# Remove the large parent Biosphere Reserve hulls before area calculations
# because their separate zones are already represented.
PAs_for_intersection <- PAs %>%
  dplyr::filter(!(DESIG %in% c("Biosphärenreservat", "UNESCO-MAB Biosphere Reserve")))

# Save the clean PA layer used for all intersections.
sf::st_write(PAs_for_intersection, file.path(OUT_DIR, "range_spatial.gpkg"), layer = "protected_areas_clean", append = TRUE, quiet = TRUE)


# 5) Grid x PA intersection ---------------------------------------------------
# this takes a while

intersected <- suppressWarnings(sf::st_intersection(grid_cells, PAs_for_intersection)) %>%
  sf::st_make_valid()

sf::st_write(intersected, file.path(OUT_DIR, "range_spatial.gpkg"), layer = "grid_pa_intersection", append = TRUE, quiet = TRUE)


# 6) Coverage helper ----------------------------------------------------------

coverage_from_subset <- function(intersection_sf, grid_sf, keep, out_name) { #outname is the name of the newly created column (e.g. cov_V)

  # Keep only the selected protected-area subset
  filtered_pa <- intersection_sf |> dplyr::filter({{ keep }})

  # Calculate protected area per grid cell
  area_table <- filtered_pa |>
    dplyr::group_by(MTB_Q) |>
    dplyr::summarise(geometry = sf::st_union(geometry), .groups = "drop") |>
    sf::st_make_valid() |>
    dplyr::mutate(group_area_m2 = as.numeric(sf::st_area(geometry))) |>
    sf::st_drop_geometry() |>
    dplyr::select(MTB_Q, group_area_m2)

  # Join protected area back to all grid cells and calculate coverage
  result <- grid_sf |>
    sf::st_drop_geometry() |>
    dplyr::select(MTB_Q, cell_area_m2) |>
    dplyr::left_join(area_table, by = "MTB_Q") |>
    dplyr::mutate(
      group_area_m2 = dplyr::coalesce(group_area_m2, 0),
      !!out_name := pmin( # gives the smaller value between 1 and x (in this case pmax(0,...) which gives the bigger value between 0 and the ratio)
        # and therefore is gives values between 0 and 1
        # !! says "create a column with the name saved in out_name
        1,
        pmax(0, group_area_m2 / cell_area_m2)
      )
    ) |>
    dplyr::select(MTB_Q, !!out_name)

  result
}

# All PAs: union ALL retained terrestrial PA polygons per cell.
cov_any <- intersected %>%
  dplyr::group_by(MTB_Q) %>%
  dplyr::summarise(geometry = sf::st_union(geometry), .groups = "drop") %>%
  sf::st_make_valid() %>%
  dplyr::mutate(group_area_m2 = as.numeric(sf::st_area(geometry))) %>%
  sf::st_drop_geometry() %>%
  dplyr::select(MTB_Q, group_area_m2)

cov_any <- grid_cells %>%
  sf::st_drop_geometry() %>%
  dplyr::select(MTB_Q, cell_area_m2) %>%
  dplyr::left_join(cov_any, by = "MTB_Q") %>%
  dplyr::mutate(
    group_area_m2 = dplyr::coalesce(group_area_m2, 0),
    cov_anyPA = pmin(1, pmax(0, group_area_m2 / cell_area_m2))
  ) %>%
  dplyr::select(MTB_Q, cov_anyPA)

cov_strict <- coverage_from_subset(intersected, grid_cells, IUCN_CAT_numeric %in% c(1L, 2L, 3L, 4L), "cov_strict")
cov_IV <- coverage_from_subset(intersected, grid_cells, IUCN_CAT_numeric == 5L, "cov_IV")
cov_V <- coverage_from_subset(intersected, grid_cells, IUCN_CAT_numeric == 6L, "cov_V")
cov_N2k <- coverage_from_subset(intersected, grid_cells, IUCN_CAT_numeric == 7L, "cov_n2000_habitats_only")

grid_cov <- grid_cells %>%
  sf::st_drop_geometry() %>%
  dplyr::select(MTB_Q, cell_area_m2, cell_area_km2) %>%
  dplyr::left_join(cov_any,    by = "MTB_Q") %>%
  dplyr::left_join(cov_strict, by = "MTB_Q") %>%
  dplyr::left_join(cov_IV,     by = "MTB_Q") %>%
  dplyr::left_join(cov_V,      by = "MTB_Q") %>%
  dplyr::left_join(cov_N2k,    by = "MTB_Q")

coverage_cols <- c("cov_anyPA", "cov_strict", "cov_IV", "cov_V", "cov_n2000_habitats_only")

# Check that all final PA coverage values are finite and on the expected 0-1 scale
coverage_matrix <- as.matrix(grid_cov[, coverage_cols])

if (any(!is.finite(coverage_matrix))) {
  stop("Non-finite values were produced in the final PA coverage table.")
}

if (any(coverage_matrix < 0 | coverage_matrix > 1)) {
  stop("PA coverage values outside the expected 0-1 range were produced.")
}

saveRDS(grid_cov, file.path(OUT_DIR, "grid_PA_coverage_present.rds"))
readr::write_csv(grid_cov, file.path(OUT_DIR, "grid_PA_coverage_present.csv"))

# 7) Optional time flag for PA designation year (not used by the main range models, which use present-day PA coverage) ---

intersected_time <- intersected %>%
  dplyr::mutate(
    time_bin = dplyr::case_when(
      !is.na(STATUS_YR_int) & STATUS_YR_int <= 1987L ~ "T1_or_earlier",
      is.na(STATUS_YR_int) ~ "year_unknown",
      TRUE ~ "after_1987"
    )
  )

sf::st_write(intersected_time, file.path(OUT_DIR, "range_spatial.gpkg"), layer = "grid_pa_intersection_with_time", append = TRUE, quiet = TRUE)


# 8) Publication tables: areas and overlap among focal PA groups -----------------

focal_specs <- list(
  "Strict"      = c(1L, 2L, 3L, 4L),
  "IUCN IV"     = 5L,
  "IUCN V"      = 6L,
  "Natura 2000" = 7L
)

analysis_area_km2 <- sum(grid_cells$cell_area_km2, na.rm = TRUE)

focal_geom <- purrr::imap(focal_specs, function(ids, nm) {
  z <- intersected %>%
    dplyr::filter(IUCN_CAT_numeric %in% ids) %>%
    dplyr::summarise(geometry = sf::st_union(geometry)) %>%
    sf::st_make_valid()
  # st_make_valid() guarantees a topologically valid geometry, but after
  # unioning many small intersection polygons it can still leave coordinate
  # precision noise (near-duplicate vertices) that trips GEOS's noding
  # algorithm in a later st_intersection() call ("found non-noded
  # intersection..."). A zero-width buffer forces GEOS to re-node the
  # geometry and removes that noise; it does not change the area.
  z <- sf::st_buffer(z, 0)
  z$PA_category <- nm
  z
})

focal_areas <- purrr::map_dfr(focal_geom, function(z) {
  tibble::tibble(PA_category = z$PA_category, area_km2 = as.numeric(sf::st_area(z)) / 1e6)
}) %>%
  dplyr::mutate(pct_analysis_area = 100 * area_km2 / analysis_area_km2)

pairwise_overlap <- purrr::map_dfr(names(focal_geom), function(a) {
  purrr::map_dfr(setdiff(names(focal_geom), a), function(b) {
    ga <- focal_geom[[a]]
    gb <- focal_geom[[b]]
    shared <- suppressWarnings(sf::st_intersection(ga, gb))
    shared_km2 <- if (nrow(shared) == 0) 0 else as.numeric(sf::st_area(shared)) / 1e6
    a_km2 <- as.numeric(sf::st_area(ga)) / 1e6

    tibble::tibble(
      focal = a,
      overlaps_with = b,
      focal_area_km2 = a_km2,
      shared_area_km2 = shared_km2,
      pct_of_focal = 100 * shared_km2 / a_km2
    )
  })
})

# Table S4 "% total shared area": for each focal category, the share of its
# area that overlaps with the UNION of all other focal categories (not the
# sum of the pairwise overlaps above, which would double-count area shared
# with more than one other category at once).
total_shared <- purrr::map_dfr(names(focal_geom), function(a) {
  message("total_shared: processing ", a, " ...")
  
  ga <- focal_geom[[a]]
  a_km2 <- as.numeric(sf::st_area(ga)) / 1e6
  
  overlaps <- purrr::map(setdiff(names(focal_geom), a), function(b) {
    suppressWarnings(sf::st_intersection(ga, focal_geom[[b]]))
  })
  overlaps <- overlaps[purrr::map_lgl(overlaps, ~ nrow(.x) > 0)]
  
  if (length(overlaps) == 0) {
    shared_km2 <- 0
  } else {
    shared <- overlaps %>%
      dplyr::bind_rows() %>%
      sf::st_union() %>%
      sf::st_make_valid() %>%
      sf::st_buffer(0)
    shared_km2 <- as.numeric(sf::st_area(shared)) / 1e6
  }
  
  out <- tibble::tibble(
    PA_category = a,
    focal_area_km2 = a_km2,
    total_shared_area_km2 = shared_km2,
    pct_total_shared_area = 100 * shared_km2 / a_km2
  )
  
  # checkpoint, so a later crash doesn't lose this category's result
  saveRDS(out, file.path("diagnostics", paste0("total_shared_", a, ".rds")))
  message("total_shared: done with ", a)
  
  out
})

focal_areas <- focal_areas %>%
  dplyr::left_join(
    total_shared %>% dplyr::select(PA_category, pct_total_shared_area),
    by = "PA_category"
  )

readr::write_csv(focal_areas, file.path(OUT_DIR, "QC_PA_area_by_focal_category.csv"))
readr::write_csv(pairwise_overlap, file.path(OUT_DIR, "QC_PA_pairwise_overlap.csv"))
readr::write_csv(total_shared, file.path(OUT_DIR, "QC_PA_total_shared_area.csv"))


# 9) Final summary ---------------------------------------------------------

qc <- tibble::tibble(
  metric = c("unique_OP_MTBQ_centroids", "centroids_in_Germany", "retained_MTBQ_polygons", "analysis_area_km2"),
  value = c(nrow(centroid_candidates), length(valid_MTBQ), nrow(grid_cells), analysis_area_km2)
)

print(qc)
