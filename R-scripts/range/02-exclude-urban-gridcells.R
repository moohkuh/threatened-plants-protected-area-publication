# 02-exclude-urban-gridcells.R
# HILDA+ 1960 urban exclusion.
#
# PURPOSE
#   Calculate 1960 urban cover for exactly the MTBQ cells retained
#   by the ORIGINAL sMon centroid-in-Germany criterion from script 01.
#
# RULE USED IN THE PAPER
#   urban cell = > 50% urban cover in 1960.
#
# OUTPUT
#   ./data/derived/range_repro/grid_urban_1960.rds

suppressPackageStartupMessages({
  library(sf)
  library(terra)
  library(dplyr)
  library(readr)
})

sf::sf_use_s2(FALSE)

OUT_DIR <- "./data/derived/range_repro"
dir.create(OUT_DIR, recursive = TRUE, showWarnings = FALSE)

PATH_SPATIAL <- file.path(OUT_DIR, "range_spatial.gpkg")
PATH_HILDA <- paste0(
  "./data/spatial/hildap_vGLOB-1.0_geotiff/",
  "hildap_GLOB-v1.0_lulc-states/",
  "hilda_plus_1960_states_GLOB-v1-0_wgs84-nn.tif"
)

stopifnot(file.exists(PATH_SPATIAL))
stopifnot(file.exists(PATH_HILDA))

grid_cells <- sf::st_read(PATH_SPATIAL, layer = "grid_centroid_in_DE", quiet = TRUE) %>%
  sf::st_transform(25832)

# HILDA+ code used in the historical script:
# 11 = urban.

# We can build a cache file, as the raster operations take a lot of time
# Path for the cached, reprojected HILDA raster
PATH_HILDA_CACHE <- file.path(OUT_DIR, "hilda_1960_urban_utm25832_crop.tif")

if (file.exists(PATH_HILDA_CACHE)) {

  # Load the cached raster if it already exists
  hilda_crop <- terra::rast(PATH_HILDA_CACHE)

} else {

  # Load original HILDA raster
  hilda <- terra::rast(PATH_HILDA)

  # Reproject categorical raster using nearest-neighbour resampling
  hilda_utm <- terra::project(hilda, y = "EPSG:25832", method = "near", res = c(1000, 1000))

  # Crop to the analysis grid extent
  hilda_crop <- terra::crop(hilda_utm, terra::ext(grid_cells))

  # Save cached raster for future runs
  terra::writeRaster(hilda_crop, PATH_HILDA_CACHE, overwrite = TRUE)

  rm(hilda, hilda_utm)
}

urban_raster <- hilda_crop
urban_raster[urban_raster != 11] <- NA
rm(hilda_crop)

urban_sf <- terra::as.polygons(urban_raster, dissolve = TRUE) %>%
  sf::st_as_sf() %>%
  sf::st_make_valid() %>%
  sf::st_transform(sf::st_crs(grid_cells))

urban_area_per_cell <- suppressWarnings(
  sf::st_intersection(grid_cells %>% dplyr::select(MTB_Q, cell_area_m2), urban_sf)
) %>%
  dplyr::mutate(urban_area_m2 = as.numeric(sf::st_area(geom))) %>%
  sf::st_drop_geometry() %>%
  dplyr::group_by(MTB_Q) %>%
  dplyr::summarise(urban_area_m2 = sum(urban_area_m2, na.rm = TRUE), .groups = "drop")

# Join urban area back to all retained grid cells and calculate urban proportion
grid_urban <- grid_cells %>%
  sf::st_drop_geometry() %>%
  dplyr::select(MTB_Q, cell_area_m2) %>%
  dplyr::left_join(urban_area_per_cell, by = "MTB_Q") %>%
  dplyr::mutate(
    # Cells without an intersection with urban land have zero urban area
    urban_area_m2 = dplyr::coalesce(urban_area_m2, 0),
    # Proportion of each grid cell covered by urban land
    urban_prop = urban_area_m2 / cell_area_m2
  )

# Allow negligible numerical deviations from the theoretical 0-1 range
tolerance <- 1e-8

if (any(grid_urban$urban_prop < -tolerance | grid_urban$urban_prop > 1 + tolerance, na.rm = TRUE)) {
  stop("Urban proportions substantially outside [0, 1] were produced.")
}

# Clamp negligible floating-point deviations to the valid range
grid_urban <- grid_urban %>%
  dplyr::mutate(urban_prop = pmin(1, pmax(0, urban_prop)))

# Classify cells according to the threshold used in the paper
grid_urban <- grid_urban %>%
  dplyr::mutate(urban_class = dplyr::if_else(urban_prop > 0.5, "urban", "non-urban"))

saveRDS(grid_urban, file.path(OUT_DIR, "grid_urban_1960.rds"))
readr::write_csv(grid_urban, file.path(OUT_DIR, "grid_urban_1960.csv"))

qc <- tibble::tibble(
  metric = c("retained_cells_before_urban_filter", "urban_cells_gt_50pct", "nonurban_cells_final"),
  value = c(
    nrow(grid_urban),
    sum(grid_urban$urban_class == "urban"),
    sum(grid_urban$urban_class == "non-urban")
  )
)

readr::write_csv(qc, file.path(OUT_DIR, "QC_urban_summary.csv"))
print(qc)
