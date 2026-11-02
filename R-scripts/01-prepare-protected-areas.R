# 01-prepare-protected-areas.R
# Shared protected-area preparation for all analyses.
#
# PURPOSE
#   Prepare one common protected-area polygon layer that is used
#   by both the range-scale and ReSurvey analyses.
#
# INPUT
#   WDPA / Protected Planet polygon data for Germany,
#   downloaded in January 2025.
#
# OUTPUT
#   data/derived/protected_areas/PA_polygons_slim.gpkg
#
# IMPORTANT
#   The protected-area classification in this script defines the
#   categories used throughout the entire publication.


# 0) Setup ---------------------------------------------------------------

suppressPackageStartupMessages({
  library(sf)
  library(dplyr)
  library(stringr)
  library(purrr)
  library(readr)
})

sf::sf_use_s2(FALSE)


# 1) Paths ---------------------------------------------------------------

WDPA_FILES <- c(
  "./data/spatial/Protected_areas/WDPA_WDOECM_Jan2025_Public_DEU_shp/WDPA_WDOECM_Jan2025_Public_DEU_shp_0/WDPA_WDOECM_Jan2025_Public_DEU_shp-polygons.shp",
  "./data/spatial/Protected_areas/WDPA_WDOECM_Jan2025_Public_DEU_shp/WDPA_WDOECM_Jan2025_Public_DEU_shp_1/WDPA_WDOECM_Jan2025_Public_DEU_shp-polygons.shp",
  "./data/spatial/Protected_areas/WDPA_WDOECM_Jan2025_Public_DEU_shp/WDPA_WDOECM_Jan2025_Public_DEU_shp_2/WDPA_WDOECM_Jan2025_Public_DEU_shp-polygons.shp"
)

# Shared protected-area output
OUT_PA_DIR <- "./data/derived/protected_areas"
dir.create(OUT_PA_DIR, recursive = TRUE, showWarnings = FALSE)

PATH_PA_GPKG <- file.path(OUT_PA_DIR, "PA_polygons_slim.gpkg")
LAYER_PA <- "pa_polygons_slim"


# 2) Read WDPA polygon files ---------------------------------------------

pa_parts <- purrr::map(WDPA_FILES, ~ sf::st_read(.x, quiet = TRUE))


# 3) Remove fully marine polygons ---------------------------------------

# Keep terrestrial and coastal protected areas, but remove
# polygons classified as fully marine.
pa_parts <- purrr::map(pa_parts, ~ .x %>% dplyr::filter(MARINE != "2"))

# Remove the known World Heritage polygon that was excluded
# in the original analysis.
pa_parts[[2]] <- pa_parts[[2]] %>% dplyr::filter(WDPAID != 903141)


# 4) Combine and prepare geometries -------------------------------------

pa <- dplyr::bind_rows(pa_parts) %>%
  sf::st_transform(25832) %>%
  sf::st_make_valid() %>%
  dplyr::mutate(
    STATUS_YR_int = suppressWarnings(as.integer(STATUS_YR)),
    # Combine German and English designation text for classification.
    desig_text = stringr::str_to_lower(paste(
      dplyr::coalesce(as.character(DESIG), ""),
      dplyr::coalesce(as.character(DESIG_ENG), "")
    ))
  )


# 5) Identify special protected-area types -------------------------------

pa <- pa %>%
  dplyr::mutate(
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


# 6) Assign common publication categories --------------------------------
#
# Numeric category definitions:
#
#   1  = IUCN Ia/Ib
#   2  = IUCN II
#   3  = Biosphere Reserve core zone
#   4  = IUCN III / National Natural Monument
#   5  = IUCN IV
#   6  = IUCN V
#   7  = Natura 2000 Habitats Directive
#   8  = Natura 2000 Birds Directive
#   9  = Ramsar
#   10 = Biosphere Reserve buffer/transition
#   11 = Other / unassigned IUCN category
#   12 = Mixed / undefined

pa <- pa %>%
  dplyr::mutate(
    IUCN_CAT_numeric = dplyr::case_when(
      is_national_natural_monument ~ 4L,
      IUCN_CAT %in% c("Ia", "Ib", "Ia/Ib") ~ 1L,
      IUCN_CAT == "II" ~ 2L,
      is_biosphere_core ~ 3L,
      # Some landscape-protection records are coded as IUCN III
      # in the source data but belong to the IUCN V group used here.
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
      IUCN_CAT_numeric == 1L ~ "Ia/Ib",
      IUCN_CAT_numeric == 2L ~ "II",
      IUCN_CAT_numeric == 3L ~ "BR core",
      IUCN_CAT_numeric == 4L ~ "III / National Natural Monument",
      IUCN_CAT_numeric == 5L ~ "IV",
      IUCN_CAT_numeric == 6L ~ "V",
      IUCN_CAT_numeric == 7L ~ "N2000 Habitats",
      IUCN_CAT_numeric == 8L ~ "N2000 Birds",
      IUCN_CAT_numeric == 9L ~ "Ramsar",
      IUCN_CAT_numeric == 10L ~ "BR buffer/transition",
      IUCN_CAT_numeric == 11L ~ "Other",
      TRUE ~ "Mixed/Undefined"
    )
  )


# 7) Remove parent Biosphere Reserve polygons ----------------------------
#
# Biosphere Reserve zones are represented as separate polygons.
# The large parent polygons are therefore removed to avoid double
# representation of the same protected area.

pa_clean <- pa %>%
  dplyr::filter(!DESIG %in% c("Biosphärenreservat", "UNESCO-MAB Biosphere Reserve"))


# 8) Calculate protected-area size ---------------------------------------

pa_clean <- pa_clean %>%
  dplyr::mutate(area_km2 = as.numeric(sf::st_area(geometry)) / 1e6)


# 9) Keep variables required by downstream analyses ---------------------

pa_clean <- pa_clean %>%
  dplyr::select(
    WDPAID, NAME, DESIG, DESIG_ENG, IUCN_CAT, STATUS_YR, STATUS_YR_int,
    IUCN_CAT_numeric, IUCN_label, area_km2, geometry
  )


# 10) Save common protected-area layer -----------------------------------

sf::st_write(pa_clean, PATH_PA_GPKG, layer = LAYER_PA, delete_dsn = TRUE, quiet = TRUE)


# 11) Simple QC -----------------------------------------------------------

pa_qc <- pa_clean %>%
  sf::st_drop_geometry() %>%
  dplyr::count(IUCN_CAT_numeric, IUCN_label, name = "n_polygons") %>%
  dplyr::arrange(IUCN_CAT_numeric)

# pa_qc is a plain data.frame here (st_drop_geometry() drops the tibble
# class too), so print()'s tibble-only `n` argument is not available -
# passing it partial-matches to print.default's `na.print` and errors
# ("invalid specification for 'na.print'"). Not needed anyway: pa_qc has
# at most a dozen rows (one per IUCN_CAT_numeric), so it always prints
# in full.
print(pa_qc)

readr::write_csv(pa_qc, file.path(OUT_PA_DIR, "QC_protected_area_categories.csv"))

cat("\nProtected-area preparation completed.\n")
cat("Output:", PATH_PA_GPKG, "\n")
