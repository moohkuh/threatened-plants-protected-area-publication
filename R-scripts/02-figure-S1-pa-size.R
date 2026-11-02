# 02-figure-S1-pa-size.R
# Figure S1: size distribution of protected areas in Germany.
#   a) map of all PA polygons, coloured by size (< 25 km² vs. >= 25 km²)
#   b) share of the summed PA area in each size class
#   c) share of PAs per protected-area category in each size class,
#      with the median PA size per category
# Needs the PA layer written by 01-prepare-protected-areas.R.


# 0) Setup ---------------------------------------------------------------

suppressPackageStartupMessages({
  library(sf)
  library(dplyr)
  library(ggplot2)
  library(patchwork)
  library(scales)
  library(sysfonts)
  library(showtext)
})

sf::sf_use_s2(FALSE)

font_add_google("Roboto Condensed", "roboto_condensed")
showtext_auto()
base_family <- "roboto_condensed"

PATH_PA_GPKG <- "./data/derived/protected_areas/PA_polygons_slim.gpkg"
PATH_DE      <- "./data/spatial/DE-boundaries/VG2500_STA.shp"
OUT_FIG_DIR  <- "./figures/protected_areas"

dir.create(OUT_FIG_DIR, recursive = TRUE, showWarnings = FALSE)

size_cols <- c("< 25 km²" = "#74c476", "≥ 25 km²" = "#006837")


# 1) Data ----------------------------------------------------------------

pa_plot <- sf::st_read(PATH_PA_GPKG, layer = "pa_polygons_slim", quiet = TRUE) %>%
  sf::st_transform(25832)

pa_plot$area_km2 <- as.numeric(sf::st_area(pa_plot)) / 1e6
pa_plot$size_class <- dplyr::if_else(pa_plot$area_km2 < 25, "< 25 km²", "≥ 25 km²")

# most recent German boundary polygon
de_poly <- sf::st_read(PATH_DE, quiet = TRUE) %>%
  dplyr::mutate(BEGINN = as.Date(BEGINN)) %>%
  dplyr::arrange(dplyr::desc(BEGINN)) %>%
  dplyr::slice(1) %>%
  sf::st_transform(25832) %>%
  sf::st_make_valid()

# categories for panel c: the focal categories of the analysis, taken from
# IUCN_label assigned in 01-prepare-protected-areas.R
pa_all <- pa_plot %>%
  sf::st_drop_geometry() %>%
  dplyr::mutate(
    category = dplyr::case_when(
      IUCN_label %in% c("Ia/Ib", "II", "BR core",
                        "III / National Natural Monument") ~ "Strict",
      IUCN_label == "IV" ~ "IUCN IV",
      IUCN_label == "V" ~ "IUCN V",
      IUCN_label == "N2000 Habitats" ~ "Natura 2000",
      IUCN_label == "N2000 Birds" ~ "N2000 Birds",
      TRUE ~ "Other"
    ),
    category = factor(category, levels = c("Strict", "IUCN IV", "IUCN V",
                                           "Natura 2000", "N2000 Birds", "Other"))
  )

# median sizes cited in the text (Discussion, Notes S1)
size_check <- pa_all %>%
  dplyr::group_by(IUCN_label) %>%
  dplyr::summarise(
    n_polygons = dplyr::n(),
    median_area_km2 = median(area_km2, na.rm = TRUE),
    .groups = "drop"
  )

print(size_check)


# 2) Panel a: map --------------------------------------------------------

p_map <- ggplot() +
  geom_sf(data = de_poly, fill = "grey90", color = "grey60", linewidth = 0.4) +
  geom_sf(data = pa_plot, aes(fill = size_class), color = NA, alpha = 0.7) +
  scale_fill_manual(values = size_cols, name = "PA size") +
  coord_sf(datum = NA, expand = FALSE) +
  theme_void(base_family = base_family) +
  theme(
    legend.position = "right",
    legend.title = element_text(size = 60, face = "bold"),
    legend.text = element_text(size = 55)
  )


# 3) Panel b: share of total PA area -------------------------------------

# polygon areas are summed, so areas covered by several overlapping
# PAs are counted more than once
pie_data <- pa_all %>%
  dplyr::group_by(size_class) %>%
  dplyr::summarise(total_area_km2 = sum(area_km2, na.rm = TRUE), .groups = "drop") %>%
  dplyr::mutate(pct = 100 * total_area_km2 / sum(total_area_km2))

print(pie_data)

p_pie <- ggplot(pie_data, aes(x = "", y = total_area_km2, fill = size_class)) +
  geom_col(width = 1, color = "white") +
  coord_polar(theta = "y") +
  geom_text(aes(label = paste0(round(pct, 1), "%")),
            position = position_stack(vjust = 0.5), size = 20) +
  scale_fill_manual(values = size_cols, name = "PA size") +
  labs(title = "Share of total PA area") +
  theme_void(base_family = base_family) +
  theme(
    legend.position = "bottom",
    legend.title = element_text(size = 45, face = "bold"),
    legend.text = element_text(size = 40),
    plot.title = element_text(hjust = 0.5, size = 50)
  )


# 4) Panel c: share of PA count per category -----------------------------

bar_data <- pa_all %>%
  dplyr::count(category, size_class) %>%
  dplyr::group_by(category) %>%
  dplyr::mutate(pct = 100 * n / sum(n)) %>%
  dplyr::ungroup()

median_area <- pa_all %>%
  dplyr::group_by(category) %>%
  dplyr::summarise(median_area = median(area_km2, na.rm = TRUE), .groups = "drop") %>%
  dplyr::mutate(label = paste0(round(median_area, 1), " km²"))

print(median_area)

p_bar <- ggplot(bar_data, aes(x = category, y = pct, fill = size_class)) +
  geom_col(color = "white", width = 0.7) +
  geom_hline(yintercept = 50, linetype = "dashed", color = "grey30", linewidth = 0.5) +
  geom_text(
    data = median_area,
    aes(x = category, y = 110, label = label),
    inherit.aes = FALSE,
    angle = 45, size = 13, hjust = 0.3, color = "grey20"
  ) +
  scale_fill_manual(values = size_cols, name = "PA size") +
  scale_y_continuous(
    labels = scales::percent_format(scale = 1),
    limits = c(0, 125),
    breaks = seq(0, 100, by = 25),
    expand = expansion(mult = c(0, 0))
  ) +
  labs(x = NULL, y = "Proportion of PAs [%]", title = "Share of PA count per category") +
  coord_cartesian(clip = "off") +
  theme_classic(base_family = base_family) +
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1, size = 40),
    axis.text.y = element_text(size = 40),
    axis.title.y = element_text(size = 45),
    plot.title = element_text(size = 45, hjust = 0.5),
    legend.position = "none"
  )


# 5) Combine and save ----------------------------------------------------

fig_s1 <- p_map / (p_pie | p_bar) +
  plot_layout(heights = c(1.5, 0.8)) +
  plot_annotation(tag_levels = "a") &
  theme(plot.tag = element_text(face = "bold", size = 60))

ggsave(
  file.path(OUT_FIG_DIR, "FigureS1_PA_size.png"),
  plot = fig_s1,
  width = 9, height = 12, dpi = 300
)
