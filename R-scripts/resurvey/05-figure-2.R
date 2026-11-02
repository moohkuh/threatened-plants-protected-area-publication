# 05-figure-2.R
# Figure 2: Violin plots + forest plot + maps with violin insets.
#
# Depends on:
#   00-preamble-resurvey.R
#   03-analyse-cover-change.R (RDS cache)
#
# Plot locations use the public approximate ReSurvey coordinates.
# Protected-area classifications use the frozen plot-level PA flags
# derived previously from the exact coordinates used in the analysis.
#
# Author: Marie Englert
# Date: 2026-04-20

source("./R-scripts/resurvey/00-preamble-resurvey.R")

suppressPackageStartupMessages({
  library(patchwork)
  library(geodata)
  library(openxlsx)
})


# ------ SETTINGS --------------------------------------------


OUT_DIR_FIG <- OUT_FIG_RS_DIR
RDS_FILE <- file.path(OUT_DIR_RDS, "species_cover_means_objects.rds")
CACHE_FILE <- file.path(OUT_DIR_RDS, "resurvey_map_cache.gpkg")

dir.create(OUT_DIR_FIG, showWarnings = FALSE, recursive = TRUE)

assert_files_exist(c(RDS_FILE, PATH_PA_GPKG, OUT_RS_CLEAN))




# STEP 1: Load results from 04 ------------------------------


cache_04 <- readRDS(RDS_FILE)

sp_all    <- cache_04$sp_all
sp_strict <- cache_04$sp_strict
sp_iucn4  <- cache_04$sp_iucn4
sp_iucn5  <- cache_04$sp_iucn5
sp_n2k    <- cache_04$sp_n2k

effects_tbl <- cache_04$effects_tbl

mm_all    <- cache_04$emm_tbl %>% dplyr::filter(scenario == "All PAs")
mm_strict <- cache_04$emm_tbl %>% dplyr::filter(scenario == "Strict")
mm_iucn4  <- cache_04$emm_tbl %>% dplyr::filter(scenario == "IUCN IV")
mm_iucn5  <- cache_04$emm_tbl %>% dplyr::filter(scenario == "IUCN V")
mm_n2k    <- cache_04$emm_tbl %>% dplyr::filter(scenario == "Natura 2000")

sp_strict_keep <- unique(cache_04$dat_strict_long$wcvp_name)
sp_n2k_keep    <- unique(cache_04$dat_n2k_long$wcvp_name)
sp_iucn4_keep  <- unique(cache_04$dat_iucn4_long$wcvp_name)
sp_iucn5_keep  <- unique(cache_04$dat_iucn5_long$wcvp_name)


# Colors
cols_forest <- c(
  "All PAs"     = "#5D7388FF",
  "Strict"      = "#8AAE8AFF",
  "IUCN IV"     = "#4F4D1DFF",
  "IUCN V"      = "#A08F27FF",
  "Natura 2000" = "#E5A729FF"
)
col_outside       <- "grey75"
col_all_inside    <- unname(cols_forest["All PAs"])
col_strict_inside <- unname(cols_forest["Strict"])
col_iucn4_inside  <- unname(cols_forest["IUCN IV"])
col_iucn5_inside  <- unname(cols_forest["IUCN V"])
col_n2k_inside    <- unname(cols_forest["Natura 2000"])



# STEP 2: Violin panels --------------------------------------


#  ----------- (a) All PAs ------------
y_all      <- range(c(sp_all$mean_rel, mm_all$lower, mm_all$upper), na.rm = TRUE)


p_box_all <- ggplot(sp_all, aes(x = side, y = mean_rel, fill = side, color = side)) +
  geom_violin(trim = TRUE, alpha = 0.8) +
  geom_errorbar(
    data = mm_all,
    aes(x = side, ymin = lower, ymax = upper),
    inherit.aes = FALSE,
    width = 0.1, linewidth = 0.7, color = "black"
  ) +
  geom_point(
    data = mm_all,
    aes(x = side, y = mean),
    inherit.aes = FALSE,
    size = 1.7, shape = 21, fill = "black", color = "black"
  ) +
  geom_hline(yintercept = 0, linetype = "dashed", linewidth = 0.5) +
  scale_color_manual(values = c(inside = col_all_inside, outside = col_outside)) +
  scale_fill_manual(values  = c(inside = col_all_inside, outside = col_outside)) +
  theme_classic(base_family = base_family) +
  labs(x = "Plot location", y = "Mean relative cover change per year") +
  theme(
    plot.margin     = margin(t = 18, r = 8, b = 8, l = 25),
    plot.title      = element_text(size = 50, family = base_family),
    legend.position = "none",
    axis.text.x     = element_text(size = 35, family = base_family),
    axis.title.x    = element_text(size = 45, family = base_family),
    axis.text.y     = element_text(size = 35, family = base_family),
    axis.title.y    = element_text(size = 45, family = base_family)
  )

# --------- (b) Strict -------
y_strict      <- range(c(sp_strict$mean_rel, mm_strict$lower, mm_strict$upper), na.rm = TRUE)

p_box_strict <- ggplot(sp_strict, aes(x = side, y = mean_rel, fill = side, color = side)) +
  geom_violin(trim = TRUE, alpha = 0.8) +
  geom_errorbar(
    data = mm_strict,
    aes(x = side, ymin = lower, ymax = upper),
    inherit.aes = FALSE,
    width = 0.1, linewidth = 0.7, color = "black"
  ) +
  geom_point(
    data = mm_strict,
    aes(x = side, y = mean),
    inherit.aes = FALSE,
    size = 1.7, shape = 21, fill = "black", color = "black"
  ) +
  geom_hline(yintercept = 0, linetype = "dashed", linewidth = 0.5) +
  scale_color_manual(values = c(inside = col_strict_inside, outside = col_outside)) +
  scale_fill_manual(values  = c(inside = col_strict_inside, outside = col_outside)) +
  theme_classic(base_family = base_family) +
  ggtitle("Strict protection (IUCN I-III + BR cores)") +
  labs(x = "Plot location", y = "Mean relative cover change per year") +
  theme(
    plot.title      = element_text(size = 35, family = base_family),
    legend.position = "none",
    axis.text.x     = element_text(size = 30, family = base_family),
    axis.title.x    = element_text(size = 35, family = base_family),
    axis.text.y     = element_text(size = 30, family = base_family),
    axis.title.y    = element_text(size = 35, family = base_family)
  )

# ------- (c) IUCN IV -------
y_i4      <- range(c(sp_iucn4$mean_rel, mm_iucn4$lower, mm_iucn4$upper), na.rm = TRUE)


p_box_iucn4 <- ggplot(sp_iucn4, aes(x = side, y = mean_rel, fill = side, color = side)) +
  geom_violin(trim = TRUE, alpha = 0.8) +
  geom_errorbar(
    data = mm_iucn4,
    aes(x = side, ymin = lower, ymax = upper),
    inherit.aes = FALSE,
    width = 0.1, linewidth = 0.7, color = "black"
  ) +
  geom_point(
    data = mm_iucn4,
    aes(x = side, y = mean),
    inherit.aes = FALSE,
    size = 1.7, shape = 21, fill = "black", color = "black"
  ) +
  geom_hline(yintercept = 0, linetype = "dashed", linewidth = 0.5) +
  scale_color_manual(values = c(inside = col_iucn4_inside, outside = col_outside)) +
  scale_fill_manual(values  = c(inside = col_iucn4_inside, outside = col_outside)) +
  theme_classic(base_family = base_family) +
  ggtitle("IUCN IV protection") +
  labs(x = "Plot location", y = "Mean relative cover change per year") +
  theme(
    plot.title      = element_text(size = 35, family = base_family),
    legend.position = "none",
    axis.text.x     = element_text(size = 30, family = base_family),
    axis.title.x    = element_text(size = 35, family = base_family),
    axis.text.y     = element_text(size = 30, family = base_family),
    axis.title.y    = element_text(size = 35, family = base_family)
  )

# ---------- (d) IUCN V --------------
y_i5      <- range(c(sp_iucn5$mean_rel, mm_iucn5$lower, mm_iucn5$upper), na.rm = TRUE)


p_box_iucn5 <- ggplot(sp_iucn5, aes(x = side, y = mean_rel, fill = side, color = side)) +
  geom_violin(trim = TRUE, alpha = 0.8) +
  geom_errorbar(
    data = mm_iucn5,
    aes(x = side, ymin = lower, ymax = upper),
    inherit.aes = FALSE,
    width = 0.1, linewidth = 0.7, color = "black"
  ) +
  geom_point(
    data = mm_iucn5,
    aes(x = side, y = mean),
    inherit.aes = FALSE,
    size = 1.7, shape = 21, fill = "black", color = "black"
  ) +
  geom_hline(yintercept = 0, linetype = "dashed", linewidth = 0.5) +
  scale_color_manual(values = c(inside = col_iucn5_inside, outside = col_outside)) +
  scale_fill_manual(values  = c(inside = col_iucn5_inside, outside = col_outside)) +
  theme_classic(base_family = base_family) +
  ggtitle("IUCN V protection") +
  labs(x = "Plot location", y = "Mean relative cover change per year") +
  theme(
    plot.title      = element_text(size = 35, family = base_family),
    legend.position = "none",
    axis.text.x     = element_text(size = 30, family = base_family),
    axis.title.x    = element_text(size = 35, family = base_family),
    axis.text.y     = element_text(size = 30, family = base_family),
    axis.title.y    = element_text(size = 35, family = base_family)
  )

# --------- (e) Natura 2000 ------------
y_n2k      <- range(c(sp_n2k$mean_rel, mm_n2k$lower, mm_n2k$upper), na.rm = TRUE)


p_box_n2k <- ggplot(sp_n2k, aes(x = side, y = mean_rel, fill = side, color = side)) +
  geom_violin(trim = TRUE, alpha = 0.8) +
  geom_errorbar(
    data = mm_n2k,
    aes(x = side, ymin = lower, ymax = upper),
    inherit.aes = FALSE,
    width = 0.1, linewidth = 0.7, color = "black"
  ) +
  geom_point(
    data = mm_n2k,
    aes(x = side, y = mean),
    inherit.aes = FALSE,
    size = 1.7, shape = 21, fill = "black", color = "black"
  ) +
  geom_hline(yintercept = 0, linetype = "dashed", linewidth = 0.5) +
  scale_color_manual(values = c(inside = col_n2k_inside, outside = col_outside)) +
  scale_fill_manual(values  = c(inside = col_n2k_inside, outside = col_outside)) +
  theme_classic(base_family = base_family) +
  ggtitle("N2k protection") +
  labs(x = "Plot location", y = "Mean relative cover change per year") +
  theme(
    plot.title      = element_text(size = 35, family = base_family),
    legend.position = "none",
    axis.text.x     = element_text(size = 30, family = base_family),
    axis.title.x    = element_text(size = 35, family = base_family),
    axis.text.y     = element_text(size = 30, family = base_family),
    axis.title.y    = element_text(size = 35, family = base_family)
  )



# STEP 3: Forest plot ---------------------------------------------

# Number of species included in each paired cover-change comparison
n_species <- c(
  "All PAs"     = dplyr::n_distinct(sp_all$wcvp_name),
  "Strict"      = dplyr::n_distinct(sp_strict$wcvp_name),
  "IUCN IV"     = dplyr::n_distinct(sp_iucn4$wcvp_name),
  "IUCN V"      = dplyr::n_distinct(sp_iucn5$wcvp_name),
  "Natura 2000" = dplyr::n_distinct(sp_n2k$wcvp_name)
)

# Labels shown on the forest-plot y-axis
scenario_labels <- c(
  "All PAs"     = paste0("All (n = ", n_species["All PAs"], ")"),
  "Strict"      = paste0("Strict (n = ", n_species["Strict"], ")"),
  "IUCN IV"     = paste0("IUCN IV (n = ", n_species["IUCN IV"], ")"),
  "IUCN V"      = paste0("IUCN V (n = ", n_species["IUCN V"], ")"),
  "Natura 2000" = paste0("Natura 2000 (n = ", n_species["Natura 2000"], ")")
)


# Prepare contrast table for plotting
effects_plot <- effects_tbl %>%
  dplyr::mutate(
    scenario = factor(scenario, levels = c("Natura 2000", "IUCN V", "IUCN IV", "Strict", "All PAs")),
    sig_flag = dplyr::case_when(p_value < 0.05 ~ "significant", TRUE ~ "not significant"),
    sig_flag = factor(sig_flag, levels = c("significant", "not significant"))
  )

p_forest <- ggplot(effects_plot) +
  geom_vline(xintercept = 0, color = "black", linetype = "dashed", linewidth = 0.6) +

  # Draw confidence intervals first
  geom_errorbar(
    aes(y = scenario, x = estimate, xmin = lower, xmax = upper, color = scenario, linetype = sig_flag),
    orientation = "y", height = 0.15, linewidth = 0.8
  ) +

  # Draw hollow points on top of confidence intervals
  geom_point(
    aes(y = scenario, x = estimate, color = scenario),
    shape = 21, fill = "white", size = 3, stroke = 1.1
  ) +

  scale_color_manual(values = cols_forest) +
  scale_y_discrete(labels = scenario_labels) +
  scale_linetype_manual(values = c("significant" = "solid", "not significant" = "dashed"), guide = "none") +
  labs(x = "Effect of protection", y = NULL) +
  theme_classic(base_family = base_family) +
  theme(
    plot.margin = margin(t = 18, r = 25, b = 8, l = 8),
    legend.position = "none",
    plot.title = element_text(size = 50, family = base_family),
    axis.text.x = element_text(size = 35, family = base_family),
    axis.title.x = element_text(size = 45, family = base_family),
    axis.text.y = element_text(size = 35, family = base_family),
    axis.title.y = element_text(size = 45, family = base_family)
  )



# STEP 4: Map cache (run once) --------------------------------------


build_map_cache <- function(force = FALSE) {
  
  if (file.exists(CACHE_FILE) && !force) return(invisible(CACHE_FILE))
  
  pa <- sf::st_read(PATH_PA_GPKG, layer = LAYER_PA, quiet = TRUE)
  pa_use <- pa %>%
    dplyr::select(WDPAID, NAME, DESIG, IUCN_CAT, IUCN_label, STATUS_YR_int, geom)
  if ("geom" %in% names(pa_use)) sf::st_geometry(pa_use) <- "geom"
  rm(pa)

  de_states_25832 <- geodata::gadm(country = "DEU", level = 1, path = "./data/spatial/DE-boundaries/GADM") %>%
    sf::st_as_sf() %>%
    sf::st_transform(25832) %>%
    sf::st_make_valid() %>%
    dplyr::transmute(state = NAME_1, geometry) %>%
    dplyr::distinct(state, .keep_all = TRUE)

  rs_clean <- readr::read_csv(OUT_RS_CLEAN, show_col_types = FALSE)

  # Build plot points using public coordinates and frozen PA flags
  pts_class <- rs_clean %>%
    dplyr::filter(is_nat_rl == TRUE, is.finite(LONGITUDE), is.finite(LATITUDE)) %>%
    dplyr::group_by(RS_PROJECT_PLOT, wcvp_name, LONGITUDE, LATITUDE) %>%
    dplyr::summarise(
      any_PA = dplyr::first(any_PA),
      any_strict = dplyr::first(any_strict),
      any_n2k = dplyr::first(any_n2k),
      any_iucn4 = dplyr::first(any_iucn4),
      any_iucn5 = dplyr::first(any_iucn5),
      .groups = "drop"
    ) %>%
    dplyr::mutate(
      side_strict = dplyr::case_when(any_strict ~ "inside", !any_PA ~ "outside", TRUE ~ NA_character_),
      side_n2k = dplyr::case_when(any_n2k ~ "inside", !any_PA ~ "outside", TRUE ~ NA_character_),
      side_iucn4 = dplyr::case_when(any_iucn4 ~ "inside", !any_PA ~ "outside", TRUE ~ NA_character_),
      side_iucn5 = dplyr::case_when(any_iucn5 ~ "inside", !any_PA ~ "outside", TRUE ~ NA_character_)
    ) %>%
    sf::st_as_sf(coords = c("LONGITUDE", "LATITUDE"), crs = 4326, remove = FALSE) %>%
    sf::st_transform(25832)
  
  pa_groups_simple <- pa_all_25832 %>%
    sf::st_make_valid() %>%
    dplyr::group_by(pa_group) %>%
    dplyr::summarise(.groups = "drop") %>%
    sf::st_as_sf() %>%
    sf::st_simplify(dTolerance = 1000, preserveTopology = TRUE)
  
  pa_strict_simple <- pa_groups_simple %>% dplyr::filter(pa_group == "strict")
  pa_n2k_simple    <- pa_groups_simple %>% dplyr::filter(pa_group == "N2000")
  pa_iucn4_simple  <- pa_groups_simple %>% dplyr::filter(pa_group == "IUCN IV")
  pa_iucn5_simple  <- pa_groups_simple %>% dplyr::filter(pa_group == "IUCN V")
  
  pts_strict_map <- pts_class %>% dplyr::filter(wcvp_name %in% sp_strict_keep, !is.na(side_strict))
  pts_n2k_map    <- pts_class %>% dplyr::filter(wcvp_name %in% sp_n2k_keep,    !is.na(side_n2k))
  pts_iucn4_map  <- pts_class %>% dplyr::filter(wcvp_name %in% sp_iucn4_keep,  !is.na(side_iucn4))
  pts_iucn5_map  <- pts_class %>% dplyr::filter(wcvp_name %in% sp_iucn5_keep,  !is.na(side_iucn5))
  
  if (file.exists(CACHE_FILE)) file.remove(CACHE_FILE)
  
  sf::st_write(de_states_25832,  CACHE_FILE, layer = "de_states_25832",  quiet = TRUE)
  sf::st_write(pa_strict_simple, CACHE_FILE, layer = "pa_strict_simple", quiet = TRUE)
  sf::st_write(pa_n2k_simple,    CACHE_FILE, layer = "pa_n2k_simple",    quiet = TRUE)
  sf::st_write(pa_iucn4_simple,  CACHE_FILE, layer = "pa_iucn4_simple",  quiet = TRUE)
  sf::st_write(pa_iucn5_simple,  CACHE_FILE, layer = "pa_iucn5_simple",  quiet = TRUE)
  sf::st_write(pts_strict_map,   CACHE_FILE, layer = "pts_strict_map",   quiet = TRUE)
  sf::st_write(pts_n2k_map,      CACHE_FILE, layer = "pts_n2k_map",      quiet = TRUE)
  sf::st_write(pts_iucn4_map,    CACHE_FILE, layer = "pts_iucn4_map",    quiet = TRUE)
  sf::st_write(pts_iucn5_map,    CACHE_FILE, layer = "pts_iucn5_map",    quiet = TRUE)
  
  invisible(CACHE_FILE)
}

if (!file.exists(CACHE_FILE)) {
  build_map_cache()
}



# STEP 5: Load map cache ------------------------------------------


stopifnot(file.exists(CACHE_FILE))
de_border <- geodata::gadm(
  country = "DEU",
  level = 0,
  path = "./data/spatial/DE-boundaries/GADM"
) %>%
  sf::st_as_sf() %>%
  sf::st_transform(25832) %>%
  sf::st_make_valid() %>%
  dplyr::transmute(geometry) %>%
  # dplyr::distinct(state, .keep_all = TRUE) %>%
  sf::st_simplify(dTolerance = 2000)
#de_states_25832  <- sf::st_read(CACHE_FILE, layer = "de_states_25832",  quiet = TRUE)
pa_strict_simple <- sf::st_read(CACHE_FILE, layer = "pa_strict_simple", quiet = TRUE)
pa_n2k_simple    <- sf::st_read(CACHE_FILE, layer = "pa_n2k_simple",    quiet = TRUE)
pa_iucn4_simple  <- sf::st_read(CACHE_FILE, layer = "pa_iucn4_simple",  quiet = TRUE)
pa_iucn5_simple  <- sf::st_read(CACHE_FILE, layer = "pa_iucn5_simple",  quiet = TRUE)
pts_strict_map   <- sf::st_read(CACHE_FILE, layer = "pts_strict_map",   quiet = TRUE)
pts_n2k_map      <- sf::st_read(CACHE_FILE, layer = "pts_n2k_map",      quiet = TRUE)
pts_iucn4_map    <- sf::st_read(CACHE_FILE, layer = "pts_iucn4_map",    quiet = TRUE)
pts_iucn5_map    <- sf::st_read(CACHE_FILE, layer = "pts_iucn5_map",    quiet = TRUE)

bb               <- sf::st_bbox(de_border)
xrange           <- bb$xmax - bb$xmin
bb_inset         <- bb
bb_inset["xmax"] <- bb["xmax"] + 0.6 * xrange

coord_de_inset <- ggplot2::coord_sf(
  xlim   = c(bb_inset["xmin"], bb_inset["xmax"]),
  ylim   = c(bb_inset["ymin"], bb_inset["ymax"]),
  datum  = NA, expand = FALSE, clip = "off"
)
# "#6E6E6E"
#col_outside <- "#C3BF9DFF"
col_strict  <- "#8AAE8AFF"
col_n2k     <- "#E5A729FF"
col_iucn4   <- "#4F4D1DFF"
col_iucn5   <- "#A08F27FF"



# STEP 6: Maps ---------------------------------------------------------------


p_map_strict <- ggplot() +
  geom_sf(data = de_border,  fill = "grey95", color = "grey70", linewidth = 0.4) +
  geom_sf(data = pa_strict_simple, fill = col_strict, color = NA, alpha = 0.7) +
  geom_sf(data = pts_strict_map, aes(fill = side_strict),
          shape = 21, size = 3, color = "grey30", stroke = 0.3, alpha = 0.8) +
  scale_fill_manual(values = c(inside = col_strict, outside = col_outside),
                    name = "Plot location") +
  coord_sf(datum = NA, expand = FALSE) +
  labs(x = NULL, y = NULL) +
  theme_void(base_family = base_family) +
  theme(plot.margin = unit(c(0,0,0,0), "cm"), legend.position = "right",
        panel.grid.major = element_blank(),
        legend.title = element_text(size = 45, family = base_family),
        legend.text  = element_text(size = 40, family = base_family))

p_map_n2k <- ggplot() +
  geom_sf(data = de_border, fill = "grey95", color = "grey70", linewidth = 0.4) +
  geom_sf(data = pa_n2k_simple,   fill = col_n2k, color = NA, alpha = 0.7) +
  geom_sf(data = pts_n2k_map, aes(fill = side_n2k),
          shape = 21, size = 3, color = "grey30", stroke = 0.3, alpha = 0.8) +
  scale_fill_manual(values = c(inside = col_n2k, outside = col_outside),
                    name = "Plot location") +
  coord_sf(datum = NA, expand = FALSE) +
  labs(x = NULL, y = NULL) +
  theme_void(base_family = base_family) +
  theme(plot.margin = unit(c(0,0,0,0), "cm"), legend.position = "right",
        panel.grid.major = element_blank(),
        legend.title = element_text(size = 45, family = base_family),
        legend.text  = element_text(size = 40, family = base_family))

p_map_iucn4 <- ggplot() +
  geom_sf(data = de_border,  fill = "grey95", color = "grey70", linewidth = 0.4) +
  geom_sf(data = pa_iucn4_simple,  fill = col_iucn4, color = NA, alpha = 0.7) +
  geom_sf(data = pts_iucn4_map, aes(fill = side_iucn4),
          shape = 21, size = 3, color = "grey30", stroke = 0.3, alpha = 0.6) +
  scale_fill_manual(values = c(inside = col_iucn4, outside = col_outside),
                    name = "Plot location") +
  coord_sf(datum = NA, expand = FALSE) +
  labs(x = NULL, y = NULL) +
  theme_void(base_family = base_family) +
  theme(plot.margin = unit(c(0,0,0,0), "cm"), legend.position = "right",
        panel.grid.major = element_blank(),
        legend.title = element_text(size = 45, family = base_family),
        legend.text  = element_text(size = 40, family = base_family))

p_map_iucn5 <- ggplot() +
  geom_sf(data = de_border,  fill = "grey95", color = "grey70", linewidth = 0.4) +
  geom_sf(data = pa_iucn5_simple,  fill = col_iucn5, color = NA, alpha = 0.7) +
  geom_sf(data = pts_iucn5_map, aes(fill = side_iucn5),
          shape = 21, size = 3, color = "grey30", stroke = 0.3, alpha = 0.8) +
  scale_fill_manual(values = c(inside = col_iucn5, outside = col_outside),
                    name = "Plot location") +
  coord_sf(datum = NA, expand = FALSE) +
  labs(x = NULL, y = NULL) +
  theme_void(base_family = base_family) +
  theme(plot.margin = unit(c(0,0,0,0), "cm"), legend.position = "right",
        panel.grid.major = element_blank(),
        legend.title = element_text(size = 45, family = base_family),
        legend.text  = element_text(size = 40, family = base_family))



# STEP 7: Helper functions for insets -----------------------------------


add_inset_space <- function(p) {
  p + coord_de_inset
}

make_inset <- function(p) {
  p +
    labs(tag = NULL) +
    theme(
      legend.position  = "none",
      axis.title       = element_blank(),
      axis.text        = element_text(size = 7, family = base_family),
      axis.ticks       = element_line(linewidth = 0.2),
      panel.grid       = element_blank(),
      plot.title       = element_blank(),
      plot.subtitle    = element_blank(),
      plot.caption     = element_blank(),
      plot.margin      = margin(0, 0, 0, 0),
      panel.background = element_blank(),
      plot.background  = element_blank(),
      plot.tag         = element_blank()
    )
}

add_inset <- function(p_map, p_inset,
                      left = 0.65, bottom = 0.02,
                      right = 0.985, top = 0.98) {
  p_map +
    inset_element(make_inset(p_inset),
                  left = left, bottom = bottom,
                  right = right, top = top,
                  align_to = "panel")
}

theme_map_legend_top <- function() {
  theme(legend.position = "none")
}



# STEP 8: Assemble maps with violin insets -------------------------------


p_box_strict <- p_box_strict + labs(tag = NULL, title = NULL)
p_box_iucn4  <- p_box_iucn4  + labs(tag = NULL, title = NULL)
p_box_iucn5  <- p_box_iucn5  + labs(tag = NULL, title = NULL)
p_box_n2k    <- p_box_n2k    + labs(tag = NULL, title = NULL)

p_map_strict_space <- add_inset_space(p_map_strict) + theme_map_legend_top()
p_map_iucn4_space  <- add_inset_space(p_map_iucn4)  + theme_map_legend_top()
p_map_iucn5_space  <- add_inset_space(p_map_iucn5)  + theme_map_legend_top()
p_map_n2k_space    <- add_inset_space(p_map_n2k)    + theme_map_legend_top()

m_strict <- add_inset(p_map_strict_space, p_box_strict)
m_i4     <- add_inset(p_map_iucn4_space,  p_box_iucn4)
m_i5     <- add_inset(p_map_iucn5_space,  p_box_iucn5)
m_n2k    <- add_inset(p_map_n2k_space,    p_box_n2k)


# ------------------------------------------------------------
# STEP 9: Tag themes + final assembly
# ------------------------------------------------------------

tag_theme_big <- theme(
  plot.margin       = margin(0, 0, 0, 0),
  plot.tag          = element_text(size = rel(4), family = base_family,
                                   face = "plain"),
  plot.tag.position = "topleft",
  plot.tag.location = "plot"
)

tag_theme_top <- theme(
  plot.margin       = margin(0, 0, 0, 0),
  plot.tag.padding  = unit(0.2, "lines"),
  plot.tag          = element_text(size = rel(4), face = "plain",
                                   family = base_family),
  plot.tag.position = "topleft",
  plot.tag.location = "plot"
)

map_header_theme <- theme(
  aspect.ratio = NULL,
  plot.title   = element_text(size = 50, hjust = 0, family = base_family),
  plot.margin  = margin(2, 0, 0, 0)
)

reserve_title_space <- theme(
  plot.title = element_text(size = 50, colour = "transparent",
                            family = base_family, margin = margin(b = 8))
)

# m_strict <- patchwork::wrap_elements(m_strict) +
#   labs(tag = "c) Strict protection") + tag_theme_big
#library(ggtext)

# m_strict <- patchwork::wrap_elements(m_strict) +
# labs(
#    tag = paste0(
#      "<span style='font-size:44pt;'>c) Strict protection</span>",
#       "<br><br>",
#      "<span style='font-size:30pt;'>",
#       "\u2003\u2003",  # Two em spaces for indentation
#       "IUCN I–III + BR cores",
#       "</span>"
#     )
#   ) +
# theme(
#     # Render HTML/markdown formatting inside the panel tag
#     plot.tag = ggtext::element_markdown(
#        family = base_family,
#        hjust = 0,
#        vjust = 1,
#        lineheight = 0.70
#      ),
# 
#      # Position the complete two-line tag inside the existing panel
#      plot.tag.position = "topleft",
#      #plot.tag.location = "topleft",
#      plot.margin = margin(0, 0, 0, 0)
#    )
# Main panel title
m_strict <- patchwork::wrap_elements(m_strict) +
  labs(tag = "c) Strict protection") +
  tag_theme_big
# Create the subtitle as an independent transparent plot
strict_subtitle <- ggplot() +
  ggtext::geom_richtext(
    aes(
      x = 0,
      y = 1,
      label = "IUCN I–III + BR cores"
    ),
    hjust = 0,
    vjust = 1,
    size = 33 / ggplot2::.pt,  # Approximately 30 pt
    family = base_family,
    fill = NA,
    label.color = NA,
    label.padding = grid::unit(rep(0, 4), "pt")
  ) +
  xlim(0, 1) +
  ylim(0, 1) +
  theme_void() +
  theme(
    plot.background = element_rect(
      fill = NA,
      color = NA
    ),
    panel.background = element_rect(
      fill = NA,
      color = NA
    ),
    plot.margin = margin(0, 0, 0, 0)
  )

# Overlay the subtitle without changing the map layout
m_strict <- m_strict +
  patchwork::inset_element(
    strict_subtitle,
    left = 0.097,
    bottom = 0.9,
    right = 0.55,
    top = 0.95,
    align_to = "full",
    on_top = TRUE
  )

m_i4 <- patchwork::wrap_elements(m_i4) +
  labs(tag = "d) IUCN IV") + tag_theme_big

m_i5 <- patchwork::wrap_elements(m_i5) +
  labs(tag = "e) IUCN V") + map_header_theme + tag_theme_big

m_n2k <- patchwork::wrap_elements(m_n2k) +
  labs(tag = "f) Natura 2000") + tag_theme_big

middle_row <- (m_strict | m_i4) +
  plot_layout(widths = c(1, 1), guides = "keep") &
  theme(plot.margin = margin(t = 0, r = 0, b = 10, l = -20))

bottom_row <- (m_i5 | m_n2k) +
  plot_layout(widths = c(1, 1), guides = "keep") &
  theme(plot.margin = margin(t = 0, r = 0, b = 0, l = -20))

p_forest_clean <- p_forest +
  theme(axis.text.y = element_text(margin = margin(r = 0),
                                   family = base_family),
        plot.margin = margin(8, 0, 0, 0)) +
  coord_cartesian(clip = "on")

top_row <- (
  (p_box_all +
     labs(tag = "a) All protected areas") +
     tag_theme_top + ggtitle(" ") + reserve_title_space) |
    (p_forest_clean +
       labs(tag = "b) Inside–outside contrast by category") +
       tag_theme_top + ggtitle(" ") + reserve_title_space)
) +
  plot_layout(widths = c(1, 1), guides = "keep") &
  theme(plot.margin = margin(t = 8, r = 18, b = 18, l = 18))

p_final_space <- (top_row / middle_row / bottom_row) +
  plot_layout(heights = c(1.6, 2.4, 2.4), widths = c(1, 1, 1)) &
  theme(text = element_text(family = base_family))

ggsave(
  filename = file.path(OUT_DIR_FIG, "Figure2_with_maps_insets.png"),
  plot     = p_final_space,
  width    = 12,
  height   = 15,
  units    = "in",
  dpi      = 280
)

