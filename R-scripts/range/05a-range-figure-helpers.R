# 05a-range-figure-helpers.R
# Shared helpers for the range-scale publication figures.
#
# PURPOSE
#   Keep the visual design of the previous working range-figure
#   helper and only replace the data-loading part with outputs from
#   the new publication pipeline.
#
# IMPORTANT
#   Plot construction, panel proportions, margins, text sizes,
#   line styles and raincloud geometry intentionally follow the
#   previous working helper.


# 0) Packages ---------------------------------------------------------------

suppressPackageStartupMessages({
  library(dplyr)
  library(purrr)
  library(tibble)
  library(ggplot2)
  library(sf)
  library(see)
  library(patchwork)
})

if (!requireNamespace("ggrastr", quietly = TRUE)) {
  stop("Package 'ggrastr' is required for the range maps.")
}

sf::sf_use_s2(FALSE)


# 1) Paths --------------------------------------------------------------------

RANGE_OUT_DIR <- "./data/derived/range_repro"
RANGE_MODEL_DIR <- file.path(RANGE_OUT_DIR, "models")
RANGE_FIG_DIR <- "./figures/range"

PATH_MODEL_OBJECTS <- file.path(RANGE_MODEL_DIR, "range_lm_objects.rds")
PATH_SPATIAL <- file.path(RANGE_OUT_DIR, "range_spatial.gpkg")
PATH_PA <- file.path(RANGE_OUT_DIR, "grid_PA_coverage_present.rds")
PATH_URBAN <- file.path(RANGE_OUT_DIR, "grid_urban_1960.rds")
PATH_DE <- "./data/spatial/DE-boundaries/VG2500_STA.shp"

dir.create(RANGE_FIG_DIR, recursive = TRUE, showWarnings = FALSE)

required_files <- c(PATH_MODEL_OBJECTS, PATH_SPATIAL, PATH_PA, PATH_URBAN, PATH_DE)
missing_files <- required_files[!file.exists(required_files)]

if (length(missing_files) > 0) {
  stop("Required figure input(s) are missing:\n", paste(missing_files, collapse = "\n"))
}


# 2) Font ----------------------------------------------------------------------
#
# The previous working figures used the Roboto Condensed alias
# "roboto_condensed". Register the same alias here.

base_family_big <- "roboto_condensed"

if (requireNamespace("sysfonts", quietly = TRUE) && requireNamespace("showtext", quietly = TRUE)) {

  if (!base_family_big %in% sysfonts::font_families()) {
    sysfonts::font_add_google("Roboto Condensed", base_family_big)
  }

  showtext::showtext_auto()

} else {
  stop("Packages 'sysfonts' and 'showtext' are required to use Roboto Condensed.")
}


# 3) Figure defaults from the previous working helper --------------------------

na_col <- "#E5E5E5"

pt_col <- "#2F436C"
line_col <- "#00336F"
rib_col <- "#E3CE5C"

urban_col <- "grey75"

high_col <- "#2F436C"


# 4) Load model outputs ---------------------------------------------------------

load_range_models <- function(analysis_id) {

  all_objects <- readRDS(PATH_MODEL_OBJECTS)

  if (!analysis_id %in% names(all_objects)) {
    stop(
      "Analysis not found in range_lm_objects.rds: ", analysis_id,
      "\nAvailable analyses: ", paste(names(all_objects), collapse = ", ")
    )
  }

  all_objects[[analysis_id]]
}

get_pa_result <- function(pa_object) {

  out <- pa_object$results %>%
    dplyr::filter(term == "range protected")

  if (nrow(out) != 1) {
    stop("Expected exactly one 'range protected' result row.")
  }

  out
}

get_pa_p <- function(pa_object) {
  get_pa_result(pa_object)$p_value[[1]]
}

get_lm_label <- function(pa_object, digits_r2 = 2, digits_p = 3) {

  z <- get_pa_result(pa_object)

  beta <- z$beta[[1]]
  p <- z$p_value[[1]]
  r2 <- z$r_squared[[1]]

  p_txt <- if (is.na(p)) {
    "p = NA"
  } else if (p < 0.001) {
    "p < 0.001"
  } else {
    paste0("p = ", signif(p, digits_p))
  }

  paste(
    paste0("R\u00b2 = ", sprintf(paste0("%.", digits_r2, "f"), r2)),
    paste0("\u03b2 = ", round(beta, 3)),
    p_txt,
    sep = ", "
  )
}


# 5) Load map data from the new pipeline ----------------------------------------

load_range_map_data <- function() {

  coverage_cols <- c("cov_anyPA", "cov_strict", "cov_IV", "cov_V", "cov_n2000_habitats_only")

  # Load the same grid polygons used by the new range pipeline.
  grid <- sf::st_read(PATH_SPATIAL, layer = "grid_centroid_in_DE", quiet = TRUE) %>%
    sf::st_transform(25832) %>%
    dplyr::mutate(MTB_Q = as.character(MTB_Q))

  # Load protected-area coverage.
  grid_cov <- readRDS(PATH_PA) %>%
    dplyr::mutate(MTB_Q = as.character(MTB_Q)) %>%
    dplyr::distinct(MTB_Q, .keep_all = TRUE)

  # Load the 1960 urban classification.
  grid_urban <- readRDS(PATH_URBAN) %>%
    dplyr::mutate(MTB_Q = as.character(MTB_Q)) %>%
    dplyr::distinct(MTB_Q, .keep_all = TRUE)

  # Combine geometry, PA coverage and urban classification.
  grid_sop <- grid %>%
    dplyr::left_join(grid_cov %>% dplyr::select(MTB_Q, dplyr::all_of(coverage_cols)), by = "MTB_Q") %>%
    dplyr::left_join(grid_urban %>% dplyr::select(MTB_Q, urban_prop, urban_class), by = "MTB_Q")

  if (any(is.na(grid_sop$urban_class))) {
    stop("Some grid cells are missing an urban classification.")
  }

  # Load the same Germany boundary used in the new analysis.
  de_border <- sf::st_read(PATH_DE, quiet = TRUE, options = "ENCODING=LATIN1") %>%
    dplyr::mutate(BEGINN = as.Date(BEGINN)) %>%
    dplyr::arrange(dplyr::desc(BEGINN)) %>%
    dplyr::slice(1) %>%
    sf::st_transform(25832) %>%
    sf::st_make_valid() %>%
    dplyr::summarise(geometry = sf::st_union(geometry))

  list(grid_sop = grid_sop, de_border = de_border)
}


# 6) Map helper ------------------------------------------------------------------
#
# This follows the active map function from the previous helper.
# Grid-cell boundary lines are intentionally NOT drawn.

make_range_map <- function(map_data, fill_var, legend = FALSE, top_panel = FALSE) {

  grid_sop <- map_data$grid_sop
  de_border <- map_data$de_border

  grid_nonurban <- grid_sop %>% dplyr::filter(urban_class != "urban")
  grid_urban <- grid_sop %>% dplyr::filter(urban_class == "urban")

  bb <- sf::st_bbox(grid_sop)
  urban_legend_dummy <- data.frame(x = bb[["xmin"]], y = bb[["ymin"]], urban = "Urban area > 50%")

  g <- ggplot2::ggplot() +

    ggrastr::rasterise(
      ggplot2::geom_sf(data = grid_nonurban, ggplot2::aes(fill = .data[[fill_var]]), colour = NA, linewidth = 0),
      dpi = 300
    ) +

    ggrastr::rasterise(
      ggplot2::geom_sf(data = grid_urban, fill = urban_col, colour = NA, linewidth = 0),
      dpi = 300
    ) +

    ggrastr::rasterise(
      ggplot2::geom_sf(data = de_border, fill = NA, colour = "grey60", linewidth = 0.3),
      dpi = 200
    ) +

    ggplot2::geom_point(
      data = urban_legend_dummy,
      ggplot2::aes(x = x, y = y, shape = urban),
      inherit.aes = FALSE, alpha = 0
    ) +

    ggplot2::scale_fill_gradient(
      low = "grey95", high = high_col, limits = c(0, 1), breaks = seq(0, 1, by = 0.25),
      labels = scales::label_number(scale = 100, accuracy = 1),
      na.value = "grey95", name = "% of protected grid cell area",
      guide = ggplot2::guide_colorbar(
        direction = "vertical", barheight = grid::unit(8, "cm"), barwidth = grid::unit(0.6, "cm"),
        title.vjust = 4, order = 1
      )
    ) +

    ggplot2::scale_shape_manual(
      values = c("Urban area > 50%" = 22),
      name = NULL,
      guide = ggplot2::guide_legend(
        order = 2,
        override.aes = list(fill = urban_col, colour = urban_col, alpha = 1, size = 8),
        theme = ggplot2::theme(legend.margin = ggplot2::margin(l = -3))
      )
    ) +

    ggplot2::coord_sf(datum = NA, expand = FALSE) +
    ggplot2::theme_void(base_family = base_family_big, base_size = 13) +

    ggplot2::theme(
      text = ggplot2::element_text(family = base_family_big),
      legend.title = ggplot2::element_text(size = if (top_panel) 55 else 50),
      legend.text = ggplot2::element_text(size = 50),
      axis.title = ggplot2::element_blank(),
      axis.text = ggplot2::element_blank(),
      plot.margin = if (top_panel) { ggplot2::margin(5, 12, 30, 5) } else { ggplot2::margin(35, 10, 10, 5) },
      plot.title = ggplot2::element_text(size = 60, family = base_family_big, hjust = 0),
      plot.title.position = "plot"
    )

  if (!legend) {
    g <- g + ggplot2::theme(legend.position = "none")
  }

  g
}


# 7) Scatter-panel helper ---------------------------------------------------------

make_range_panel <- function(pa_object, y_lim, axis_title_size = 50, x_limits = NULL) {

  model_data <- pa_object$data
  pred_df <- pa_object$prediction

  p_val <- get_pa_p(pa_object)
  is_sig <- !is.na(p_val) && p_val < 0.05
  lab_txt <- get_lm_label(pa_object)

  g <- ggplot2::ggplot() +
    ggplot2::geom_point(
      data = model_data, ggplot2::aes(x = range_pct, y = trend),
      colour = pt_col, size = 2.1, alpha = 0.30
    )

  if (is_sig) {

    g <- g +
      ggplot2::geom_ribbon(
        data = pred_df, ggplot2::aes(x = range_pct, ymin = conf_low, ymax = conf_high),
        fill = rib_col, alpha = 0.18, colour = NA
      ) +
      ggplot2::geom_line(
        data = pred_df, ggplot2::aes(x = range_pct, y = predicted),
        colour = line_col, linewidth = 0.8, linetype = "solid"
      )

  } else {

    g <- g +
      ggplot2::geom_line(
        data = pred_df, ggplot2::aes(x = range_pct, y = predicted),
        colour = line_col, linewidth = 0.8, linetype = "dashed"
      )
  }

  g <- g +
    ggplot2::geom_hline(yintercept = 0, colour = "black", linetype = "dashed", linewidth = 0.7) +
    ggplot2::coord_cartesian(ylim = y_lim) +
    ggplot2::annotate(
      "text", x = Inf, y = y_lim[1], label = lab_txt,
      hjust = 1.05, vjust = -0.2, size = 20, family = base_family_big, colour = "grey10"
    ) +
    ggplot2::labs(x = "Protected range [%]", y = "Range retention") +
    ggplot2::theme_classic(base_family = base_family_big) +
    ggplot2::theme(
      axis.text.x = ggplot2::element_text(size = 45),
      axis.text.y = ggplot2::element_text(size = 45),
      axis.title.x = ggplot2::element_text(size = axis_title_size),
      axis.title.y = ggplot2::element_text(size = axis_title_size),
      plot.margin = ggplot2::margin(20, 12, 30, 5)
    )

  if (!is.null(x_limits)) {
    g <- g + ggplot2::xlim(x_limits[1], x_limits[2])
  }

  g
}


# 8) Raincloud helper --------------------------------------------------------------

make_range_raincloud <- function(overall_object, y_lim) {

  model_data <- overall_object$data
  mean_trend <- mean(model_data$trend, na.rm = TRUE)

  ggplot2::ggplot(model_data, ggplot2::aes(x = "All species", y = trend)) +

    see::geom_violinhalf(
      fill = "#5D7388FF", colour = "black", linewidth = 0.7, alpha = 0.7, width = 0.8,
      position = ggplot2::position_nudge(x = 0.08), trim = FALSE
    ) +

    ggplot2::geom_jitter(width = 0.05, alpha = 0.5, size = 2.1, colour = pt_col) +

    ggplot2::geom_boxplot(
      width = 0.15, alpha = 0.5, fill = "#D7C463", colour = "black", outlier.shape = NA,
      position = ggplot2::position_nudge(x = -0.15)
    ) +

    ggplot2::geom_hline(yintercept = 0, colour = "black", linetype = "dashed", linewidth = 0.7) +

    ggplot2::geom_point(
      data = data.frame(x = "All species", trend = mean_trend),
      ggplot2::aes(x = x, y = trend),
      shape = 24, size = 5, fill = "#B0A473", colour = "black",
      position = ggplot2::position_nudge(x = -0.15)
    ) +

    ggplot2::scale_x_discrete(labels = NULL) +
    ggplot2::labs(x = "Density", y = "Range retention") +
    ggplot2::coord_cartesian(clip = "off", ylim = y_lim) +
    ggplot2::theme_classic(base_family = base_family_big) +
    ggplot2::theme(
      axis.text.x = ggplot2::element_blank(),
      axis.title.x = ggplot2::element_text(size = 60),
      axis.text.y = ggplot2::element_text(size = 45),
      axis.title.y = ggplot2::element_text(size = 60),
      plot.margin = ggplot2::margin(20, 12, 30, 5)
    )
}


# 9) Assemble the multipanel figure ------------------------------------------------

build_range_multipanel <- function(
    analysis_id, include_natura = TRUE, top_height = 1.2,
    category_map_height = 1.5, category_plot_height = 0.9) {

  fits <- load_range_models(analysis_id)
  map_data <- load_range_map_data()

  required <- c("All PAs", "Strict", "IUCN IV", "IUCN V")
  if (include_natura) {
    required <- c(required, "Natura 2000")
  }

  missing <- setdiff(required, names(fits))
  if (length(missing) > 0) {
    stop("Missing PA model(s) for ", analysis_id, ": ", paste(missing, collapse = ", "))
  }

  y_lim <- range(fits[["All PAs"]]$data$trend, na.rm = TRUE)

  # Top row: all-PA map, all-PA model and trend distribution.
  map_all <- make_range_map(map_data, fill_var = "cov_anyPA", legend = TRUE, top_panel = TRUE) +
    ggplot2::ggtitle("a) All protected area categories")

  p_overall <- make_range_panel(fits[["All PAs"]], y_lim = y_lim, axis_title_size = 60)
  p_raincloud <- make_range_raincloud(fits[["All PAs"]], y_lim)

  row_top <- (map_all | p_overall | p_raincloud) +
    patchwork::plot_layout(widths = c(1.1, 1.1, 0.6), axes = "collect")

  category_specs <- tibble::tribble(
    ~pa_category,  ~fill_var,                 ~title,        ~subtitle,                 ~x_min, ~x_max,
    "Strict",      "cov_strict",              "Strict protection",  "IUCN I-III + BR cores",      0,     70,
    "IUCN IV",     "cov_IV",                  "IUCN IV",     NA_character_,                 0,     70,
    "IUCN V",      "cov_V",                   "IUCN V",      NA_character_,                 0,     70,
    "Natura 2000", "cov_n2000_habitats_only", "Natura 2000", NA_character_,               NA,     NA
  )

  if (!include_natura) {
    category_specs <- category_specs %>% dplyr::filter(pa_category != "Natura 2000")
  }

  panel_letters <- letters[2:(nrow(category_specs) + 1)]

  maps <- purrr::map2(seq_len(nrow(category_specs)), panel_letters, function(j, letter) {

    spec <- category_specs[j, ]

    g <- make_range_map(map_data, fill_var = spec$fill_var, legend = FALSE, top_panel = FALSE) +
      ggplot2::ggtitle(paste0(letter, ") ", spec$title))

    if (!is.na(spec$subtitle)) {
      g <- g +
        ggplot2::labs(subtitle = spec$subtitle) +
        ggplot2::theme(plot.subtitle = ggplot2::element_text(size = 50, family = base_family_big))
    }

    g
  })

  panels <- purrr::map(seq_len(nrow(category_specs)), function(j) {

    spec <- category_specs[j, ]
    x_limits <- if (is.na(spec$x_min) || is.na(spec$x_max)) { NULL } else { c(spec$x_min, spec$x_max) }

    make_range_panel(fits[[spec$pa_category]], y_lim = y_lim, axis_title_size = 50, x_limits = x_limits)
  })

  n_categories <- nrow(category_specs)

  row_maps <- patchwork::wrap_plots(maps, ncol = n_categories, widths = rep(1, n_categories)) +
    ggplot2::theme(plot.margin = ggplot2::margin(15, 0, 5, 0))

  row_plots <- patchwork::wrap_plots(panels, ncol = n_categories, widths = rep(1, n_categories), axes = "collect") &
    ggplot2::theme(
      axis.text.x = ggplot2::element_text(size = 50),
      axis.title.x = ggplot2::element_text(size = 65),
      axis.text.y = ggplot2::element_text(size = 50),
      axis.title.y = ggplot2::element_text(size = 65)
    )

  plot_categories <- row_maps / row_plots +
    patchwork::plot_layout(ncol = 1, heights = c(category_map_height, category_plot_height))

  full_figure <- row_top / plot_categories +
    patchwork::plot_layout(ncol = 1, heights = c(top_height, 1.5)) &
    ggplot2::theme(plot.margin = ggplot2::margin(8, 8, 8, 8))

  qc <- purrr::map_dfr(required, function(pa) {

    z <- get_pa_result(fits[[pa]])

    tibble::tibble(
      analysis_id = analysis_id,
      pa_category = pa,
      beta = z$beta[[1]],
      p_value = z$p_value[[1]],
      r_squared = z$r_squared[[1]],
      line_style = ifelse(z$p_value[[1]] < 0.05, "solid + ribbon", "dashed")
    )
  })

  list(figure = full_figure, qc = qc)
}
