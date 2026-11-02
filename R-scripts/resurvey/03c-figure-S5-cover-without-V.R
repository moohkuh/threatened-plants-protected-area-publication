# 03c-figure-S5-cover-without-V.R
# Figure S5: inside-outside contrasts of species-level cover change,
# all Red List categories vs. without near-threatened species (V).
# Needs the cache written by 03-analyse-cover-change.R.
# Strict is not shown because the model without V did not converge.


# 0) Setup ---------------------------------------------------------------

source("./R-scripts/resurvey/00-preamble-resurvey.R")

RDS_FILE <- file.path(OUT_DIR_RDS, "species_cover_means_objects.rds")
assert_files_exist(RDS_FILE)

cache <- readRDS(RDS_FILE)


# 1) Combine main and without-V contrasts --------------------------------

# number of species in the main models
species_means_main <- list(
  "All PAs"     = cache$sp_all,
  "IUCN IV"     = cache$sp_iucn4,
  "IUCN V"      = cache$sp_iucn5,
  "Natura 2000" = cache$sp_n2k
)

effects_main <- cache$effects_tbl %>%
  dplyr::mutate(scenario = as.character(scenario)) %>%
  dplyr::filter(scenario != "Strict")

effects_main$n_species <- NA_integer_
for (i in seq_len(nrow(effects_main))) {
  sp <- species_means_main[[effects_main$scenario[i]]]
  effects_main$n_species[i] <- dplyr::n_distinct(sp$wcvp_name)
}
effects_main$subset <- "All red list categories"

effects_no_v <- cache$effects_no_v %>%
  dplyr::mutate(
    scenario = as.character(scenario),
    subset = "Without near threatened species (V)"
  )

plot_data <- dplyr::bind_rows(
  effects_main %>% dplyr::select(scenario, subset, n_species, estimate, lower, upper, p_value),
  effects_no_v %>% dplyr::select(scenario, subset, n_species, estimate, lower, upper, p_value)
) %>%
  dplyr::mutate(
    scenario = dplyr::recode(scenario, "All PAs" = "All"),
    # last factor level is drawn at the top
    scenario = factor(scenario, levels = c("Natura 2000", "IUCN V", "IUCN IV", "All")),
    subset = factor(subset, levels = c("All red list categories",
                                       "Without near threatened species (V)")),
    sig_flag = dplyr::if_else(p_value < 0.05, "significant", "not significant")
  )

print(plot_data)

# position of the sample-size labels: right of the widest CI
x_range <- range(c(plot_data$lower, plot_data$upper), na.rm = TRUE)
plot_data$label_x <- x_range[2] + 0.15 * diff(x_range)


# 2) Figure --------------------------------------------------------------

cols_scenario <- c(
  "All"         = col_allpa,
  "IUCN IV"     = col_iucn4,
  "IUCN V"      = col_iucn5,
  "Natura 2000" = col_n2k
)

dodge <- ggplot2::position_dodge(width = 0.5)

fig_s5 <- ggplot2::ggplot(
  plot_data,
  ggplot2::aes(
    y = scenario, x = estimate,
    color = scenario, group = subset, shape = subset, linetype = sig_flag
  )
) +
  ggplot2::geom_vline(xintercept = 0, linetype = "dashed",
                      color = "grey50", linewidth = 0.6) +
  ggplot2::geom_errorbar(
    ggplot2::aes(xmin = lower, xmax = upper),
    orientation = "y", height = 0.2, linewidth = 0.7, position = dodge
  ) +
  ggplot2::geom_point(fill = "white", size = 3, stroke = 2, position = dodge) +
  ggplot2::geom_text(
    ggplot2::aes(x = label_x, label = paste0("n = ", n_species)),
    position = dodge, hjust = 0.4, vjust = 0.54,
    size = 13, color = "grey20", show.legend = FALSE
  ) +
  ggplot2::scale_x_continuous(expand = ggplot2::expansion(mult = c(0.05, 0.1))) +
  ggplot2::scale_color_manual(values = cols_scenario, guide = "none") +
  ggplot2::scale_shape_manual(
    values = c("All red list categories" = 21,
               "Without near threatened species (V)" = 24),
    name = NULL
  ) +
  ggplot2::scale_linetype_manual(
    values = c("significant" = "solid", "not significant" = "32"),
    guide = "none"
  ) +
  ggplot2::labs(x = "Effect inside vs. outside", y = NULL) +
  ggplot2::theme_classic(base_family = base_family) +
  ggplot2::theme(
    legend.position = "bottom",
    legend.text = ggplot2::element_text(size = 40),
    axis.text.x = ggplot2::element_text(size = 40),
    axis.text.y = ggplot2::element_text(size = 40),
    axis.title = ggplot2::element_text(size = 45)
  )

ggplot2::ggsave(
  file.path(OUT_FIG_RS_DIR, "FigureS5_cover_change_without_V.png"),
  plot = fig_s5,
  width = 10, height = 7, units = "in", dpi = 300
)
