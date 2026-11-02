# 05c-range-size-protection-figure.R
# Baseline range size vs. protected range share (Figure S2).
# Reproduces the original Figure S2 aesthetics, but fixes the old
# construction bug so that Natura 2000 is actually included as the
# fifth facet.

source("./R-scripts/range/05a-range-figure-helpers.R")

fits <- load_range_models("T1_T3_main")

pa_order <- c("All PAs", "Strict", "IUCN IV", "IUCN V", "Natura 2000")

missing <- setdiff(pa_order, names(fits))
if (length(missing) > 0) {
  stop("Missing T1-T3 model(s): ", paste(missing, collapse = ", "))
}

cor_plot_data <- purrr::map_dfr(pa_order, function(pa) {
  fits[[pa]]$data %>%
    dplyr::transmute(pa_cat = pa, range_pct = range_pct, log_SOP_T1 = log_SOP_start)
}) %>%
  dplyr::mutate(pa_cat = factor(pa_cat, levels = pa_order))

# Spearman rho per category, printed as a quick sanity check against the text
spearman_tbl <- cor_plot_data %>%
  dplyr::group_by(pa_cat) %>%
  dplyr::summarise(
    rho = stats::cor(range_pct, log_SOP_T1, method = "spearman", use = "complete.obs"),
    .groups = "drop"
  )

print(spearman_tbl)

p_cor <- ggplot2::ggplot(cor_plot_data, ggplot2::aes(x = log_SOP_T1, y = range_pct)) +
  ggplot2::geom_point(alpha = 0.25, size = 1.2, colour = "#2F436C") +
  ggplot2::geom_smooth(
    method = "loess", formula = y ~ x,
    colour = "#00336F", fill = "#E3CE5C", linewidth = 0.8, alpha = 0.2
  ) +
  ggplot2::facet_wrap(~ pa_cat, ncol = 5) +
  ggplot2::labs(x = "Baseline range size (log SOP₁)", y = "Protected range share [%]") +
  ggplot2::theme_bw(base_size = 11, base_family = "roboto_condensed") +
  ggplot2::theme(
    strip.text       = ggplot2::element_text(size = 30),
    strip.background = ggplot2::element_rect(fill = "grey96", colour = "grey80"),
    axis.text        = ggplot2::element_text(size = 25),
    axis.title       = ggplot2::element_text(size = 30)
  )

ggplot2::ggsave(
  filename = file.path(RANGE_FIG_DIR, "FigureS2_range_size_protection.png"),
  plot = p_cor,
  width = 14, height = 8, dpi = 200
)
