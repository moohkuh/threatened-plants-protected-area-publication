# 05b-figure1-range-main.R
# Main T1 -> T3 range figure (Figure 1).

source("./R-scripts/range/05a-range-figure-helpers.R")

fig_main <- build_range_multipanel(
  analysis_id = "T1_T3_main",
  include_natura = TRUE,
  top_height = 1.3,
  category_map_height = 1.5,
  category_plot_height = 0.9
)

print(fig_main$qc)

ggplot2::ggsave(
  filename = file.path(RANGE_FIG_DIR, "Figure1_range_T1_T3.png"),
  plot = fig_main$figure,
  width = 17, height = 17, units = "in", dpi = 280
)
