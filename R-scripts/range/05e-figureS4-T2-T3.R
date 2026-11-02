# 05e-figureS4-T2-T3.R
# Temporal sensitivity check: T2 -> T3 (Figure S4).
# Uses the species set excluding Red List category V, consistent with the
# main analysis 

source("./R-scripts/range/05a-range-figure-helpers.R")

fig_obj <- build_range_multipanel(
  analysis_id = "T2_T3_temporal",
  include_natura = TRUE,
  top_height = 1.3
)

print(fig_obj$qc)

# Save as png.
ggplot2::ggsave(
  filename = file.path(RANGE_FIG_DIR, "FigureS4_T2_T3.png"),
  plot = fig_obj$figure,
  width = 17, height = 17, dpi = 280
)
