# 05d-figureS3-withV.R
# Sensitivity check: T1 -> T3 including Red List category V (Figure S3).

source("./R-scripts/range/05a-range-figure-helpers.R")

fig_obj <- build_range_multipanel(
  analysis_id = "T1_T3_with_V_sensitivity",
  include_natura = TRUE,
  top_height = 1.3
)

print(fig_obj$qc)

ggplot2::ggsave(
  filename = file.path(RANGE_FIG_DIR, "FigureS3_T1_T3_with_V.png"),
  plot = fig_obj$figure,
  width = 17, height = 17, dpi = 280
)
