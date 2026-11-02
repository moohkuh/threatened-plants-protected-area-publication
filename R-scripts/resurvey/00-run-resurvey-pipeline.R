# 00-run-resurvey-pipeline.R
# Run the current ReSurvey Germany publication pipeline.
#
# Run this script from the project root in a fresh R session.
#
# Included analyses:
#   01 Clean ReSurvey data and protected-area overlay
#   02 Colonization and final-extinction models
#      - main analysis excluding Red List category V
#      - sensitivity analysis including category V
#   03 Species-level annual cover-change analysis
#      - paired LMMs
#      - paired Wilcoxon sensitivity
#      - sensitivity excluding Red List category V
#   03b Robust plot-level cover-change models (Table S12)
#   03c Figure S5 (cover change with and without category V)
#   04 Protected-area-size analyses
#   05 Figure 2
#
# Development-only scripts for single-vs-multiple survey intervals,
# distance to protected areas and exploratory model testing are not
# part of the publication pipeline.

RUN_PA_SIZE <- TRUE
RUN_FIGURE_2 <- TRUE

source("./R-scripts/resurvey/01-build-resurvey-data.R")
source("./R-scripts/resurvey/02-fit-event-models.R")
source("./R-scripts/resurvey/03-analyse-cover-change.R")
source("./R-scripts/resurvey/03b-robust-cover-change.R")
source("./R-scripts/resurvey/03c-figure-S5-cover-without-V.R")

if (RUN_PA_SIZE) {
  source("./R-scripts/resurvey/04-analyse-pa-size.R")
}

if (RUN_FIGURE_2) {
  source("./R-scripts/resurvey/05-figure-2.R")
}

cat("\n============================================================\n")
cat("ReSurvey publication pipeline finished.\n")
cat("Run 99-check-resurvey-results.R to check key published values.\n")
cat("============================================================\n")
