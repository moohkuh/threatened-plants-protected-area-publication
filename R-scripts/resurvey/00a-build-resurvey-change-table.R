# 00a-build-resurvey-change-table.R
# Build rs_raw from the original ReSurvey Germany data.
#
# SOURCE
#   This script is based on the code provided with ReSurvey Germany
#   and Jandt et al. (2022).
#
#   The procedure for deriving species-level changes between subsequent
#   vegetation surveys is taken from the ReSurvey Germany analysis code
#   and is not original code developed for the present study.
#
#   Only file paths and the output location were adapted for the
#   reproducible workflow of the present publication.
#
# REFERENCE
#   Jandt et al. (2022). More losses than gains during one century
#   of plant biodiversity change in Germany. Nature.
#
# OUTPUT
#   Species x plot x survey-interval change table used as input
#   for the protected-area analyses.


# 0) Packages ---------------------------------------------------------------

library(data.table)
library(reshape2)
library(tibble)
library(readr)


# 1) Paths --------------------------------------------------------------------

PATH_RS_VEGETATION <- "./data/raw/resurvey/ReSurveyGermany.csv"
PATH_RS_HEADER_PUBLIC <- "./data/raw/resurvey/Header_ReSurveyGermany.csv"
PATH_RS_RAW <- "./data/derived/resurvey/RS_raw.csv"

dir.create(dirname(PATH_RS_RAW), recursive = TRUE, showWarnings = FALSE)


# 2) Function to combine cover across vegetation layers -----------------------

combine_cover <- function(x) {
  # Combine cover values assuming random overlap between layers
  while (length(x) > 1) {
    x[2] <- x[1] + (100 - x[1]) * x[2] / 100
    x <- x[-1]
  }
  x
}


# 3) Load ReSurvey Germany data ------------------------------------------------

DT <- fread(PATH_RS_VEGETATION)
header <- fread(PATH_RS_HEADER_PUBLIC)

# Create a cross-link table between project names and PROJECT_ID
Projects <- unique(header[, .(PROJECT_ID, RS_PROJECT)])


# 4) Combine cover across vegetation layers ------------------------------------

DT2 <- DT[, .(Cover_Perc = combine_cover(Cover_Perc)), by = .(PROJECT_ID_RELEVE_NR, TaxonName)]

# Split PROJECT_ID_RELEVE_NR into PROJECT_ID and RELEVE_NR
DT2[, PROJECT_ID := as.numeric(tstrsplit(PROJECT_ID_RELEVE_NR, "_")[[1]])]
DT2[, RELEVE_NR := as.numeric(tstrsplit(PROJECT_ID_RELEVE_NR, "_")[[2]])]


# 5) Build species-change table ------------------------------------------------

out_list <- vector("list", length = 0L)
idx <- 0L
project_ids <- sort(unique(DT2$PROJECT_ID))

for (pid in project_ids) {

  message("Processing PROJECT_ID = ", pid)

  # Species cover data for this project
  species_long <- DT2[PROJECT_ID == pid]
  if (nrow(species_long) == 0) next

  # Header data for this project
  env <- header[PROJECT_ID == pid]
  if (nrow(env) == 0) next

  # Create RELEVE_NR x TaxonName matrix
  species_mat <- reshape2::acast(species_long, RELEVE_NR ~ TaxonName, value.var = "Cover_Perc", fill = 0)

  # Keep only releves that occur in the header
  common_releves <- base::intersect(as.numeric(rownames(species_mat)), env$RELEVE_NR)
  if (length(common_releves) == 0) next

  species_mat <- species_mat[as.character(common_releves), , drop = FALSE]

  # Match ordering between species matrix and header
  env2 <- env[match(as.numeric(rownames(species_mat)), env$RELEVE_NR)]
  RS_PLOT <- env2$RS_PLOT
  YEAR <- env2$YEAR

  # Only time series with more than one record
  plot_list <- sort(names(which(table(RS_PLOT) > 1)))
  if (length(plot_list) == 0) next

  # Project name used in output
  proj_name <- Projects[PROJECT_ID == pid, RS_PROJECT][1]


  # 6) Loop over plot time series -------------------------------------------

  for (pl in plot_list) {

    sel <- which(RS_PLOT == pl)
    target_year <- YEAR[sel]

    # At least two different survey years are required
    if (length(unique(target_year)) <= 1) next

    sp_tp <- species_mat[sel, , drop = FALSE]

    # Remove species absent throughout the time series
    sp_tp <- sp_tp[, colSums(sp_tp) != 0, drop = FALSE]
    if (ncol(sp_tp) == 0) next

    # Average replicate releves within the same year
    sp_tp_dt <- as.data.table(sp_tp)
    sp_tp_dt[, year_group := target_year]
    sp_year <- sp_tp_dt[, lapply(.SD, mean), by = year_group, .SDcols = setdiff(names(sp_tp_dt), "year_group")]

    sp_year[, year_num := as.numeric(year_group)]

    # Order from latest to earliest survey
    setorder(sp_year, -year_num)

    years_sorted <- sp_year$year_num
    sp_year_mat <- as.matrix(sp_year[, !c("year_group", "year_num")])

    # Relative cover scaled by the maximum cover of each species
    max_by_species <- apply(sp_year_mat, 2, max)
    sp_rel <- sweep(sp_year_mat, 2, max_by_species, FUN = "/")

    # Relative abundance within each survey
    sp_rel2 <- sweep(sp_year_mat, 1, rowSums(sp_year_mat), FUN = "/")
    sp_rel2[is.nan(sp_rel2)] <- 0

    # Relative species ranks within each survey
    if (ncol(sp_rel2) == 1) {
      # With only one species, its rank is always 1
      ranks <- matrix(1, nrow = nrow(sp_rel2), ncol = 1)
    } else {
      ranks <- t(apply(-sp_rel2, 1, rank))
    }
    rel_ranks <- sweep(ranks, 1, apply(ranks, 1, max), FUN = "/")


    # 7) Compare subsequent survey intervals ---------------------------------

    for (k in 1:(nrow(sp_year_mat) - 1)) {

      to_year <- years_sorted[k]
      from_year <- years_sorted[k + 1]

      # Absolute change
      abs_change <- sp_year_mat[k, ] - sp_year_mat[k + 1, ]
      abs_change[sp_year_mat[k, ] == 0 & sp_year_mat[k + 1, ] == 0] <- NA_real_

      # Relative change
      rel_change <- sp_rel[k, ] - sp_rel[k + 1, ]
      rel_change[sp_rel[k, ] == 0 & sp_rel[k + 1, ] == 0] <- NA_real_

      # Relative rank change
      rr_change <- rel_ranks[k, ] - rel_ranks[k + 1, ]
      rr_change[rel_ranks[k, ] == 0 & rel_ranks[k + 1, ] == 0] <- NA_real_

      # Colonizers: present later, absent earlier
      abs_col <- abs_change
      abs_col[sp_year_mat[k + 1, ] != 0] <- NA_real_

      # Extinct: present earlier, absent later
      abs_ext <- abs_change
      abs_ext[sp_year_mat[k, ] != 0] <- NA_real_

      # Append interval to output
      idx <- idx + 1L

      out_list[[idx]] <- data.table(
        RS_PROJECT = proj_name,
        RS_PLOT = pl,
        from.n = as.numeric(sum(RS_PLOT == pl & YEAR == from_year)),
        to.n = as.numeric(sum(RS_PLOT == pl & YEAR == to_year)),
        from = from_year,
        to = to_year,
        species = colnames(sp_year_mat),
        absolute.change = as.numeric(abs_change),
        relative.change = as.numeric(rel_change),
        relative.rank.change = as.numeric(rr_change),
        absolute.change.colonizer = as.numeric(abs_col),
        absolute.change.extinct = as.numeric(abs_ext)
      )
    }
  }
}


# 8) Combine results -----------------------------------------------------------

rs_raw <- rbindlist(out_list, use.names = TRUE, fill = TRUE)

# Keep the same filtering as in the original workflow
rs_raw <- rs_raw[!is.na(absolute.change)]

# Add combined project-plot identifier
rs_raw[, RS_PROJECT_PLOT := paste(RS_PROJECT, RS_PLOT, sep = "_")]

# Convert to tibble
rs_raw <- as_tibble(rs_raw)


# 9) Save derived ReSurvey table ------------------------------------------------

write_csv(rs_raw, PATH_RS_RAW)
