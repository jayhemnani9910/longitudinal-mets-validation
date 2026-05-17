#!/usr/bin/env Rscript
# scripts/04_compute_scores.R
#
# Computes all 5 risk scores on the analysis cohort:
#   RMRS, B9 decision tree, PCE, Framingham 2008, FINDRISC.
#
# Output: data/processed/cohort_with_scores.rds
#
# Run after scripts/03_apply_inclusion.R has produced
# data/processed/analysis_cohort.rds.

.libPaths("/home/po/projects/work/longitudinal-mets-validation/renv/library/R-4.3/x86_64-pc-linux-gnu")

suppressMessages({
  library(dplyr)
})

# Source order matters: rmrs first (defines elliot_sigmoid + scale_factor)
source("R/scores/rmrs.R")
source("R/scores/b9_features.R")
source("R/scores/b9_tree.R")
source("R/scores/pce.R")
source("R/scores/framingham.R")
source("R/scores/findrisc.R")

message("Loading analysis cohort ...")
df <- readRDS("data/processed/analysis_cohort.rds")
message(sprintf("  N = %d", nrow(df)))

# Some PCE inputs may be NA for subjects who only have non-fasting cholesterol.
# Note these losses but compute anyway.
n_missing_tchol <- sum(is.na(df$total_chol))
message(sprintf("  Subjects without total cholesterol (PCE will be NA): %d",
                n_missing_tchol))

# ---- RMRS ------------------------------------------------------------------
message("Computing RMRS ...")
df$rmrs_score <- mapply(
  rmrs,
  df$waist_cm, df$sbp, df$dbp, df$hdl, df$fasting_glucose,
  df$triglycerides, as.character(df$sex),
  MoreArgs = list(ethnicity = "us")
)

# ---- B9 decision tree -------------------------------------------------------
message("Computing B9 tree predictions ...")
df$b9_features <- mapply(
  b9_features,
  df$waist_cm, df$sbp, df$dbp, as.character(df$sex),
  MoreArgs = list(ethnicity = "us"),
  SIMPLIFY = FALSE
)
df$BPWC_add <- sapply(df$b9_features, `[[`, "BPWC_add")
df$BPWC_mul <- sapply(df$b9_features, `[[`, "BPWC_mul")
df$BPWC_dif <- sapply(df$b9_features, `[[`, "BPWC_dif")
df$b9_features <- NULL  # drop list-col after extraction

df$b9_score <- mapply(b9_tree_predict, df$BPWC_add, df$BPWC_mul, df$BPWC_dif)

# ---- PCE -------------------------------------------------------------------
message("Computing PCE ...")
# PCE supports white and black only; assign white as default approximation
# for Hispanic, Asian, Other per ACC/AHA app convention. Document in manuscript.
df$pce_race <- case_when(
  df$race == "NHBlack" ~ "black",
  TRUE                 ~ "white"
)
# PCE only valid for 40-79; assign NA outside
df$pce_score <- ifelse(
  df$age >= 40 & df$age <= 79 & !is.na(df$total_chol),
  mapply(
    pce,
    df$age, as.character(df$sex), df$pce_race,
    df$total_chol, df$hdl, df$sbp,
    df$bp_treatment, df$current_smoker, df$prior_t2d
  ),
  NA_real_
)

# ---- Framingham 2008 -------------------------------------------------------
message("Computing Framingham 2008 ...")
df$framingham_score <- ifelse(
  !is.na(df$total_chol),
  mapply(
    framingham_2008,
    df$age, as.character(df$sex),
    df$total_chol, df$hdl, df$sbp,
    df$bp_treatment, df$current_smoker, df$prior_t2d
  ),
  NA_real_
)

# ---- FINDRISC --------------------------------------------------------------
message("Computing FINDRISC ...")
# Use simple proxies for NHANES variables not directly captured:
#   physical_active: would come from PAQ; for now FALSE (most NHANES adults)
#   daily_vegetables: would come from DR1TOT; for now FALSE
#   family_history_dm: would come from MCQ300C; for now "none"
# These proxies will be refined in a future iteration (Phase 6 sensitivity).
df$physical_active <- FALSE
df$daily_vegetables <- FALSE
df$family_history_dm <- "none"
df$prior_high_glucose <- df$prior_t2d  # use prior T2D as proxy

df$findrisc_score <- mapply(
  findrisc,
  df$age, df$bmi, df$waist_cm, as.character(df$sex),
  df$physical_active, df$daily_vegetables,
  df$bp_treatment, df$prior_high_glucose,
  df$family_history_dm
)

# ---- Save + summary --------------------------------------------------------
saveRDS(df, "data/processed/cohort_with_scores.rds")

message("\n=== Score distributions ===")
for (s in c("rmrs_score", "b9_score", "pce_score",
            "framingham_score", "findrisc_score")) {
  vals <- df[[s]]
  vals <- vals[!is.na(vals)]
  message(sprintf("  %-20s n=%6d  median=%.4f  IQR=[%.4f, %.4f]  range=[%.4f, %.4f]",
                  s, length(vals), median(vals),
                  quantile(vals, 0.25), quantile(vals, 0.75),
                  min(vals), max(vals)))
}

message(sprintf("\n=== Saved: data/processed/cohort_with_scores.rds (%d subjects) ===",
                nrow(df)))
