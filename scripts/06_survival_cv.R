#!/usr/bin/env Rscript
# scripts/06_survival_cv.R
#
# CV mortality survival analysis. Fine-Gray subdistribution hazards model for
# competing risks, time-dependent AUC at 5 and 10 years.
#
# Inputs:  data/processed/cohort_with_scores.rds
# Outputs: results/cache/cv_survival.rds, results/cv_summary.csv

.libPaths("/home/po/projects/work/longitudinal-mets-validation/renv/library/R-4.3/x86_64-pc-linux-gnu")

suppressMessages({
  library(survival)
  library(riskRegression)
  library(cmprsk)
  library(timeROC)
  library(prodlim)
  library(dplyr)
})

source("R/utils/survival_helpers.R")

message("Loading cohort with scores ...")
df <- readRDS("data/processed/cohort_with_scores.rds")

# Cause-specific CV mortality is restricted to cycles with full leading-cause
# coding (1999-2014); the 2015-2018 public-use file folds stroke into "all
# other", so those cycles are excluded from the CV analysis.
df <- df[df$cause_coded, ]
df$cluster_id <- paste(df$sdmvstra, df$sdmvpsu, sep = "_")
message(sprintf("  Cause-coded cohort N = %d", nrow(df)))

# Scores valid for CV outcome (skip FINDRISC which is diabetes-specific)
scores <- c("rmrs_score", "b9_score", "pce_score", "framingham_score")

dir.create("results/cache", recursive = TRUE, showWarnings = FALSE)
results <- list()

for (s in scores) {
  message(sprintf("\n--- %s vs CV mortality ---", s))
  rows <- !is.na(df[[s]])
  df_s <- df[rows, ]

  # Survey-weighted Fine-Gray subdistribution hazard ratio (per-SD-free, raw
  # score scale). finegray() + weighted Cox carries the MEC weight; robust SE
  # is clustered on the PSU.
  fg <- weighted_finegray_hr(df_s, s, "followup_years", "competing_cv",
                             "wt_mec", "cluster_id", cause = 1)
  fit <- fg

  # Time-dependent AUC for cause-specific outcome. 10y horizon operationalized
  # as t=9.5 to avoid IPCW degeneracy at the follow-up cap (see script 05 note).
  roc_obj <- tryCatch(
    timeROC(
      T      = df_s$followup_years,
      delta  = df_s$competing_cv,
      marker = df_s[[s]],
      cause  = 1,
      times  = c(5, 9.5),
      weighting = "marginal",
      iid    = FALSE   # iid=TRUE OOM; bootstrap CIs deferred
    ),
    error = function(e) {
      message(sprintf("  timeROC failed: %s", e$message))
      NULL
    }
  )

  results[[s]] <- list(rows = sum(rows), fit = fit, roc = roc_obj)

  if (!is.null(roc_obj)) {
    # For ipcwcompetingrisksROC: AUC_1[t1,t2] = AUC for cause 1 at each timepoint
    message(sprintf("  n=%d  AUC_cause1(5y)=%.3f  AUC_cause1(9.5y)=%s  wFG HR=%s",
                    sum(rows), roc_obj$AUC_1[1],
                    ifelse(is.na(roc_obj$AUC_1[2]), "NA", sprintf("%.3f", roc_obj$AUC_1[2])),
                    ifelse(is.na(fg$hr), "NA", sprintf("%.3f (%.3f, %.3f)", fg$hr, fg$lo, fg$hi))))
  }
}

saveRDS(results, "results/cache/cv_survival.rds")

summary_tbl <- data.frame(
  score    = scores,
  n        = sapply(results, function(r) r$rows),
  auc_5y   = sapply(results, function(r) if (is.null(r$roc)) NA else r$roc$AUC_1[1]),
  auc_9_5y = sapply(results, function(r) if (is.null(r$roc)) NA else r$roc$AUC_1[2]),
  wfg_hr   = sapply(results, function(r) if (is.null(r$fit)) NA else r$fit$hr),
  wfg_lo   = sapply(results, function(r) if (is.null(r$fit)) NA else r$fit$lo),
  wfg_hi   = sapply(results, function(r) if (is.null(r$fit)) NA else r$fit$hi)
)
write.csv(summary_tbl, "results/cv_summary.csv", row.names = FALSE)

message("\n=== CV mortality summary ===")
print(summary_tbl)
