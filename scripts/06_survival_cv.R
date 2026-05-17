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

message("Loading cohort with scores ...")
df <- readRDS("data/processed/cohort_with_scores.rds")

# Scores valid for CV outcome (skip FINDRISC which is diabetes-specific)
scores <- c("rmrs_score", "b9_score", "pce_score", "framingham_score")

dir.create("results/cache", recursive = TRUE, showWarnings = FALSE)
results <- list()

for (s in scores) {
  message(sprintf("\n--- %s vs CV mortality ---", s))
  rows <- !is.na(df[[s]])
  df_s <- df[rows, ]

  # Fine-Gray for subdistribution hazards
  # Note: get(s) fails inside formula for FGR; use as.formula() instead.
  fgr_formula <- as.formula(sprintf("Hist(followup_years, competing_cv) ~ %s", s))
  fit <- tryCatch(
    FGR(
      fgr_formula,
      data  = df_s,
      cause = 1
    ),
    error = function(e) {
      message(sprintf("  FGR failed: %s", e$message))
      NULL
    }
  )

  # Time-dependent AUC for cause-specific outcome
  roc_obj <- tryCatch(
    timeROC(
      T      = df_s$followup_years,
      delta  = df_s$competing_cv,
      marker = df_s[[s]],
      cause  = 1,
      times  = c(5, 10),
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
    message(sprintf("  n=%d  AUC_cause1(5y)=%.3f  AUC_cause1(10y)=%s",
                    sum(rows), roc_obj$AUC_1[1],
                    ifelse(is.na(roc_obj$AUC_1[2]), "NA", sprintf("%.3f", roc_obj$AUC_1[2]))))
  }
}

saveRDS(results, "results/cache/cv_survival.rds")

summary_tbl <- data.frame(
  score   = scores,
  n       = sapply(results, function(r) r$rows),
  auc_5y  = sapply(results, function(r) if (is.null(r$roc)) NA else r$roc$AUC_1[1]),
  auc_10y = sapply(results, function(r) if (is.null(r$roc)) NA else r$roc$AUC_1[2])
)
write.csv(summary_tbl, "results/cv_summary.csv", row.names = FALSE)

message("\n=== CV mortality summary ===")
print(summary_tbl)
