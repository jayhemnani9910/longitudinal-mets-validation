#!/usr/bin/env Rscript
# scripts/07_survival_diabetes.R
#
# Diabetes-related mortality survival analysis. Same Fine-Gray approach as CV,
# but scores: RMRS, B9 tree, FINDRISC (skip PCE / Framingham since those
# are CVD-specific).
#
# Inputs:  data/processed/cohort_with_scores.rds
# Outputs: results/cache/dm_survival.rds, results/dm_summary.csv

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

scores <- c("rmrs_score", "b9_score", "findrisc_score")

dir.create("results/cache", recursive = TRUE, showWarnings = FALSE)
results <- list()

for (s in scores) {
  message(sprintf("\n--- %s vs diabetes mortality ---", s))
  rows <- !is.na(df[[s]])
  df_s <- df[rows, ]

  # Diabetes outcome uses the 15y follow-up frame and broadened cause
  # definition built in script 03 (followup_years_dm + competing_dm).
  fgr_formula <- as.formula(sprintf("Hist(followup_years_dm, competing_dm) ~ %s", s))
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

  # Time-dependent AUC at 5y, 10y, and the 15y horizon operationalized as
  # t=14.5 to avoid IPCW degeneracy at the 15y follow-up cap.
  roc_obj <- tryCatch(
    timeROC(
      T      = df_s$followup_years_dm,
      delta  = df_s$competing_dm,
      marker = df_s[[s]],
      cause  = 1,
      times  = c(5, 10, 14.5),
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
    fmt_auc <- function(x) ifelse(is.na(x), "NA", sprintf("%.3f", x))
    message(sprintf("  n=%d  AUC(5y)=%s  AUC(10y)=%s  AUC(14.5y)=%s",
                    sum(rows),
                    fmt_auc(roc_obj$AUC_1[1]),
                    fmt_auc(roc_obj$AUC_1[2]),
                    fmt_auc(roc_obj$AUC_1[3])))
  }
}

saveRDS(results, "results/cache/dm_survival.rds")

summary_tbl <- data.frame(
  score     = scores,
  n         = sapply(results, function(r) r$rows),
  auc_5y    = sapply(results, function(r) if (is.null(r$roc)) NA else r$roc$AUC_1[1]),
  auc_10y   = sapply(results, function(r) if (is.null(r$roc)) NA else r$roc$AUC_1[2]),
  auc_14_5y = sapply(results, function(r) if (is.null(r$roc)) NA else r$roc$AUC_1[3])
)
write.csv(summary_tbl, "results/dm_summary.csv", row.names = FALSE)

message("\n=== Diabetes mortality summary ===")
print(summary_tbl)
