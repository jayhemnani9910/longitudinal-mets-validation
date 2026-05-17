#!/usr/bin/env Rscript
# scripts/05_survival_allcause.R
#
# All-cause mortality survival analysis. Cox proportional hazards model with
# survey weights, time-dependent ROC AUC at 5 and 10 years (IPCW-corrected),
# and Brier score.
#
# Inputs:  data/processed/cohort_with_scores.rds
# Outputs: results/cache/allcause_survival.rds, results/allcause_summary.csv

.libPaths("/home/po/projects/work/longitudinal-mets-validation/renv/library/R-4.3/x86_64-pc-linux-gnu")

suppressMessages({
  library(survey)
  library(survival)
  library(timeROC)
  library(riskRegression)
  library(dplyr)
})

source("R/utils/survey_design.R")

message("Loading cohort with scores ...")
df <- readRDS("data/processed/cohort_with_scores.rds")
des <- build_survey_design(df)

scores <- c("rmrs_score", "b9_score", "pce_score",
            "framingham_score", "findrisc_score")

dir.create("results/cache", recursive = TRUE, showWarnings = FALSE)
results <- list()

for (s in scores) {
  message(sprintf("\n--- %s ---", s))
  rows <- !is.na(df[[s]])
  df_s <- df[rows, ]

  # Cox model with survey weights
  fit <- tryCatch(
    svycoxph(
      as.formula(sprintf("Surv(followup_years, event_allcause) ~ %s", s)),
      design = build_survey_design(df_s)
    ),
    error = function(e) {
      message(sprintf("  svycoxph failed: %s", e$message))
      NULL
    }
  )

  # Time-dependent ROC at 5 and 10 years (IPCW)
  roc_obj <- tryCatch(
    timeROC(
      T      = df_s$followup_years,
      delta  = df_s$event_allcause,
      marker = df_s[[s]],
      cause  = 1,
      times  = c(5, 10),
      weighting = "marginal",
      iid    = TRUE
    ),
    error = function(e) {
      message(sprintf("  timeROC failed: %s", e$message))
      NULL
    }
  )

  results[[s]] <- list(rows = sum(rows), fit = fit, roc = roc_obj)

  if (!is.null(roc_obj)) {
    message(sprintf("  n=%d  AUC(5y)=%.3f  AUC(10y)=%.3f",
                    sum(rows), roc_obj$AUC[1], roc_obj$AUC[2]))
  }
}

saveRDS(results, "results/cache/allcause_survival.rds")

# Summary CSV
summary_tbl <- data.frame(
  score   = scores,
  n       = sapply(results, function(r) r$rows),
  auc_5y  = sapply(results, function(r) if (is.null(r$roc)) NA else r$roc$AUC[1]),
  auc_10y = sapply(results, function(r) if (is.null(r$roc)) NA else r$roc$AUC[2])
)
write.csv(summary_tbl, "results/allcause_summary.csv", row.names = FALSE)

message("\n=== All-cause mortality summary ===")
print(summary_tbl)
