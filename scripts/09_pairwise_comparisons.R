#!/usr/bin/env Rscript
# scripts/09_pairwise_comparisons.R
#
# Pre-registered primary pairwise comparisons:
#   1. RMRS vs PCE on CV mortality
#   2. RMRS vs Framingham on CV mortality
#   3. B9 tree vs FINDRISC on diabetes-related mortality
#   4. RMRS vs B9 tree on all-cause mortality
#
# DeLong test on time-dependent AUC at 10y via timeROC::compare.
#
# Inputs:  data/processed/cohort_with_scores.rds
# Outputs: results/pairwise_comparisons.{rds,csv}

.libPaths("/home/po/projects/work/longitudinal-mets-validation/renv/library/R-4.3/x86_64-pc-linux-gnu")

suppressMessages({
  library(timeROC)
  library(survival)
})

df <- readRDS("data/processed/cohort_with_scores.rds")

pairs <- list(
  list(label = "RMRS vs PCE on CV mortality",
       score1 = "rmrs_score", score2 = "pce_score",
       outcome_col = "event_cv"),
  list(label = "RMRS vs Framingham on CV mortality",
       score1 = "rmrs_score", score2 = "framingham_score",
       outcome_col = "event_cv"),
  list(label = "B9 vs FINDRISC on diabetes mortality",
       score1 = "b9_score", score2 = "findrisc_score",
       outcome_col = "event_dm"),
  list(label = "RMRS vs B9 on all-cause mortality",
       score1 = "rmrs_score", score2 = "b9_score",
       outcome_col = "event_allcause")
)

results <- list()

for (p in pairs) {
  message(sprintf("\n--- %s ---", p$label))
  rows <- !is.na(df[[p$score1]]) & !is.na(df[[p$score2]])
  df_p <- df[rows, ]

  roc1 <- tryCatch(
    timeROC(T = df_p$followup_years, delta = df_p[[p$outcome_col]],
            marker = df_p[[p$score1]], cause = 1, times = 10, iid = TRUE),
    error = function(e) NULL
  )
  roc2 <- tryCatch(
    timeROC(T = df_p$followup_years, delta = df_p[[p$outcome_col]],
            marker = df_p[[p$score2]], cause = 1, times = 10, iid = TRUE),
    error = function(e) NULL
  )

  if (is.null(roc1) || is.null(roc2)) {
    message("  Skipping pair due to timeROC failure")
    next
  }

  delong <- tryCatch(
    compare(roc1, roc2, abseps = 1e-06),
    error = function(e) NULL
  )

  results[[p$label]] <- list(
    pair     = p,
    n        = sum(rows),
    auc1     = roc1$AUC[1],
    auc2     = roc2$AUC[1],
    delong_p = if (is.null(delong)) NA else delong$p_values_AUC[2]
  )

  message(sprintf("  n=%d  AUC(%s)=%.3f  AUC(%s)=%.3f  DeLong p=%.4g",
                  sum(rows), p$score1, roc1$AUC[1],
                  p$score2, roc2$AUC[1],
                  if (is.null(delong)) NA else delong$p_values_AUC[2]))
}

dir.create("results", showWarnings = FALSE)
saveRDS(results, "results/pairwise_comparisons.rds")

summary_tbl <- do.call(rbind, lapply(names(results), function(nm) {
  r <- results[[nm]]
  data.frame(
    comparison = nm,
    n          = r$n,
    auc_score1 = r$auc1,
    auc_score2 = r$auc2,
    delong_p   = r$delong_p
  )
}))
write.csv(summary_tbl, "results/pairwise_comparisons.csv", row.names = FALSE)

message("\n=== Pairwise comparisons ===")
print(summary_tbl)
