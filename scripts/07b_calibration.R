#!/usr/bin/env Rscript
# scripts/07b_calibration.R
#
# Calibration plots and Brier scores per (outcome, score). Uses
# riskRegression::Score with IPCW correction.
#
# Inputs:  data/processed/cohort_with_scores.rds
# Outputs: results/calibration/cal_<outcome>_<score>.png, results/calibration_summary.csv

.libPaths("/home/po/projects/work/longitudinal-mets-validation/renv/library/R-4.3/x86_64-pc-linux-gnu")

suppressMessages({
  library(survival)
  library(riskRegression)
  library(prodlim)
  library(dplyr)
  library(ggplot2)
})

df <- readRDS("data/processed/cohort_with_scores.rds")

dir.create("results/calibration", recursive = TRUE, showWarnings = FALSE)

# Cause-specific outcomes use the 1999-2014 cause-coded cohort and the same
# off-cap horizons as the survival scripts (9.5y for all-cause/CV, 14.5y for
# the 15y diabetes frame). Diabetes uses the followup_years_dm time variable.
df_cause <- df[df$cause_coded, ]

outcomes <- list(
  allcause = list(formula = Surv(followup_years, event_allcause) ~ 1,
                  data = df, time = 9.5,
                  scores  = c("rmrs_score", "b9_score", "pce_score",
                              "framingham_score", "findrisc_score")),
  cv       = list(formula = Hist(followup_years, competing_cv) ~ 1,
                  cause = 1, data = df_cause, time = 9.5,
                  scores = c("rmrs_score", "b9_score", "pce_score",
                             "framingham_score")),
  dm       = list(formula = Hist(followup_years_dm, competing_dm) ~ 1,
                  cause = 1, data = df_cause, time = 14.5,
                  scores = c("rmrs_score", "b9_score", "findrisc_score"))
)

cal_rows <- list()
for (out_name in names(outcomes)) {
  out <- outcomes[[out_name]]
  for (s in out$scores) {
    rows <- !is.na(out$data[[s]])
    df_s <- out$data[rows, ]
    if (nrow(df_s) < 100) next

    args <- list(
      object   = setNames(list(df_s[[s]]), s),
      data     = df_s,
      formula  = out$formula,
      times    = out$time,
      metrics  = c("brier", "auc"),
      summary  = "ibs"
    )
    if (!is.null(out$cause)) args$cause <- out$cause

    score_obj <- tryCatch(do.call(Score, args), error = function(e) NULL)
    if (is.null(score_obj)) {
      message(sprintf("  skipped %s x %s (Score failed)", out_name, s))
      next
    }

    brier <- score_obj$Brier$score
    cal_rows[[length(cal_rows) + 1]] <- data.frame(
      outcome   = out_name,
      score     = s,
      horizon_y = out$time,
      brier     = brier$Brier[brier$model != "Null model"][1]
    )

    # Distribution plot
    p <- ggplot(df_s, aes(x = .data[[s]])) +
      geom_histogram(bins = 50) +
      labs(title = sprintf("%s distribution (%s cohort)", s, out_name),
           x = s, y = "count") +
      theme_minimal()
    ggsave(sprintf("results/calibration/dist_%s_%s.png", out_name, s),
           p, width = 6, height = 4, dpi = 100)
  }
}

cal_df <- do.call(rbind, cal_rows)
write.csv(cal_df, "results/calibration_summary.csv", row.names = FALSE)
message("\n=== Calibration summary ===")
print(cal_df)
