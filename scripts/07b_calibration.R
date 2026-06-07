#!/usr/bin/env Rscript
# scripts/07b_calibration.R
#
# Calibration plots and Brier scores per (outcome, score). Uses
# riskRegression::Score with IPCW correction.
#
# Brier is computed on each score's RECALIBRATED horizon-specific absolute risk
# (Fine-Gray CIF for the competing-risks outcomes, Cox complement for all-cause),
# the same mapping used by the decision curve analysis (R/utils/dca_competing.R).
# Brier on the raw score is meaningless because several scores are not on a
# probability scale (FINDRISC is an integer points count; RMRS sits in a narrow
# band), which previously produced impossible Brier values. The recalibration is
# fit and evaluated on the same cohort, so the Brier is apparent in-sample,
# matching the calibration curves and the DCA.
#
# Inputs:  data/processed/cohort_with_scores.rds
# Outputs: results/calibration/dist_<outcome>_<score>.png, results/calibration_summary.csv

.libPaths("renv/library/R-4.3/x86_64-pc-linux-gnu")

suppressMessages({
  library(survival)
  library(riskRegression)
  library(prodlim)
  library(dplyr)
  library(ggplot2)
})

source("R/utils/dca_competing.R")   # recalibrated_risk()

df <- readRDS("data/processed/cohort_with_scores.rds")

dir.create("results/calibration", recursive = TRUE, showWarnings = FALSE)

# Cause-specific outcomes use the 1999-2014 cause-coded cohort and the same
# off-cap horizons as the survival scripts (9.5y for all-cause/CV, 14.5y for
# the 15y diabetes frame). Diabetes uses the followup_years_dm time variable.
df_cause <- df[df$cause_coded, ]

outcomes <- list(
  allcause = list(formula = Surv(followup_years, event_allcause) ~ 1,
                  data = df, time = 9.5,
                  time_col = "followup_years", status_col = "event_allcause",
                  cause = 1, competing = FALSE,
                  scores  = c("rmrs_score", "b9_score", "pce_score",
                              "framingham_score", "findrisc_score")),
  cv       = list(formula = Hist(followup_years, competing_cv) ~ 1,
                  cause = 1, data = df_cause, time = 9.5,
                  time_col = "followup_years", status_col = "competing_cv",
                  competing = TRUE,
                  scores = c("rmrs_score", "b9_score", "pce_score",
                             "framingham_score")),
  dm       = list(formula = Hist(followup_years_dm, competing_dm) ~ 1,
                  cause = 1, data = df_cause, time = 14.5,
                  time_col = "followup_years_dm", status_col = "competing_dm",
                  competing = TRUE,
                  scores = c("rmrs_score", "b9_score", "findrisc_score"))
)

cal_rows <- list()
for (out_name in names(outcomes)) {
  out <- outcomes[[out_name]]
  for (s in out$scores) {
    rows <- !is.na(out$data[[s]])
    df_s <- out$data[rows, ]
    if (nrow(df_s) < 100) next

    # Map the raw score to a horizon-specific absolute risk before scoring.
    risk <- recalibrated_risk(df_s, s, out$time_col, out$status_col,
                              out$cause, out$time, out$competing)
    keep <- !is.na(risk)
    if (sum(keep) < 100) {
      message(sprintf("  skipped %s x %s (recalibration produced no risk)",
                      out_name, s))
      next
    }
    df_k   <- df_s[keep, ]
    risk_k <- risk[keep]

    args <- list(
      object   = setNames(list(risk_k), s),
      data     = df_k,
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
