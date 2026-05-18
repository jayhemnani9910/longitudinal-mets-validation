#!/usr/bin/env Rscript
# scripts/08_dca.R
#
# Decision Curve Analysis per outcome over threshold 0.01-0.30. Net benefit
# relative to "treat all" / "treat none" reference strategies.
#
# Horizons follow the Phase 2 survival scripts: 9.5y for all-cause and CV
# (10y shifted off the IPCW boundary), 14.5y for diabetes (15y horizon
# shifted off the boundary). The diabetes outcome uses the broadened
# definition built in script 03 (event_dm + followup_years_dm).
#
# Note on calibration: dcurves::dca needs each marker to be a predicted
# probability in [0, 1]. RMRS / B9 / PCE / Framingham already produce
# probability-style outputs (range checked at load), but FINDRISC is a
# raw integer count (0-26). We convert each score to a horizon-calibrated
# predicted risk via a per-outcome Cox model with the score as the only
# covariate, then evaluate at the horizon. This recalibrates every score
# to NHANES so the DCA reads as net benefit of using that score's ranking,
# not of its absolute risk scale.
#
# Inputs:  data/processed/cohort_with_scores.rds
# Outputs: results/dca_<outcome>.{rds,png}, results/plots/dca_<outcome>.png

.libPaths("/home/po/projects/work/longitudinal-mets-validation/renv/library/R-4.3/x86_64-pc-linux-gnu")

suppressMessages({
  library(dcurves)
  library(survival)
  library(dplyr)
  library(ggplot2)
})

df <- readRDS("data/processed/cohort_with_scores.rds")

outcomes <- list(
  allcause = list(time_col  = "followup_years",
                  event_col = "event_allcause",
                  horizon   = 9.5,
                  scores    = c("rmrs_score", "b9_score", "pce_score",
                                "framingham_score", "findrisc_score")),
  cv       = list(time_col  = "followup_years",
                  event_col = "event_cv",
                  horizon   = 9.5,
                  scores    = c("rmrs_score", "b9_score", "pce_score",
                                "framingham_score")),
  dm       = list(time_col  = "followup_years_dm",
                  event_col = "event_dm",
                  horizon   = 14.5,
                  scores    = c("rmrs_score", "b9_score", "findrisc_score"))
)

dir.create("results/plots", recursive = TRUE, showWarnings = FALSE)

# Convert a raw score column into a horizon-calibrated predicted risk via a
# Cox model. Returns a numeric vector of predicted P(event by horizon | score)
# values, NA for rows outside `rows`.
cox_risk_at_horizon <- function(df, score, time_col, event_col, horizon, rows) {
  df_s <- df[rows, ]
  fit <- coxph(as.formula(sprintf("Surv(%s, %s) ~ %s",
                                  time_col, event_col, score)),
               data = df_s)
  sf <- survfit(fit, newdata = df_s)
  # survfit returns a matrix of survival probabilities indexed by event time;
  # we extract the row at the horizon.
  idx <- findInterval(horizon, sf$time)
  if (idx < 1) idx <- 1
  surv_at_h <- sf$surv[idx, ]
  risk <- 1 - surv_at_h
  out <- rep(NA_real_, nrow(df))
  out[rows] <- pmin(pmax(risk, 0), 1)
  out
}

for (out_name in names(outcomes)) {
  out <- outcomes[[out_name]]
  message(sprintf("\n--- DCA for %s (horizon %.1fy) ---",
                  out_name, out$horizon))

  rows <- complete.cases(df[, c(out$scores, out$time_col, out$event_col),
                            drop = FALSE])
  df_s <- df[rows, ]

  # Build per-score calibrated risk columns
  risk_cols <- c()
  for (s in out$scores) {
    rcol <- paste0(s, "_p_", out_name)
    df_s[[rcol]] <- cox_risk_at_horizon(
      df, s, out$time_col, out$event_col, out$horizon, rows
    )[rows]
    risk_cols <- c(risk_cols, rcol)
  }

  formula <- as.formula(sprintf(
    "Surv(%s, %s) ~ %s",
    out$time_col, out$event_col,
    paste(risk_cols, collapse = " + ")
  ))

  dca_result <- tryCatch(
    dcurves::dca(
      formula,
      data       = df_s,
      time       = out$horizon,
      thresholds = seq(0.01, 0.30, by = 0.01)
    ),
    error = function(e) {
      message(sprintf("  DCA failed: %s", e$message))
      NULL
    }
  )

  if (is.null(dca_result)) next

  saveRDS(dca_result, sprintf("results/dca_%s.rds", out_name))

  p <- tryCatch(plot(dca_result, smooth = TRUE),
                error = function(e) plot(dca_result))
  ggsave(sprintf("results/dca_%s.png", out_name),
         p, width = 8, height = 6, dpi = 120)
  ggsave(sprintf("results/plots/dca_%s.png", out_name),
         p, width = 8, height = 6, dpi = 120)

  # Persist a tidy net-benefit table at clinically relevant thresholds for
  # the manuscript. dcurves stores a data frame in $dca after extraction via
  # as_tibble; we just snapshot the long-format internal data.
  nb_tbl <- tryCatch(as.data.frame(dca_result$dca), error = function(e) NULL)
  if (!is.null(nb_tbl)) {
    write.csv(nb_tbl,
              sprintf("results/dca_%s_netbenefit.csv", out_name),
              row.names = FALSE)
  }

  message(sprintf("  Saved results/dca_%s.{rds,png,csv} and results/plots/dca_%s.png",
                  out_name, out_name))
}
