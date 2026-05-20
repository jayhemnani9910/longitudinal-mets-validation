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
  library(riskRegression)
  library(prodlim)
  library(survival)
  library(dplyr)
  library(ggplot2)
})

source("R/utils/dca_competing.R")

df <- readRDS("data/processed/cohort_with_scores.rds")
# Cause-specific outcomes use the 1999-2014 cause-coded cohort.
df_cause <- df[df$cause_coded, ]

# competing = TRUE outcomes recalibrate via the Fine-Gray CIF and score net
# benefit with the competing-risks (Aalen-Johansen) incidence; all-cause uses a
# Cox model and Kaplan-Meier complement. status_col is the competing-risk status
# for competing outcomes, or the 0/1 event indicator for all-cause.
outcomes <- list(
  allcause = list(time_col   = "followup_years",
                  status_col = "event_allcause",
                  competing  = FALSE,
                  horizon    = 9.5, data = df,
                  scores     = c("rmrs_score", "b9_score", "pce_score",
                                 "framingham_score", "findrisc_score")),
  cv       = list(time_col   = "followup_years",
                  status_col = "competing_cv",
                  competing  = TRUE,
                  horizon    = 9.5, data = df_cause,
                  scores     = c("rmrs_score", "b9_score", "pce_score",
                                 "framingham_score")),
  dm       = list(time_col   = "followup_years_dm",
                  status_col = "competing_dm",
                  competing  = TRUE,
                  horizon    = 14.5, data = df_cause,
                  scores     = c("rmrs_score", "b9_score", "findrisc_score"))
)

dir.create("results/plots", recursive = TRUE, showWarnings = FALSE)
thresholds <- seq(0.01, 0.30, by = 0.01)

for (out_name in names(outcomes)) {
  out <- outcomes[[out_name]]
  message(sprintf("\n--- DCA for %s (horizon %.1fy) ---",
                  out_name, out$horizon))

  dat   <- out$data
  tvec  <- dat[[out$time_col]]
  svec  <- dat[[out$status_col]]

  # Per-score calibrated absolute risk at the horizon.
  risks <- lapply(out$scores, function(s)
    recalibrated_risk(dat, s, out$time_col, out$status_col, 1,
                      out$horizon, out$competing))
  names(risks) <- out$scores

  # Net-benefit curve: each score, plus treat-all and treat-none references.
  curve_rows <- list()
  for (th in thresholds) {
    for (s in out$scores) {
      curve_rows[[length(curve_rows) + 1]] <- data.frame(
        outcome = out_name, threshold = th, strategy = s,
        net_benefit = cr_net_benefit(risks[[s]], tvec, svec, 1, out$horizon, th)
      )
    }
    curve_rows[[length(curve_rows) + 1]] <- data.frame(
      outcome = out_name, threshold = th, strategy = "treat_all",
      net_benefit = cr_net_benefit_all(tvec, svec, 1, out$horizon, th)
    )
    curve_rows[[length(curve_rows) + 1]] <- data.frame(
      outcome = out_name, threshold = th, strategy = "treat_none",
      net_benefit = 0
    )
  }
  nb_tbl <- do.call(rbind, curve_rows)

  saveRDS(nb_tbl, sprintf("results/dca_%s.rds", out_name))
  write.csv(nb_tbl, sprintf("results/dca_%s_netbenefit.csv", out_name),
            row.names = FALSE)

  p <- ggplot(nb_tbl, aes(x = threshold, y = net_benefit,
                          color = strategy, linetype = strategy)) +
    geom_line(linewidth = 0.7) +
    coord_cartesian(ylim = c(min(-0.01, min(nb_tbl$net_benefit, na.rm = TRUE)),
                             max(nb_tbl$net_benefit, na.rm = TRUE))) +
    labs(title = sprintf("Decision curve: %s mortality (%.1fy)",
                         out_name, out$horizon),
         x = "Threshold probability", y = "Net benefit") +
    theme_minimal()
  ggsave(sprintf("results/dca_%s.png", out_name), p,
         width = 8, height = 6, dpi = 120)
  ggsave(sprintf("results/plots/dca_%s.png", out_name), p,
         width = 8, height = 6, dpi = 120)

  message(sprintf("  Saved results/dca_%s.{rds,png,csv} and results/plots/dca_%s.png",
                  out_name, out_name))
}
