#!/usr/bin/env Rscript
# scripts/12_calibration_curves.R
#
# Predicted-versus-observed calibration curves on the recalibrated absolute-risk
# scale, per (outcome, score). Subjects are binned by predicted risk; within each
# bin the mean predicted risk (x) is plotted against the observed cumulative
# incidence at the horizon (y), estimated by Aalen-Johansen for competing-risk
# outcomes and 1 - Kaplan-Meier for all-cause. The 45-degree line is perfect
# calibration. Uses the same recalibration as the decision curve analysis.
#
# Inputs:  data/processed/cohort_with_scores.rds, R/utils/dca_competing.R
# Outputs: results/plots/calibration_<outcome>.png, results/calibration_curve_points.csv

.libPaths("/home/po/projects/work/longitudinal-mets-validation/renv/library/R-4.3/x86_64-pc-linux-gnu")

suppressMessages({
  library(survival)
  library(riskRegression)
  library(prodlim)
  library(ggplot2)
  library(dplyr)
})

source("R/utils/dca_competing.R")

df       <- readRDS("data/processed/cohort_with_scores.rds")
df_cause <- df[df$cause_coded, ]

score_labels <- c(rmrs_score = "RMRS", b9_score = "B9 tree", pce_score = "PCE",
                  framingham_score = "Framingham", findrisc_score = "FINDRISC")

outcomes <- list(
  allcause = list(data = df, time = "followup_years", status = "event_allcause",
                  cause = 1, competing = FALSE, horizon = 9.5, bins = 10,
                  scores = c("rmrs_score", "b9_score", "pce_score",
                             "framingham_score", "findrisc_score"),
                  title = "All-cause mortality, 9.5-year horizon"),
  cv       = list(data = df_cause, time = "followup_years", status = "competing_cv",
                  cause = 1, competing = TRUE, horizon = 9.5, bins = 5,
                  scores = c("rmrs_score", "b9_score", "pce_score", "framingham_score"),
                  title = "Cardiovascular mortality, 9.5-year horizon"),
  dm       = list(data = df_cause, time = "followup_years_dm", status = "competing_dm",
                  cause = 1, competing = TRUE, horizon = 14.5, bins = 5,
                  scores = c("rmrs_score", "b9_score", "findrisc_score"),
                  title = "Diabetes-related mortality, 14.5-year horizon")
)

dir.create("results/plots", recursive = TRUE, showWarnings = FALSE)
all_points <- list()

for (out_name in names(outcomes)) {
  out <- outcomes[[out_name]]
  pts <- list()

  for (s in out$scores) {
    risk <- recalibrated_risk(out$data, s, out$time, out$status,
                              out$cause, out$horizon, out$competing)
    ok <- !is.na(risk)
    if (sum(ok) < 100) next
    d <- data.frame(risk = risk[ok],
                    time = out$data[[out$time]][ok],
                    status = out$data[[out$status]][ok])

    # Quantile bins on predicted risk; collapse if the score is near-degenerate.
    g <- out$bins
    brks <- unique(quantile(d$risk, probs = seq(0, 1, length.out = g + 1),
                            na.rm = TRUE))
    if (length(brks) < 3) next
    d$bin <- cut(d$risk, breaks = brks, include.lowest = TRUE, labels = FALSE)

    for (b in sort(unique(d$bin))) {
      sub <- d[d$bin == b, ]
      if (nrow(sub) < 20) next
      obs <- .cif_at(sub$time, sub$status, out$cause, out$horizon)
      pts[[length(pts) + 1]] <- data.frame(
        outcome   = out_name,
        score     = score_labels[s],
        bin       = b,
        n         = nrow(sub),
        predicted = mean(sub$risk),
        observed  = obs
      )
    }
  }

  if (length(pts) == 0) next
  pdf_pts <- do.call(rbind, pts)
  all_points[[out_name]] <- pdf_pts

  lim <- max(pdf_pts$predicted, pdf_pts$observed, na.rm = TRUE) * 1.05
  p <- ggplot(pdf_pts, aes(x = predicted, y = observed)) +
    geom_abline(slope = 1, intercept = 0, linetype = "dashed", color = "grey50") +
    geom_line(color = "#1B4F72") +
    geom_point(aes(size = n), color = "#1B4F72", alpha = 0.8) +
    facet_wrap(~score, scales = "free") +
    scale_size_continuous(range = c(1.5, 4), guide = "none") +
    labs(title = out$title,
         x = "Predicted risk (recalibrated)",
         y = "Observed cumulative incidence") +
    theme_minimal(base_size = 11) +
    theme(plot.title = element_text(face = "bold"))
  ggsave(sprintf("results/plots/calibration_%s.png", out_name),
         p, width = 8, height = 5.5, dpi = 110)
  message(sprintf("Wrote results/plots/calibration_%s.png (%d points)",
                  out_name, nrow(pdf_pts)))
}

cp <- do.call(rbind, all_points)
write.csv(cp, "results/calibration_curve_points.csv", row.names = FALSE)
message("\n=== Calibration curve points ===")
print(cp, row.names = FALSE)
