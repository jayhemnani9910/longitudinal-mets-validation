#!/usr/bin/env Rscript
# scripts/08_dca.R
#
# Decision Curve Analysis at 10 years per outcome over threshold 0.01-0.30.
# Net benefit relative to "treat all" / "treat none" reference strategies.
#
# Inputs:  data/processed/cohort_with_scores.rds
# Outputs: results/dca_<outcome>.{rds,png}

.libPaths("/home/po/projects/work/longitudinal-mets-validation/renv/library/R-4.3/x86_64-pc-linux-gnu")

suppressMessages({
  library(dcurves)
  library(survival)
  library(dplyr)
  library(ggplot2)
})

df <- readRDS("data/processed/cohort_with_scores.rds")

outcomes <- list(
  allcause = list(event_col = "event_allcause",
                  scores    = c("rmrs_score", "b9_score", "pce_score",
                                "framingham_score", "findrisc_score")),
  cv       = list(event_col = "event_cv",
                  scores    = c("rmrs_score", "b9_score", "pce_score",
                                "framingham_score")),
  dm       = list(event_col = "event_dm",
                  scores    = c("rmrs_score", "b9_score", "findrisc_score"))
)

for (out_name in names(outcomes)) {
  out <- outcomes[[out_name]]
  message(sprintf("\n--- DCA for %s ---", out_name))

  # Subset to non-NA in any included score
  rows <- complete.cases(df[, out$scores, drop = FALSE])
  df_s <- df[rows, ]

  formula <- as.formula(sprintf(
    "Surv(followup_years, %s) ~ %s",
    out$event_col,
    paste(out$scores, collapse = " + ")
  ))

  dca_result <- tryCatch(
    dcurves::dca(
      formula,
      data       = df_s,
      time       = 10,
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

  message(sprintf("  Saved results/dca_%s.{rds,png}", out_name))
}
