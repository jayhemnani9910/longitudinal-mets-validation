# Survey-weighted Fine-Gray subdistribution model, shared by scripts 06 and 07.
#
# riskRegression::FGR / cmprsk::crr cannot carry NHANES survey weights, so the
# subdistribution coefficient is fit the Geskus way: survival::finegray() builds
# the weighted risk set (its fgwt already folds in the prior case weights passed
# via `weights`), then a Cox model on that expanded data with weights = fgwt
# recovers the Fine-Gray coefficient. Robust standard errors clustered on the
# PSU give design-consistent inference.

suppressMessages({
  library(survival)
})

#' Survey-weighted Fine-Gray subdistribution hazard ratio for one score.
#'
#' @param df data frame holding the columns below
#' @param score column name of the marker
#' @param time_col follow-up time column
#' @param status_col competing-risk status (0 censored, 1 cause, 2 competing)
#' @param weight_col survey weight column (e.g. wt_mec)
#' @param cluster_col PSU cluster id column for robust SE
#' @param cause integer cause of interest (default 1)
#' @return list(hr, lo, hi, se, n): per-1-SD subdistribution HR, or NAs on failure
weighted_finegray_hr <- function(df, score, time_col, status_col,
                                  weight_col, cluster_col, cause = 1) {
  ok <- !is.na(df[[score]]) & !is.na(df[[time_col]]) &
        !is.na(df[[status_col]]) & !is.na(df[[weight_col]])
  d <- df[ok, ]
  if (nrow(d) < 50) return(list(hr = NA, lo = NA, hi = NA, se = NA, n = nrow(d)))

  status_f <- factor(d[[status_col]], levels = c(0, 1, 2),
                     labels = c("censor", "cause", "competing"))
  d$.status_f <- status_f
  d$.time     <- d[[time_col]]
  # Standardize the marker so the hazard ratio is per 1 SD, interpretable and
  # comparable across scores on different native scales (0-1 vs 0-26).
  sd_s <- stats::sd(d[[score]], na.rm = TRUE)
  d$.score    <- (d[[score]] - mean(d[[score]], na.rm = TRUE)) /
                 ifelse(sd_s > 0, sd_s, 1)
  d$.w        <- d[[weight_col]]
  d$.cluster  <- d[[cluster_col]]
  etype <- if (cause == 1) "cause" else "competing"

  fg <- tryCatch(
    survival::finegray(survival::Surv(.time, .status_f) ~ .score + .cluster,
                       data = d, weights = .w, etype = etype),
    error = function(e) NULL
  )
  if (is.null(fg)) return(list(hr = NA, lo = NA, hi = NA, se = NA, n = nrow(d)))

  fit <- tryCatch(
    survival::coxph(
      survival::Surv(fgstart, fgstop, fgstatus) ~ .score,
      data = fg, weights = fgwt, cluster = .cluster, robust = TRUE
    ),
    error = function(e) NULL
  )
  if (is.null(fit)) return(list(hr = NA, lo = NA, hi = NA, se = NA, n = nrow(d)))

  beta <- as.numeric(coef(fit)[1])
  se   <- sqrt(diag(vcov(fit)))[1]
  list(hr = exp(beta),
       lo = exp(beta - 1.96 * se),
       hi = exp(beta + 1.96 * se),
       se = as.numeric(se),
       n  = nrow(d))
}
