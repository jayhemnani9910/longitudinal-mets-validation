# NHANES complex-survey design constructor.
# Per NHANES analytic guidelines, when combining N cycles divide MEC weights
# by N. Use SDMVPSU + SDMVSTRA for design.

suppressMessages({
  library(survey)
  library(dplyr)
})

#' Build a survey design object for a multi-cycle NHANES cohort.
#'
#' The analysis cohort is the morning fasting subsample (it requires fasting
#' glucose and triglycerides), so the NCHS-correct design weight is WTSAF2YR
#' (column wt_fast), not the MEC exam weight. Subjects with a zero fasting
#' weight were not selected into the fasting subsample and carry no
#' representation, so they are dropped from the weighted design.
#'
#' @param df data frame with columns sdmvpsu, sdmvstra, the chosen weight, cycle.
#' @param weight_col name of the design-weight column (default the fasting weight).
#' @param n_cycles number of 2-year cycles pooled in the dataset (default 10 for
#'   the full 1999-2018 span; pass the pooled count explicitly for a subset).
#' @return svydesign object.
build_survey_design <- function(df, weight_col = "wt_fast", n_cycles = 10L) {
  w <- df[[weight_col]] / n_cycles
  keep <- !is.na(w) & w > 0
  df <- df[keep, ]
  df$weight_combined <- w[keep]
  svydesign(
    ids      = ~sdmvpsu,
    strata   = ~sdmvstra,
    weights  = ~weight_combined,
    nest     = TRUE,
    data     = df
  )
}
