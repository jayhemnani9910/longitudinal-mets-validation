# NHANES complex-survey design constructor.
# Per NHANES analytic guidelines, when combining N cycles divide MEC weights
# by N. Use SDMVPSU + SDMVSTRA for design.

suppressMessages({
  library(survey)
  library(dplyr)
})

#' Build a survey design object for a multi-cycle NHANES cohort.
#'
#' @param df data frame with columns sdmvpsu, sdmvstra, wt_mec, cycle.
#' @return svydesign object.
build_survey_design <- function(df) {
  n_cycles <- length(unique(df$cycle))
  df$weight_combined <- df$wt_mec / n_cycles
  svydesign(
    ids      = ~sdmvpsu,
    strata   = ~sdmvstra,
    weights  = ~weight_combined,
    nest     = TRUE,
    data     = df
  )
}
