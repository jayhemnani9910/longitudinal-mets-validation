# B9 synthetic feature engineering per Shin/Shim/Oh PLoS One 2023
# (Machine learning-based predictive model for prevention of metabolic syndrome).
#
# These features feed both the B9 decision tree (R/scores/b9_tree.R) and the
# US-refit CART (Task 5.2b).

source("R/scores/rmrs.R")  # for elliot_sigmoid + scale_factor

#' Compute B9 synthetic features from raw measurements.
#'
#' B9 uses 5 derived features:
#'   WC       = elliot-sigmoid-scaled waist circumference (sex/ethnicity-specific)
#'   BP       = elliot-sigmoid-scaled blood pressure (max of sbp/dbp excess)
#'   BPWC_add = BP + WC          (joint elevation)
#'   BPWC_mul = BP * WC          (interaction)
#'   BPWC_dif = BP - WC          (imbalance)
#'
#' @param waist Waist circumference (cm)
#' @param sbp Systolic BP (mm Hg)
#' @param dbp Diastolic BP (mm Hg)
#' @param sex "male" or "female"
#' @param ethnicity "us" or "kr"
#' @return named list with WC, BP, BPWC_add, BPWC_mul, BPWC_dif
b9_features <- function(waist, sbp, dbp, sex, ethnicity = "us") {
  wc_key <- sprintf("waist_%s_%s", sex, ethnicity)
  WC <- scale_factor(waist, wc_key)
  BP <- scale_factor(c(sbp, dbp), "bp")
  list(
    WC       = WC,
    BP       = BP,
    BPWC_add = BP + WC,
    BPWC_mul = BP * WC,
    BPWC_dif = BP - WC
  )
}
