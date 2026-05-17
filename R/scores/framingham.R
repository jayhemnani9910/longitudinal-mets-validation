# Framingham 2008 general CVD risk equation per D'Agostino et al.
# Circulation 117(6):743-753 (2008).
#
# Kept as a secondary CVD baseline alongside PCE because B-cluster's reference
# papers often cite Framingham, and including both strengthens comparison.

#' Compute 10-year general CVD risk per Framingham 2008
#'
#' Inputs are the same as PCE except race is not used.
#'
#' @param age years
#' @param sex "male" or "female"
#' @param total_chol mg/dl
#' @param hdl mg/dl
#' @param sbp mm Hg
#' @param bp_treated logical
#' @param smoker logical
#' @param diabetes logical
#' @return 10-year CVD risk in [0,1]
framingham_2008 <- function(age, sex, total_chol, hdl, sbp,
                            bp_treated, smoker, diabetes) {
  if (any(is.na(c(age, total_chol, hdl, sbp)))) return(NA_real_)
  if (!sex %in% c("male", "female")) stop("sex must be 'male' or 'female'")

  ln_age <- log(age)
  ln_tc  <- log(total_chol)
  ln_hdl <- log(hdl)
  ln_sbp <- log(sbp)

  if (sex == "female") {
    # D'Agostino 2008 Table 2 female coefficients
    sum_terms <- 2.32888 * ln_age +
                 1.20904 * ln_tc +
                -0.70833 * ln_hdl +
                 2.76157 * ln_sbp * as.numeric(!bp_treated) +
                 2.82263 * ln_sbp * as.numeric(bp_treated) +
                 0.52873 * as.numeric(smoker) +
                 0.69154 * as.numeric(diabetes)
    s0 <- 0.95012
    mean_terms <- 26.1931
  } else {  # male
    sum_terms <- 3.06117 * ln_age +
                 1.12370 * ln_tc +
                -0.93263 * ln_hdl +
                 1.93303 * ln_sbp * as.numeric(!bp_treated) +
                 1.99881 * ln_sbp * as.numeric(bp_treated) +
                 0.65451 * as.numeric(smoker) +
                 0.57367 * as.numeric(diabetes)
    s0 <- 0.88936
    mean_terms <- 23.9802
  }

  ind_sum <- sum_terms - mean_terms
  risk <- 1 - s0 ^ exp(ind_sum)
  pmin(pmax(risk, 0), 1)
}

#' Vectorized Framingham 2008 for a data frame
framingham_vec <- function(df,
                            age_col = "age", sex_col = "sex",
                            tc_col = "total_chol", hdl_col = "hdl",
                            sbp_col = "sbp", bptx_col = "bp_treatment",
                            smoker_col = "current_smoker",
                            dm_col = "prior_t2d") {
  mapply(
    framingham_2008,
    df[[age_col]], df[[sex_col]],
    df[[tc_col]], df[[hdl_col]], df[[sbp_col]],
    df[[bptx_col]], df[[smoker_col]], df[[dm_col]]
  )
}
