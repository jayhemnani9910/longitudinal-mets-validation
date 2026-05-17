# ACC/AHA Pooled Cohort Equations (PCE) per Goff et al. 2014 (Circulation),
# 2018 update (Yadlowsky et al. Ann Intern Med).
#
# Implements 10-year ASCVD risk for four race/sex strata:
#   - white female, black female, white male, black male
#
# For Hispanic, Asian, and Other Hispanic NHANES groups: the ACC/AHA app uses
# the white-stratum equations as a default approximation. This is the standard
# practice in published US validation studies but is a known limitation
# documented in the manuscript.

#' Compute 10-year ASCVD risk per Pooled Cohort Equations
#'
#' @param age years (40-79 in original derivation)
#' @param sex "male" or "female"
#' @param race "white" or "black". For other races, "white" is used per ACC/AHA app
#' @param total_chol mg/dl
#' @param hdl mg/dl
#' @param sbp mm Hg
#' @param bp_treated logical
#' @param smoker logical
#' @param diabetes logical
#' @return 10-year ASCVD risk in [0,1]
pce <- function(age, sex, race, total_chol, hdl, sbp,
                bp_treated, smoker, diabetes) {
  # Validate inputs
  if (any(is.na(c(age, total_chol, hdl, sbp)))) return(NA_real_)
  if (!sex %in% c("male", "female")) stop("sex must be 'male' or 'female'")
  if (!race %in% c("white", "black")) race <- "white"  # default approximation

  ln_age <- log(age)
  ln_tc  <- log(total_chol)
  ln_hdl <- log(hdl)
  ln_sbp <- log(sbp)

  # Coefficients per Goff 2014 Table A
  if (sex == "female" && race == "white") {
    s0 <- 0.9665
    mean_terms <- -29.18
    sum_terms <- -29.799 * ln_age + 4.884 * (ln_age ^ 2) +
                  13.540 * ln_tc  + -3.114 * ln_age * ln_tc +
                 -13.578 * ln_hdl +  3.149 * ln_age * ln_hdl +
                   2.019 * ln_sbp * as.numeric(bp_treated) +
                   1.957 * ln_sbp * as.numeric(!bp_treated) +
                   7.574 * as.numeric(smoker) +
                  -1.665 * ln_age * as.numeric(smoker) +
                   0.661 * as.numeric(diabetes)
  } else if (sex == "female" && race == "black") {
    s0 <- 0.9533
    mean_terms <- 86.61
    sum_terms <- 17.114 * ln_age +
                  0.940 * ln_tc +
                -18.920 * ln_hdl + 4.475 * ln_age * ln_hdl +
                 29.291 * ln_sbp * as.numeric(bp_treated) +
                 -6.432 * ln_age * ln_sbp * as.numeric(bp_treated) +
                 27.820 * ln_sbp * as.numeric(!bp_treated) +
                 -6.087 * ln_age * ln_sbp * as.numeric(!bp_treated) +
                  0.691 * as.numeric(smoker) +
                  0.874 * as.numeric(diabetes)
  } else if (sex == "male" && race == "white") {
    s0 <- 0.9144
    mean_terms <- 61.18
    sum_terms <- 12.344 * ln_age +
                 11.853 * ln_tc + -2.664 * ln_age * ln_tc +
                 -7.990 * ln_hdl + 1.769 * ln_age * ln_hdl +
                  1.797 * ln_sbp * as.numeric(bp_treated) +
                  1.764 * ln_sbp * as.numeric(!bp_treated) +
                  7.837 * as.numeric(smoker) +
                 -1.795 * ln_age * as.numeric(smoker) +
                  0.658 * as.numeric(diabetes)
  } else if (sex == "male" && race == "black") {
    s0 <- 0.8954
    mean_terms <- 19.54
    sum_terms <- 2.469 * ln_age +
                  0.302 * ln_tc +
                 -0.307 * ln_hdl +
                  1.916 * ln_sbp * as.numeric(bp_treated) +
                  1.809 * ln_sbp * as.numeric(!bp_treated) +
                  0.549 * as.numeric(smoker) +
                  0.645 * as.numeric(diabetes)
  }

  ind_sum <- sum_terms - mean_terms
  risk <- 1 - s0 ^ exp(ind_sum)
  pmin(pmax(risk, 0), 1)
}

#' Vectorized PCE for a data frame
#'
#' @param df data frame
#' @param ... column names: age_col, sex_col, race_col, tc_col, hdl_col,
#'   sbp_col, bptx_col, smoker_col, dm_col
pce_vec <- function(df,
                    age_col = "age", sex_col = "sex", race_col = "race",
                    tc_col = "total_chol", hdl_col = "hdl",
                    sbp_col = "sbp", bptx_col = "bp_treatment",
                    smoker_col = "current_smoker", dm_col = "prior_t2d") {
  mapply(
    pce,
    df[[age_col]], df[[sex_col]], df[[race_col]],
    df[[tc_col]], df[[hdl_col]], df[[sbp_col]],
    df[[bptx_col]], df[[smoker_col]], df[[dm_col]]
  )
}
