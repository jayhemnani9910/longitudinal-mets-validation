# Finnish Diabetes Risk Score (FINDRISC) per Lindstrom & Tuomilehto 2003,
# Diabetes Care 26(3):725-731.
#
# Additive ordinal score in [0, 26]. Higher = higher 10-year T2D risk.

#' FINDRISC score (0-26 points)
#'
#' @param age years
#' @param bmi kg/m^2
#' @param waist cm
#' @param sex "male" or "female"
#' @param physical_active logical (TRUE if >=30 min daily physical activity)
#' @param daily_vegetables logical (TRUE if eats vegetables/fruits daily)
#' @param bp_treatment logical (currently on antihypertensive medication)
#' @param prior_high_glucose logical (ever told blood glucose was high)
#' @param family_history "none", "second_degree", or "first_degree"
#' @return integer score in [0, 26]
findrisc <- function(age, bmi, waist, sex,
                     physical_active, daily_vegetables,
                     bp_treatment, prior_high_glucose,
                     family_history) {
  if (any(is.na(c(age, bmi, waist)))) return(NA_integer_)
  if (!sex %in% c("male", "female")) stop("sex must be 'male' or 'female'")
  if (!family_history %in% c("none", "second_degree", "first_degree")) {
    family_history <- "none"
  }

  pts <- 0L

  # Age points (0-4)
  pts <- pts + if (age < 45) 0L
               else if (age < 55) 2L
               else if (age < 65) 3L
               else 4L

  # BMI points (0-3)
  pts <- pts + if (bmi < 25) 0L
               else if (bmi < 30) 1L
               else 3L

  # Waist points (0-4, sex-specific)
  pts <- pts + if (sex == "male") {
    if (waist < 94) 0L
    else if (waist < 102) 3L
    else 4L
  } else {
    if (waist < 80) 0L
    else if (waist < 88) 3L
    else 4L
  }

  # Physical activity (TRUE = 0 pts, FALSE = 2 pts)
  pts <- pts + if (physical_active) 0L else 2L

  # Daily vegetables/fruits (TRUE = 0 pts, FALSE = 1 pt)
  pts <- pts + if (daily_vegetables) 0L else 1L

  # BP medication (TRUE = 2 pts, FALSE = 0)
  pts <- pts + if (bp_treatment) 2L else 0L

  # Prior high glucose (TRUE = 5 pts)
  pts <- pts + if (prior_high_glucose) 5L else 0L

  # Family history (0/3/5 pts)
  pts <- pts + switch(family_history,
                      "none" = 0L,
                      "second_degree" = 3L,
                      "first_degree" = 5L)

  as.integer(pts)
}

#' Vectorized FINDRISC for a data frame
findrisc_vec <- function(df) {
  mapply(
    findrisc,
    df$age, df$bmi, df$waist_cm, df$sex,
    df$physical_active, df$daily_vegetables,
    df$bp_treatment, df$prior_high_glucose,
    df$family_history_dm
  )
}
