# Robust Metabolic Syndrome Risk Score (RMRS) per Shin, Shim, Oh 2024 PeerJ CS.
# Implements the Triangular Areal Similarity (TAS) method on 5 MetS factors.
#
# Reference: Shin H, Shim S, Oh S (2024). Robust metabolic syndrome risk score
# based on triangular areal similarity. PeerJ Computer Science 10:e2015.
# DOI: 10.7717/peerj-cs.2015

#' Elliot sigmoid transformation
#'
#' f(x) = 0.5 + 0.5 * x / (1 + |x|)
#' Maps R -> (0,1) with f(0)=0.5, f(-1)=0.25, f(1)=0.75. Used to scale
#' MetS factor measurements relative to their diagnostic thresholds.
#'
#' @param x numeric (scalar or vector)
#' @return same shape as x, values in (0,1)
elliot_sigmoid <- function(x) {
  0.5 + 0.5 * x / (1 + abs(x))
}

#' Scale a MetS risk factor to [0,1] centered at diagnostic threshold
#'
#' Implements equation 1 from Shin et al. 2024 PeerJ CS. Each factor's raw
#' measurement is normalized so the threshold maps to 0 (z-score), then
#' transformed via the Elliot sigmoid so the threshold maps to 0.5.
#'
#' @param x The measurement. For BP, length-2 numeric c(sbp, dbp).
#' @param factor One of:
#'   "glucose", "triglycerides",
#'   "waist_male_us", "waist_female_us",
#'   "waist_male_kr", "waist_female_kr",
#'   "hdl_male", "hdl_female",
#'   "bp"
#' @return Scaled value in (0,1) where 0.5 = at diagnostic threshold.
scale_factor <- function(x, factor) {
  thresholds <- list(
    glucose         = 100,
    triglycerides   = 150,
    waist_male_us   = 102,
    waist_female_us = 88,
    waist_male_kr   = 90,
    waist_female_kr = 85,
    hdl_male        = 40,
    hdl_female      = 50,
    bp_sbp          = 130,
    bp_dbp          = 85
  )
  if (factor == "bp") {
    sbp <- x[1]; dbp <- x[2]
    d_sbp <- thresholds$bp_sbp
    d_dbp <- thresholds$bp_dbp
    # Per B8 eq 1: scale by the difference d_bp = d_sbp - d_dbp = 45 (not by 0.1*d).
    # The systolic and diastolic excesses are compared and the larger taken.
    d_bp <- d_sbp - d_dbp  # 45
    z <- max(sbp - d_sbp, dbp - d_dbp) / (0.1 * d_bp)
    return(elliot_sigmoid(z))
  }
  if (factor %in% c("hdl_male", "hdl_female")) {
    d <- thresholds[[factor]]
    # HDL: lower = higher risk, so flip sign on (x - d)
    z <- -(x - d) / (0.1 * d)
    return(elliot_sigmoid(z))
  }
  d <- thresholds[[factor]]
  z <- (x - d) / (0.1 * d)
  elliot_sigmoid(z)
}

# ---- TAS components ---------------------------------------------------------

#' Area of a triangle on a three-axis radar chart (axes 120 degrees apart)
#'
#' For axis values (a, b, c) the triangle has area
#'   A = (sqrt(3) / 4) * (a*b + b*c + c*a)
#'
#' @param x length-3 numeric, axis values in [0,1]
.triangle_area <- function(x) {
  (sqrt(3) / 4) * (x[1] * x[2] + x[2] * x[3] + x[3] * x[1])
}

#' TAS score for one three-axis radar chart
#'
#' Implements the closeness + severity components per Shin et al. 2024 sec 2.2.
#' Closeness measures how close the patient's triangle is to the diagnostic
#' threshold triangle. Severity measures how far the patient exceeds it.
#'
#' @param x length-3 vector, each in [0,1] (scaled risk-factor values)
#' @return TAS score in [0,1]
tas_score <- function(x) {
  stopifnot(length(x) == 3, all(x >= 0 & x <= 1))
  internal <- pmin(x, 0.5)
  external <- x
  critical <- rep(0.5, 3)
  maximum  <- rep(1.0, 3)

  area_internal <- .triangle_area(internal)
  area_external <- .triangle_area(external)
  area_critical <- .triangle_area(critical)
  area_maximum  <- .triangle_area(maximum)

  closeness <- area_internal / area_critical
  severity  <- max(0, (area_external - area_critical) /
                     (area_maximum  - area_critical))

  (closeness + severity) / 2
}

# ---- RMRS aggregation -------------------------------------------------------

#' Compute RMRS from 5 pre-scaled MetS factors
#'
#' Generates all C(5,3) = 10 three-axis combinations, computes TAS for each,
#' returns sqrt of the mean.
#'
#' @param scaled length-5 numeric in [0,1]. Order is symmetric; any consistent
#'               ordering of (WC, BP, HDL, glucose, TG) gives the same result.
rmrs_from_scaled <- function(scaled) {
  stopifnot(length(scaled) == 5)
  combos <- combn(5, 3, simplify = FALSE)
  tas_scores <- vapply(combos, function(idx) tas_score(scaled[idx]), numeric(1))
  sqrt(mean(tas_scores))
}

#' Compute RMRS from raw measurements
#'
#' @param waist Waist circumference (cm)
#' @param sbp Systolic BP (mm Hg)
#' @param dbp Diastolic BP (mm Hg)
#' @param hdl HDL cholesterol (mg/dl)
#' @param glucose Fasting glucose (mg/dl)
#' @param triglycerides (mg/dl)
#' @param sex "male" or "female"
#' @param ethnicity "us" or "kr" for WC threshold selection
#' @return RMRS in (0, 1)
rmrs <- function(waist, sbp, dbp, hdl, glucose, triglycerides,
                 sex, ethnicity = "us") {
  wc_key  <- sprintf("waist_%s_%s", sex, ethnicity)
  hdl_key <- sprintf("hdl_%s", sex)
  scaled <- c(
    scale_factor(waist,         wc_key),
    scale_factor(c(sbp, dbp),   "bp"),
    scale_factor(hdl,           hdl_key),
    scale_factor(glucose,       "glucose"),
    scale_factor(triglycerides, "triglycerides")
  )
  rmrs_from_scaled(scaled)
}
