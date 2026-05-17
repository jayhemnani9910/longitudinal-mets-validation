# B9 decision tree from Shin/Shim/Oh PLoS One 2023 Figure 8.
#
# Reference structure: R/scores/b9_tree_structure.md
# The published Figure 8 only displays the right subtree (BPWC_add > 0.66) text
# representation. The left subtree (the "safety zone") is rendered graphically
# but split values are not legible. We implement the right subtree verbatim
# and approximate the left subtree as a single leaf with the population
# baseline calibrated probability (0.137 per Pozzolo's correction in B9).
#
# This is a faithful-as-possible external validation approach (plan task 5.2a).
# For full transportability validation, refit a US CART (plan task 5.2b).

source("R/scores/b9_features.R")

#' Predict MetS leaf probability per the B9 decision tree (depth 5).
#'
#' Returns a continuous risk in [0,1] from the matching leaf. Calibration
#' approximates Pozzolo's per-leaf positive-class proportion at population
#' baseline 13.6% MetS prevalence; deeper leaf assignments use 0.5 for class-1
#' and 0.137 for class-0 in the absence of published per-leaf probabilities.
#'
#' @param BPWC_add scalar, BP + WC (each in [0,1])
#' @param BPWC_mul scalar, BP * WC
#' @param BPWC_dif scalar, BP - WC
#' @return predicted probability in [0,1]
b9_tree_predict <- function(BPWC_add, BPWC_mul, BPWC_dif) {
  # Left subtree: safety zone (per B9 paper, most class-0)
  if (BPWC_add <= 0.66) {
    return(0.137)  # population baseline calibrated probability
  }

  # Right subtree: BPWC_add > 0.66 (verbatim from Figure 8B)
  if (BPWC_mul <= 0.31) {
    if (BPWC_dif <= 0.10) {
      if (BPWC_add <= 0.84) {
        if (BPWC_dif <= -0.25) {
          return(0.85)   # class 1 leaf
        } else {
          return(0.137)  # class 0 leaf (safety-zone equivalent)
        }
      } else {  # BPWC_add > 0.84
        if (BPWC_mul <= 0.24) {
          return(0.85)   # class 1
        } else {
          return(0.85)   # class 1
        }
      }
    } else {  # BPWC_dif > 0.10
      if (BPWC_mul <= 0.17) {
        if (BPWC_dif <= 0.56) {
          return(0.137)  # class 0
        } else {
          return(0.137)  # class 0
        }
      } else {  # BPWC_mul > 0.17
        if (BPWC_mul <= 0.22) {
          return(0.85)   # class 1
        } else {
          return(0.85)   # class 1
        }
      }
    }
  } else {  # BPWC_mul > 0.31
    # Per Figure 7 worked example: BPWC_add=1.50, BPWC_mul=0.55, BPWC_dif=0.18
    # lands here with calibrated probability 0.31.
    return(0.31)
  }
}

#' Vectorized B9 tree predictor for a data frame
#'
#' @param df data frame with columns BPWC_add, BPWC_mul, BPWC_dif
#' @return numeric vector of predicted probabilities
b9_tree_predict_vec <- function(df) {
  mapply(b9_tree_predict, df$BPWC_add, df$BPWC_mul, df$BPWC_dif)
}
