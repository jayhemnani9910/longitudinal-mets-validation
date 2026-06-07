# B9 decision tree from Shin H, Oh S, Shim S (2023), "Machine learning-based
# predictive model for prevention of metabolic syndrome", PLoS ONE 18(6):e0286635.
#
# Full 22-leaf structure (depth 5) transcribed verbatim from the paper's
# Supporting Information S1 Table (the complete final ruleset). The main-text
# Figure 8 is only partly legible; S1 Table gives the authoritative split
# variables, thresholds, and the binary class per leaf (0 = non-MetS,
# 1 = MetS). It does NOT tabulate per-leaf probabilities.
#
# Per-leaf risk. The paper expresses leaf risk as a "risk multiple" relative to
# the 13.7% MetS prevalence (calibrated probability = risk_multiple * 0.137).
# These multiples appear only on the published risk map (Figure 12A). They are
# read from that figure here, anchored to the paper's own worked example
# (Figure 7: a subject routing to Leaf 21 with stated risk 0.31 = 2.25x). Because
# the paper did not publish a numeric per-leaf probability table, these leaf
# risks are approximate to the resolution of Figure 12A; the class labels (from
# S1 Table) are exact. This limitation is reported in the manuscript, and a CART
# refit on US data (scripts/15_b9_refit.R) is the complementary fair-tree check.
#
# Inputs are the three engineered features from b9_features.R, computed with the
# KOREAN waist thresholds the tree was trained on (ethnicity = "kr"):
#   BPWC_add = BP + WC,  BPWC_mul = BP * WC,  BPWC_dif = BP - WC
# where BP and WC are Elliot-sigmoid-scaled blood pressure and waist in [0,1].
#
# NOTE: depends on elliot_sigmoid() + scale_factor() from rmrs.R and the feature
# builder in b9_features.R. Source rmrs.R then b9_features.R before this file.

B9_PREVALENCE <- 0.137  # MetS prevalence the paper calibrates leaf risk against

# Per-leaf risk multiples (relative to B9_PREVALENCE), keyed by the leaf numbers
# in S1 Table. class 0 leaves are < 1.0 (risk below the prevalence threshold);
# class 1 leaves are >= 1.0. Leaf 21 is fixed to the Figure 7 worked example
# (2.25); the others are read from the Figure 12A risk map.
B9_LEAF_RISK <- c(
  L1  = 0.7,  L2  = 0.1,  L3  = 0.05, L4  = 0.2,  L5  = 0.1,
  L6  = 0.3,  L7  = 0.3,  L8  = 0.1,  L9  = 0.2,  L10 = 1.8,
  L11 = 0.5,  L12 = 2.3,  L13 = 4.0,  L14 = 0.9,  L15 = 0.5,
  L16 = 1.7,  L17 = 2.5,  L18 = 5.0,  L19 = 1.3,  L20 = 5.7,
  L21 = 2.25, L22 = 2.2
)

#' Assign a subject to the matching B9 leaf (1..22) per the S1 Table ruleset.
#'
#' @param BPWC_add scalar, BP + WC
#' @param BPWC_mul scalar, BP * WC
#' @param BPWC_dif scalar, BP - WC
#' @return integer leaf index 1..22
b9_tree_leaf <- function(BPWC_add, BPWC_mul, BPWC_dif) {
  if (BPWC_add <= 0.66) {                       # left subtree ("safety zone")
    if (BPWC_dif <= -0.34) return(1L)
    if (BPWC_add <= 0.42) {
      if (BPWC_add <= 0.33) {
        if (BPWC_dif <= -0.13) return(2L) else return(3L)
      } else {
        if (BPWC_dif <= -0.10) return(4L) else return(5L)
      }
    } else {                                    # BPWC_add > 0.42
      if (BPWC_dif <= -0.06) {
        if (BPWC_mul <= 0.06) return(6L) else return(7L)
      } else {
        if (BPWC_add <= 0.56) return(8L) else return(9L)
      }
    }
  } else {                                      # BPWC_add > 0.66 (right subtree)
    if (BPWC_mul <= 0.31) {
      if (BPWC_dif <= 0.10) {
        if (BPWC_add <= 0.84) {
          if (BPWC_dif <= -0.25) return(10L) else return(11L)
        } else {                                # BPWC_add > 0.84
          if (BPWC_mul <= 0.24) return(12L) else return(13L)
        }
      } else {                                  # BPWC_dif > 0.10
        if (BPWC_mul <= 0.17) {
          if (BPWC_dif <= 0.56) return(14L) else return(15L)
        } else {                                # BPWC_mul > 0.17
          if (BPWC_mul <= 0.22) return(16L) else return(17L)
        }
      }
    } else {                                    # BPWC_mul > 0.31
      if (BPWC_dif <= 0.44) {
        if (BPWC_mul <= 0.38) {
          if (BPWC_dif <= 0.23) return(18L) else return(19L)
        } else {                                # BPWC_mul > 0.38
          if (BPWC_dif <= 0.08) return(20L) else return(21L)
        }
      } else {                                  # BPWC_dif > 0.44
        return(22L)
      }
    }
  }
}

#' Predict MetS risk in [0,1] from the full B9 decision tree (22 leaves).
#'
#' @param BPWC_add scalar, BP + WC
#' @param BPWC_mul scalar, BP * WC
#' @param BPWC_dif scalar, BP - WC
#' @return predicted MetS probability in [0,1]
b9_tree_predict <- function(BPWC_add, BPWC_mul, BPWC_dif) {
  leaf <- b9_tree_leaf(BPWC_add, BPWC_mul, BPWC_dif)
  risk <- B9_LEAF_RISK[[leaf]] * B9_PREVALENCE
  min(max(risk, 0), 1)
}

#' Vectorized B9 tree predictor for a data frame
#'
#' @param df data frame with columns BPWC_add, BPWC_mul, BPWC_dif
#' @return numeric vector of predicted probabilities
b9_tree_predict_vec <- function(df) {
  mapply(b9_tree_predict, df$BPWC_add, df$BPWC_mul, df$BPWC_dif)
}
