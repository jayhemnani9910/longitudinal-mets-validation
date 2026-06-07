library(testthat)
source("../../R/scores/rmrs.R")        # provides elliot_sigmoid + scale_factor
source("../../R/scores/b9_features.R") # depends on the above
source("../../R/scores/b9_tree.R")

# ---- b9_features -----------------------------------------------------------

test_that("b9_features returns the expected 5-element structure", {
  feat <- b9_features(waist = 100, sbp = 130, dbp = 80,
                       sex = "male", ethnicity = "us")
  expect_named(feat, c("WC", "BP", "BPWC_add", "BPWC_mul", "BPWC_dif"))
  expect_true(all(is.finite(unlist(feat))))
})

test_that("b9_features matches Figure 7 worked example (woman, sBP=140, dBP=90, waist=89, Korean cohort)", {
  # B9 was trained on Korean data, so the Figure 7 worked example uses the
  # Korean WC threshold (Female >= 85 cm), NOT the US threshold (Female >= 88).
  feat <- b9_features(waist = 89, sbp = 140, dbp = 90,
                       sex = "female", ethnicity = "kr")
  # Per B9 paper Fig 7: BP = 0.84, WC = 0.66
  # Then BPWC_add = 1.50, BPWC_mul = 0.55, BPWC_dif = 0.18
  expect_equal(feat$BP, 0.84, tolerance = 0.05)
  expect_equal(feat$WC, 0.66, tolerance = 0.05)
  expect_equal(feat$BPWC_add, 1.50, tolerance = 0.05)
  expect_equal(feat$BPWC_mul, 0.55, tolerance = 0.05)
  expect_equal(feat$BPWC_dif, 0.18, tolerance = 0.05)
})

# ---- b9_tree_predict (full 22-leaf S1-Table reconstruction) ----------------

test_that("b9_tree_predict returns a probability in [0,1]", {
  feat <- b9_features(waist = 100, sbp = 130, dbp = 80, sex = "male", ethnicity = "kr")
  p <- with(feat, b9_tree_predict(BPWC_add, BPWC_mul, BPWC_dif))
  expect_gte(p, 0)
  expect_lte(p, 1)
})

test_that("Figure 7 worked example routes to Leaf 21 with probability ~0.31", {
  feat <- b9_features(waist = 89, sbp = 140, dbp = 90,
                       sex = "female", ethnicity = "kr")
  leaf <- with(feat, b9_tree_leaf(BPWC_add, BPWC_mul, BPWC_dif))
  p    <- with(feat, b9_tree_predict(BPWC_add, BPWC_mul, BPWC_dif))
  expect_equal(leaf, 21L)
  # B9 Figure 7 reports diagnosis: MetS, risk: 0.31 (= 2.25 * 0.137)
  expect_equal(p, 0.31, tolerance = 0.02)
})

test_that("safety zone (BPWC_add <= 0.66) returns a low risk below the MetS prevalence", {
  # A low-WC + low-BP subject lands in the left subtree.
  feat <- b9_features(waist = 70, sbp = 110, dbp = 65, sex = "female", ethnicity = "kr")
  expect_lte(feat$BPWC_add, 0.66)
  p <- with(feat, b9_tree_predict(BPWC_add, BPWC_mul, BPWC_dif))
  expect_lt(p, 0.137)
})

test_that("a clearly high-risk subject routes into the high-risk red leaves", {
  feat <- b9_features(waist = 120, sbp = 175, dbp = 110, sex = "male", ethnicity = "kr")
  leaf <- with(feat, b9_tree_leaf(BPWC_add, BPWC_mul, BPWC_dif))
  p    <- with(feat, b9_tree_predict(BPWC_add, BPWC_mul, BPWC_dif))
  expect_true(leaf %in% c(18L, 20L, 21L, 22L, 13L, 17L))
  expect_gt(p, 0.30)
})

test_that("b9_tree_predict gives higher risk for clearly abnormal inputs", {
  feat_low  <- b9_features(waist = 70,  sbp = 110, dbp = 65,  sex = "male", ethnicity = "kr")
  feat_high <- b9_features(waist = 110, sbp = 165, dbp = 105, sex = "male", ethnicity = "kr")
  p_low  <- with(feat_low,  b9_tree_predict(BPWC_add, BPWC_mul, BPWC_dif))
  p_high <- with(feat_high, b9_tree_predict(BPWC_add, BPWC_mul, BPWC_dif))
  expect_gt(p_high, p_low)
})

test_that("the full reconstruction produces many distinct risk levels (not 3)", {
  # The earlier partial reconstruction emitted only 3 values; the full S1 tree
  # spans up to 17 distinct leaf risks. Sweep the BP x WC feature plane.
  g <- expand.grid(BP = seq(0, 1, 0.04), WC = seq(0, 1, 0.04))
  preds <- mapply(function(bp, wc) b9_tree_predict(bp + wc, bp * wc, bp - wc),
                  g$BP, g$WC)
  expect_gte(length(unique(round(preds, 5))), 10)
})

# ---- vectorized form -------------------------------------------------------

test_that("b9_tree_predict_vec works on a data.frame", {
  df <- data.frame(
    BPWC_add = c(0.30, 1.50, 1.20),
    BPWC_mul = c(0.05, 0.55, 0.20),
    BPWC_dif = c(0.10, 0.18, 0.05)
  )
  preds <- b9_tree_predict_vec(df)
  expect_length(preds, 3)
  expect_true(all(preds >= 0 & preds <= 1))
  expect_lt(preds[1], 0.137)               # safety zone -> low risk
  expect_equal(preds[2], 0.31, tolerance = 0.02)  # Figure 7 example, Leaf 21
})
