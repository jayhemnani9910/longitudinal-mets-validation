library(testthat)
source("R/scores/b9_features.R")
source("R/scores/b9_tree.R")

# ---- b9_features -----------------------------------------------------------

test_that("b9_features returns the expected 5-element structure", {
  feat <- b9_features(waist = 100, sbp = 130, dbp = 80,
                       sex = "male", ethnicity = "us")
  expect_named(feat, c("WC", "BP", "BPWC_add", "BPWC_mul", "BPWC_dif"))
  expect_true(all(is.finite(unlist(feat))))
})

test_that("b9_features matches Figure 7 worked example (woman, sBP=140, dBP=90, waist=89)", {
  feat <- b9_features(waist = 89, sbp = 140, dbp = 90,
                       sex = "female", ethnicity = "us")
  # Per B9 paper Fig 7: BP = 0.84, WC = 0.66
  # Then BPWC_add = 1.50, BPWC_mul = 0.55, BPWC_dif = 0.18
  expect_equal(feat$BP, 0.84, tolerance = 0.05)
  expect_equal(feat$WC, 0.66, tolerance = 0.05)
  expect_equal(feat$BPWC_add, 1.50, tolerance = 0.05)
  expect_equal(feat$BPWC_mul, 0.55, tolerance = 0.05)
  expect_equal(feat$BPWC_dif, 0.18, tolerance = 0.05)
})

# ---- b9_tree_predict --------------------------------------------------------

test_that("b9_tree_predict returns a probability in [0,1]", {
  feat <- b9_features(waist = 100, sbp = 130, dbp = 80, sex = "male")
  p <- with(feat, b9_tree_predict(BPWC_add, BPWC_mul, BPWC_dif))
  expect_gte(p, 0)
  expect_lte(p, 1)
})

test_that("Figure 7 worked example yields probability ~0.31", {
  feat <- b9_features(waist = 89, sbp = 140, dbp = 90,
                       sex = "female", ethnicity = "us")
  p <- with(feat, b9_tree_predict(BPWC_add, BPWC_mul, BPWC_dif))
  # B9 Figure 7 reports diagnosis: MetS, risk: 0.31
  expect_equal(p, 0.31, tolerance = 0.05)
})

test_that("safety zone (BPWC_add <= 0.66) returns population baseline", {
  # Use a low-WC + low-BP subject
  feat <- b9_features(waist = 70, sbp = 110, dbp = 65, sex = "female")
  expect_lte(feat$BPWC_add, 0.66)
  p <- with(feat, b9_tree_predict(BPWC_add, BPWC_mul, BPWC_dif))
  expect_equal(p, 0.137, tolerance = 1e-6)
})

test_that("b9_tree_predict gives higher risk for clearly abnormal inputs", {
  feat_low  <- b9_features(waist = 70,  sbp = 110, dbp = 65, sex = "male")
  feat_high <- b9_features(waist = 110, sbp = 160, dbp = 100, sex = "male")
  p_low  <- with(feat_low,  b9_tree_predict(BPWC_add, BPWC_mul, BPWC_dif))
  p_high <- with(feat_high, b9_tree_predict(BPWC_add, BPWC_mul, BPWC_dif))
  expect_gt(p_high, p_low)
})

# ---- vectorized form -------------------------------------------------------

test_that("b9_tree_predict_vec works on a data.frame", {
  df <- data.frame(
    BPWC_add = c(0.3, 1.50, 1.20),
    BPWC_mul = c(0.05, 0.55, 0.20),
    BPWC_dif = c(0.1, 0.18, 0.05)
  )
  preds <- b9_tree_predict_vec(df)
  expect_length(preds, 3)
  expect_true(all(preds >= 0 & preds <= 1))
  expect_equal(preds[1], 0.137, tolerance = 1e-6)  # safety zone
  expect_equal(preds[2], 0.31,  tolerance = 0.05)  # Figure 7 example
})
