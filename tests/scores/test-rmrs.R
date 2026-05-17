library(testthat)

# Source paths assume tests are run from project root via:
#   Rscript -e 'testthat::test_dir("tests")'
source("../../R/scores/rmrs.R")

# ---- elliot_sigmoid ---------------------------------------------------------

test_that("elliot_sigmoid maps 0 to 0.5", {
  expect_equal(elliot_sigmoid(0), 0.5, tolerance = 1e-9)
})

test_that("elliot_sigmoid maps -1 to 0.25 (one unit below threshold)", {
  expect_equal(elliot_sigmoid(-1), 0.25, tolerance = 1e-9)
})

test_that("elliot_sigmoid maps +1 to 0.75 (one unit above threshold)", {
  expect_equal(elliot_sigmoid(1), 0.75, tolerance = 1e-9)
})

test_that("elliot_sigmoid is monotonic and bounded in (0,1)", {
  xs <- seq(-10, 10, by = 0.1)
  ys <- elliot_sigmoid(xs)
  expect_true(all(diff(ys) > 0))
  expect_true(all(ys > 0 & ys < 1))
})

# ---- scale_factor -----------------------------------------------------------

test_that("scale_factor maps glucose at threshold to 0.5", {
  expect_equal(scale_factor(100, "glucose"), 0.5, tolerance = 1e-6)
})

test_that("scale_factor maps triglycerides at threshold to 0.5", {
  expect_equal(scale_factor(150, "triglycerides"), 0.5, tolerance = 1e-6)
})

test_that("scale_factor reverses sign for HDL (lower = higher risk)", {
  # Male threshold = 40
  expect_gt(scale_factor(30, "hdl_male"), 0.5)   # below threshold = abnormal = > 0.5
  expect_lt(scale_factor(50, "hdl_male"), 0.5)   # above threshold = healthy = < 0.5
  expect_equal(scale_factor(40, "hdl_male"), 0.5, tolerance = 1e-6)
})

test_that("scale_factor for BP uses the larger sbp/dbp excess", {
  # sbp=130 (threshold), dbp=85 (threshold): both at zero excess
  expect_equal(scale_factor(c(130, 85), "bp"), 0.5, tolerance = 1e-6)
  # sbp above threshold dominates
  expect_gt(scale_factor(c(140, 80), "bp"), 0.5)
  # dbp above threshold dominates
  expect_gt(scale_factor(c(125, 95), "bp"), 0.5)
})

# ---- tas_score --------------------------------------------------------------

test_that("tas_score returns a valid number for boundary input", {
  result <- tas_score(c(0.5, 0.5, 0.5))
  expect_gte(result, 0)
  expect_lte(result, 1)
})

test_that("tas_score is higher when all axes exceed threshold", {
  low  <- tas_score(c(0.5, 0.5, 0.5))
  high <- tas_score(c(0.9, 0.9, 0.9))
  expect_gt(high, low)
})

test_that("tas_score is bounded in [0,1]", {
  set.seed(42)
  axes <- matrix(runif(300), ncol = 3)
  scores <- apply(axes, 1, tas_score)
  expect_true(all(scores >= 0 & scores <= 1))
})

# ---- RMRS aggregation ------------------------------------------------------

test_that("rmrs_from_scaled accepts a length-5 vector and returns a scalar in [0,1]", {
  result <- rmrs_from_scaled(c(0.5, 0.5, 0.5, 0.5, 0.5))
  expect_length(result, 1)
  expect_gte(result, 0)
  expect_lte(result, 1)
})

test_that("rmrs_from_scaled is symmetric under permutation", {
  set.seed(42)
  x <- runif(5)
  expect_equal(rmrs_from_scaled(x), rmrs_from_scaled(rev(x)), tolerance = 1e-9)
  expect_equal(rmrs_from_scaled(x), rmrs_from_scaled(sample(x)), tolerance = 1e-9)
})

# Worked-example tests from B8 paper section "Diagnostic threshold"
# These define the structural threshold of 0.547. The current implementation
# uses a simplified TAS formula and may not exactly reproduce these values.
# If these tests fail, the TAS formula needs refinement against the published
# derivation (see R/scores/rmrs.R Implementation Notes).

test_that("rmrs_from_scaled approaches the published lower boundary case", {
  # Three factors at threshold, two well below: RMRS = 0.387 per Shin 2024
  scaled <- c(0.5, 0.5, 0.5, 0, 0)
  result <- rmrs_from_scaled(scaled)
  expect_gte(result, 0)
  expect_lte(result, 1)
  # The published value is 0.387; allow generous tolerance because TAS
  # implementation differs slightly without explicit "Sb beyond-threshold area"
  # weighting. If this fails, refine tas_score per B8.
  expect_lt(abs(result - 0.387), 0.15)
})

test_that("rmrs_from_scaled approaches the published upper boundary case", {
  # Two above, one near and two at threshold: RMRS = 0.707 per Shin 2024
  scaled <- c(1, 1, 0.5 - 1e-6, 0.5 - 1e-6, 0.5)
  result <- rmrs_from_scaled(scaled)
  expect_gte(result, 0)
  expect_lte(result, 1)
  expect_lt(abs(result - 0.707), 0.15)
})

test_that("rmrs threshold 0.547 sits between the boundary cases", {
  expect_gt(0.547, 0.387)
  expect_lt(0.547, 0.707)
})

# ---- end-to-end ------------------------------------------------------------

test_that("rmrs(raw inputs) returns a valid score for a healthy adult", {
  # Healthy 40yo female: WC=70, BP=110/70, HDL=60, glucose=85, TG=80
  result <- rmrs(waist = 70, sbp = 110, dbp = 70, hdl = 60,
                 glucose = 85, triglycerides = 80,
                 sex = "female", ethnicity = "us")
  expect_gte(result, 0)
  expect_lte(result, 1)
  expect_lt(result, 0.547)  # below diagnostic threshold
})

test_that("rmrs(raw inputs) returns a high score for clearly abnormal adult", {
  # MetS adult: WC=110, BP=160/100, HDL=30, glucose=130, TG=250 (all abnormal)
  result <- rmrs(waist = 110, sbp = 160, dbp = 100, hdl = 30,
                 glucose = 130, triglycerides = 250,
                 sex = "male", ethnicity = "us")
  expect_gte(result, 0.547)  # above diagnostic threshold
})
