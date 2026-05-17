library(testthat)
source("../../R/scores/pce.R")

test_that("pce returns plausible 10-year risk for a low-risk young adult", {
  # Healthy 45yo white woman, no smoker, no diabetes, normal lipids and BP
  risk <- pce(
    age = 45, sex = "female", race = "white",
    total_chol = 180, hdl = 60, sbp = 110,
    bp_treated = FALSE, smoker = FALSE, diabetes = FALSE
  )
  expect_gte(risk, 0)
  expect_lte(risk, 1)
  expect_lt(risk, 0.05)  # under 5% is the expected ballpark for low-risk young adult
})

test_that("pce returns plausible 10-year risk for a high-risk older adult", {
  # 65yo white man, smoker, diabetic, high cholesterol, high BP
  risk <- pce(
    age = 65, sex = "male", race = "white",
    total_chol = 250, hdl = 35, sbp = 160,
    bp_treated = TRUE, smoker = TRUE, diabetes = TRUE
  )
  expect_gte(risk, 0)
  expect_lte(risk, 1)
  expect_gt(risk, 0.20)  # expect >20% 10-year ASCVD risk
})

test_that("pce handles NA gracefully", {
  expect_true(is.na(pce(
    age = NA, sex = "male", race = "white",
    total_chol = 200, hdl = 40, sbp = 130,
    bp_treated = FALSE, smoker = FALSE, diabetes = FALSE
  )))
})

test_that("pce yields different scores for black vs white men with same inputs", {
  # The PCE has separate equations for white and black race; scores should differ
  base <- list(
    age = 55, sex = "male",
    total_chol = 200, hdl = 45, sbp = 140,
    bp_treated = FALSE, smoker = FALSE, diabetes = FALSE
  )
  risk_w <- do.call(pce, c(base, race = "white"))
  risk_b <- do.call(pce, c(base, race = "black"))
  expect_false(isTRUE(all.equal(risk_w, risk_b)))
})

test_that("pce_vec works on a data frame", {
  df <- data.frame(
    age = c(45, 60),
    sex = c("female", "male"),
    race = c("white", "black"),
    total_chol = c(180, 220),
    hdl = c(60, 40),
    sbp = c(110, 145),
    bp_treatment = c(FALSE, TRUE),
    current_smoker = c(FALSE, TRUE),
    prior_t2d = c(FALSE, FALSE)
  )
  risks <- pce_vec(df)
  expect_length(risks, 2)
  expect_true(all(risks >= 0 & risks <= 1))
})
