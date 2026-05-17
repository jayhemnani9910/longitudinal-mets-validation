library(testthat)
source("../../R/scores/framingham.R")

test_that("framingham_2008 returns plausible risk for low-risk young adult", {
  risk <- framingham_2008(
    age = 40, sex = "female",
    total_chol = 180, hdl = 60, sbp = 110,
    bp_treated = FALSE, smoker = FALSE, diabetes = FALSE
  )
  expect_gte(risk, 0)
  expect_lte(risk, 1)
  expect_lt(risk, 0.05)
})

test_that("framingham_2008 returns plausible risk for high-risk older adult", {
  risk <- framingham_2008(
    age = 65, sex = "male",
    total_chol = 250, hdl = 35, sbp = 160,
    bp_treated = TRUE, smoker = TRUE, diabetes = TRUE
  )
  expect_gte(risk, 0)
  expect_lte(risk, 1)
  expect_gt(risk, 0.30)
})

test_that("framingham_2008 differentiates sex", {
  base <- list(
    age = 55,
    total_chol = 220, hdl = 45, sbp = 140,
    bp_treated = FALSE, smoker = FALSE, diabetes = FALSE
  )
  risk_m <- do.call(framingham_2008, c(base, sex = "male"))
  risk_f <- do.call(framingham_2008, c(base, sex = "female"))
  expect_false(isTRUE(all.equal(risk_m, risk_f)))
})

test_that("framingham_2008 handles NA gracefully", {
  expect_true(is.na(framingham_2008(
    age = NA, sex = "male",
    total_chol = 200, hdl = 40, sbp = 130,
    bp_treated = FALSE, smoker = FALSE, diabetes = FALSE
  )))
})

test_that("framingham_vec works on a data frame", {
  df <- data.frame(
    age = c(45, 60),
    sex = c("female", "male"),
    total_chol = c(180, 220),
    hdl = c(60, 40),
    sbp = c(110, 145),
    bp_treatment = c(FALSE, TRUE),
    current_smoker = c(FALSE, TRUE),
    prior_t2d = c(FALSE, FALSE)
  )
  risks <- framingham_vec(df)
  expect_length(risks, 2)
  expect_true(all(risks >= 0 & risks <= 1))
})
