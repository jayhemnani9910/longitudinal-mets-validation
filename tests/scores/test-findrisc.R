library(testthat)
source("../../R/scores/findrisc.R")

test_that("findrisc returns 0 for the lowest-risk profile", {
  pts <- findrisc(
    age = 30, bmi = 22, waist = 70, sex = "male",
    physical_active = TRUE, daily_vegetables = TRUE,
    bp_treatment = FALSE, prior_high_glucose = FALSE,
    family_history = "none"
  )
  expect_equal(pts, 0L)
})

test_that("findrisc returns 26 for the highest-risk profile (man)", {
  # Max: age 4 + BMI 3 + waist 4 + no exercise 2 + no veg 1 + BP-tx 2
  #      + prior hi glucose 5 + first-degree family 5 = 26
  pts <- findrisc(
    age = 70, bmi = 35, waist = 105, sex = "male",
    physical_active = FALSE, daily_vegetables = FALSE,
    bp_treatment = TRUE, prior_high_glucose = TRUE,
    family_history = "first_degree"
  )
  expect_equal(pts, 26L)
})

test_that("findrisc sex-specific waist scoring", {
  # Female waist 85 (>= 80 but < 88 = 3 pts)
  pts_f <- findrisc(
    age = 40, bmi = 23, waist = 85, sex = "female",
    physical_active = TRUE, daily_vegetables = TRUE,
    bp_treatment = FALSE, prior_high_glucose = FALSE,
    family_history = "none"
  )
  # Male waist 85 (< 94 = 0 pts)
  pts_m <- findrisc(
    age = 40, bmi = 23, waist = 85, sex = "male",
    physical_active = TRUE, daily_vegetables = TRUE,
    bp_treatment = FALSE, prior_high_glucose = FALSE,
    family_history = "none"
  )
  expect_equal(pts_f - pts_m, 3L)
})

test_that("findrisc components add as expected for a known example", {
  # 55yo woman, BMI 28, waist 90cm, no exercise, no veggies,
  # on BP meds, no prior hi glucose, no family history
  # Expected: 3 (age) + 1 (BMI) + 4 (waist >= 88 female) + 2 (no exercise)
  #         + 1 (no veg) + 2 (BP tx) + 0 + 0 = 13
  pts <- findrisc(
    age = 55, bmi = 28, waist = 90, sex = "female",
    physical_active = FALSE, daily_vegetables = FALSE,
    bp_treatment = TRUE, prior_high_glucose = FALSE,
    family_history = "none"
  )
  expect_equal(pts, 13L)
})

test_that("findrisc handles NA gracefully", {
  expect_true(is.na(findrisc(
    age = NA, bmi = 25, waist = 90, sex = "male",
    physical_active = TRUE, daily_vegetables = TRUE,
    bp_treatment = FALSE, prior_high_glucose = FALSE,
    family_history = "none"
  )))
})

test_that("findrisc_vec works on a data frame", {
  df <- data.frame(
    age = c(45, 60),
    bmi = c(24, 32),
    waist_cm = c(85, 105),
    sex = c("female", "male"),
    physical_active = c(TRUE, FALSE),
    daily_vegetables = c(TRUE, FALSE),
    bp_treatment = c(FALSE, TRUE),
    prior_high_glucose = c(FALSE, TRUE),
    family_history_dm = c("none", "first_degree")
  )
  pts <- findrisc_vec(df)
  expect_length(pts, 2)
  expect_true(all(pts >= 0 & pts <= 26))
  expect_lt(pts[1], pts[2])  # second row is clearly higher risk
})
