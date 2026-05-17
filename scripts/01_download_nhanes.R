#!/usr/bin/env Rscript
# scripts/01_download_nhanes.R
#
# Downloads NHANES continuous cycles 1999-2018 via the nhanesA package.
# Required tables per cycle: DEMO, BMX, BPX, GLU, TRIGLY (or TCHOL), HDL,
# MCQ, DIQ, SMQ, BPQ.
#
# Output: one .rds per (cycle, table) under data/raw/nhanes/.
# Idempotent: skips a file if it already exists.

# Activate renv-installed library
.libPaths("/home/po/projects/work/longitudinal-mets-validation/renv/library/R-4.3/x86_64-pc-linux-gnu")

suppressMessages({
  library(nhanesA)
  library(purrr)
})

# Cycle letter suffixes per NHANES naming convention
# A=1999-2000, B=2001-02, C=2003-04, D=2005-06, E=2007-08,
# F=2009-10, G=2011-12, H=2013-14, I=2015-16, J=2017-18
cycles <- c("A", "B", "C", "D", "E", "F", "G", "H", "I", "J")

# Tables required for MetS / cohort construction + risk-score computation
required_tables <- c(
  "DEMO",   # demographics (age, sex, race, weights, PSU, strata)
  "BMX",    # body measurements (BMI, waist, height, weight)
  "BPX",    # blood pressure
  "GLU",    # plasma fasting glucose + fasting subsample weight
  "TRIGLY", # triglycerides (fasting subsample)
  "HDL",    # HDL cholesterol
  "TCHOL",  # total cholesterol (for PCE + Framingham)
  "MCQ",    # medical conditions (prior CVD)
  "DIQ",    # diabetes questionnaire (prior T2D)
  "SMQ",    # smoking
  "BPQ"     # BP/cholesterol questionnaire (medication)
)

dir.create("data/raw/nhanes", recursive = TRUE, showWarnings = FALSE)

download_cycle <- function(cycle) {
  # Naming: cycle A has no suffix; B-J use _B through _J
  cycle_suffix <- if (cycle == "A") "" else paste0("_", cycle)
  for (table in required_tables) {
    table_name <- paste0(table, cycle_suffix)
    out_path <- file.path("data/raw/nhanes", sprintf("%s_%s.rds", table, cycle))
    if (file.exists(out_path)) {
      message(sprintf("  [skip] %s (already exists)", out_path))
      next
    }
    message(sprintf("Fetching %s ...", table_name))
    df <- tryCatch(
      nhanes(table_name),
      error = function(e) {
        message(sprintf("  FAILED: %s", e$message))
        NULL
      }
    )
    if (!is.null(df)) {
      saveRDS(df, out_path)
      message(sprintf("  -> %s (n=%d)", out_path, nrow(df)))
    }
  }
}

walk(cycles, download_cycle)
message("\nDone. Files in data/raw/nhanes/:")
files <- list.files("data/raw/nhanes")
message(sprintf("  %d files total", length(files)))
