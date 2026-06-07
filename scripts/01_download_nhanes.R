#!/usr/bin/env Rscript
# scripts/01_download_nhanes.R
#
# Downloads NHANES continuous cycles 1999-2018 via the nhanesA package.
# Required tables per cycle: DEMO, BMX, BPX, GLU, TRIGLY, HDL, TCHOL,
# MCQ, DIQ, SMQ, BPQ, PAQ.
#
# Output: one .rds per (cycle, table) under data/raw/nhanes/, named with the
# standardized "<TABLE>_<CYCLE>.rds" convention the rest of the pipeline expects.
# Idempotent: skips a file that already exists.
#
# Early-cycle note: NHANES standardized the fasting-lab table names
# (GLU, TRIGLY, HDL, TCHOL) only from 2005-2006 (cycle D). The 1999-2004 cycles
# (A, B, C) publish these labs under legacy names (mapped in legacy_lab_names).
# Requesting the modern names for those cycles returns NULL from nhanesA and
# would silently drop every fasting-subsample subject from the first three
# cycles, shrinking the analytic cohort from 17,031 to about 12,272 with no
# error. The integrity check at the end fails loudly if any expected file is
# missing, and fetch_table retries transient NULL returns.

.libPaths("renv/library/R-4.3/x86_64-pc-linux-gnu")

suppressMessages({
  library(nhanesA)
  library(purrr)
})

# Cycle letter suffixes per NHANES naming convention.
# A=1999-2000, B=2001-02, C=2003-04, D=2005-06, E=2007-08,
# F=2009-10, G=2011-12, H=2013-14, I=2015-16, J=2017-18
cycles <- c("A", "B", "C", "D", "E", "F", "G", "H", "I", "J")

# Tables required for MetS / cohort construction + risk-score computation.
required_tables <- c(
  "DEMO",   # demographics (age, sex, race, weights, PSU, strata)
  "BMX",    # body measurements (BMI, waist, height, weight)
  "BPX",    # blood pressure
  "GLU",    # plasma fasting glucose + fasting subsample weight
  "TRIGLY", # triglycerides (fasting subsample)
  "HDL",    # HDL cholesterol
  "TCHOL",  # total cholesterol (for PCE + Framingham)
  "MCQ",    # medical conditions (prior CVD, family history of diabetes)
  "DIQ",    # diabetes questionnaire (prior T2D, prediabetes history)
  "SMQ",    # smoking
  "BPQ",    # BP/cholesterol questionnaire (medication)
  "PAQ"     # physical activity questionnaire (FINDRISC activity item)
)

# Legacy NHANES table names for the 1999-2004 fasting labs. HDL and total
# cholesterol share the combined lipid panel (LAB13 / L13_B / L13_C) in these
# cycles; the downstream loader selects LBXTC (total cholesterol) and
# LBDHDL/LBXHDD (HDL) as needed.
legacy_lab_names <- list(
  GLU    = c(A = "LAB10AM", B = "L10AM_B", C = "L10AM_C"),
  TRIGLY = c(A = "LAB13AM", B = "L13AM_B", C = "L13AM_C"),
  HDL    = c(A = "LAB13",   B = "L13_B",   C = "L13_C"),
  TCHOL  = c(A = "LAB13",   B = "L13_B",   C = "L13_C")
)

# Resolve the NHANES table name to request for a (table, cycle) pair.
resolve_table_name <- function(table, cycle) {
  if (table %in% names(legacy_lab_names) && cycle %in% c("A", "B", "C")) {
    return(legacy_lab_names[[table]][[cycle]])
  }
  suffix <- if (cycle == "A") "" else paste0("_", cycle)
  paste0(table, suffix)
}

# Fetch one table, retrying on a NULL or empty return. nhanesA occasionally
# returns NULL on a transient server or network hiccup; without a retry a single
# blip would leave a file missing and abort the run at the integrity check.
fetch_table <- function(nhanes_name, attempts = 4L) {
  for (i in seq_len(attempts)) {
    df <- tryCatch(nhanes(nhanes_name), error = function(e) NULL)
    if (!is.null(df) && nrow(df) > 0) return(df)
    if (i < attempts) {
      message(sprintf("  retry %d/%d for %s ...", i, attempts, nhanes_name))
      Sys.sleep(3)
    }
  }
  NULL
}

dir.create("data/raw/nhanes", recursive = TRUE, showWarnings = FALSE)

download_cycle <- function(cycle) {
  for (table in required_tables) {
    out_path <- file.path("data/raw/nhanes", sprintf("%s_%s.rds", table, cycle))
    if (file.exists(out_path)) {
      message(sprintf("  [skip] %s (already exists)", out_path))
      next
    }
    nhanes_name <- resolve_table_name(table, cycle)
    message(sprintf("Fetching %s (cycle %s, requested as %s) ...",
                    table, cycle, nhanes_name))
    df <- fetch_table(nhanes_name)
    if (!is.null(df)) {
      saveRDS(df, out_path)
      message(sprintf("  -> %s (n=%d)", out_path, nrow(df)))
    } else {
      message(sprintf("  COULD NOT FETCH %s for cycle %s (requested as %s)",
                      table, cycle, nhanes_name))
    }
  }
}

walk(cycles, download_cycle)

# Integrity check: every (table, cycle) we need must now exist on disk. A silent
# NULL from nhanesA would otherwise shrink the cohort downstream with no error,
# so fail loudly and name exactly what is missing.
expected <- as.vector(outer(required_tables, cycles,
                            function(t, c) sprintf("%s_%s.rds", t, c)))
present  <- file.exists(file.path("data/raw/nhanes", expected))
if (!all(present)) {
  missing <- expected[!present]
  stop(sprintf(paste0("Download incomplete: %d of %d expected NHANES files are ",
                      "missing. Re-run to retry (downloads are idempotent), or ",
                      "check the table-name mapping for:\n  %s"),
               length(missing), length(expected), paste(missing, collapse = "\n  ")))
}
message(sprintf("\nDone. All %d expected NHANES files present in data/raw/nhanes/.",
                length(expected)))
