#!/usr/bin/env Rscript
# scripts/02_link_mortality.R
#
# Downloads the public-use NHANES Linked Mortality File (2019 release) per
# NHANES cycle (1999-2018) and parses the fixed-width files into a combined
# tibble keyed by SEQN.
#
# Output: data/raw/lmf/lmf_combined.rds

.libPaths("/home/po/projects/work/longitudinal-mets-validation/renv/library/R-4.3/x86_64-pc-linux-gnu")

suppressMessages({
  library(readr)
})

lmf_base <- "https://ftp.cdc.gov/pub/Health_Statistics/NCHS/datalinkage/linked_mortality"

# 2019 release filenames per NHANES cycle (verify against current NCHS docs)
lmf_files <- c(
  "NHANES_1999_2000_MORT_2019_PUBLIC.dat",
  "NHANES_2001_2002_MORT_2019_PUBLIC.dat",
  "NHANES_2003_2004_MORT_2019_PUBLIC.dat",
  "NHANES_2005_2006_MORT_2019_PUBLIC.dat",
  "NHANES_2007_2008_MORT_2019_PUBLIC.dat",
  "NHANES_2009_2010_MORT_2019_PUBLIC.dat",
  "NHANES_2011_2012_MORT_2019_PUBLIC.dat",
  "NHANES_2013_2014_MORT_2019_PUBLIC.dat",
  "NHANES_2015_2016_MORT_2019_PUBLIC.dat",
  "NHANES_2017_2018_MORT_2019_PUBLIC.dat"
)

dir.create("data/raw/lmf", recursive = TRUE, showWarnings = FALSE)

# Download missing files
for (f in lmf_files) {
  url <- file.path(lmf_base, f)
  dest <- file.path("data/raw/lmf", f)
  if (file.exists(dest)) {
    message(sprintf("  [skip] %s exists", basename(dest)))
    next
  }
  message(sprintf("Downloading %s ...", basename(dest)))
  tryCatch(
    download.file(url, dest, mode = "wb", quiet = TRUE),
    error = function(e) {
      message(sprintf("  FAILED: %s", e$message))
    }
  )
}

# LMF column positions per NCHS 2019 documentation
# Reference: NCHS Public-Use Linked Mortality Files documentation, 2019 release
# https://www.cdc.gov/nchs/data-linkage/mortality-public.htm
# Verify these positions against the current PDF documentation before running.
lmf_widths <- fwf_positions(
  start = c(1,  15, 16, 17, 21, 22, 23, 24, 28, 49, 60, 71, 73, 84, 95),
  end   = c(14, 15, 16, 20, 21, 22, 23, 27, 47, 59, 70, 72, 83, 94, 95),
  col_names = c(
    "SEQN",            # NHANES Sequence ID (string, 14 char)
    "ELIGSTAT",        # Mortality eligibility (1=eligible, 2=under 18 at survey, 3=ineligible)
    "MORTSTAT",        # Vital status (0=alive, 1=deceased)
    "UCOD_LEADING",    # Underlying cause of death (3-char NCHS code)
    "DIABETES",        # Diabetes flagged as cause (0/1)
    "HYPERTEN",        # Hypertension flagged as cause (0/1)
    "DODQTR",          # Quarter of death (1-4)
    "DODYEAR",         # Year of death (4-digit)
    "WGT_NEW",         # MEC sample weight adjusted for mortality follow-up
    "SA_WGT_NEW",      # Subsample weight adjusted
    "PERMTH_INT",      # Person-months from interview to death/censor
    "PERMTH_EXM",      # Person-months from MEC exam to death/censor
    "MORTSRCE_NDI",    # Source: NDI match
    "MORTSRCE_CMS",    # Source: CMS Medicare match
    "MORTSRCE_SSA"     # Source: SSA match
  )
)

# NOTE on widths: NCHS publishes the exact column positions in a SAS or R
# example file alongside each release. If parsing produces nonsense values
# (e.g., MORTSTAT outside {0,1}), re-check positions against the current
# documentation PDF.

# Parse + stack
lmf_list <- lapply(file.path("data/raw/lmf", lmf_files), function(p) {
  if (!file.exists(p)) {
    message(sprintf("  [missing] %s", basename(p)))
    return(NULL)
  }
  tryCatch(
    read_fwf(p, col_positions = lmf_widths, col_types = cols(.default = "c")),
    error = function(e) {
      message(sprintf("  parse error %s: %s", basename(p), e$message))
      NULL
    }
  )
})
lmf_list <- Filter(Negate(is.null), lmf_list)
lmf_all <- do.call(rbind, lmf_list)

# Coerce numeric columns
for (col in c("ELIGSTAT", "MORTSTAT", "DIABETES", "HYPERTEN", "DODQTR",
              "DODYEAR", "WGT_NEW", "SA_WGT_NEW", "PERMTH_INT", "PERMTH_EXM",
              "MORTSRCE_NDI", "MORTSRCE_CMS", "MORTSRCE_SSA")) {
  lmf_all[[col]] <- as.integer(lmf_all[[col]])
}
lmf_all$SEQN <- as.integer(lmf_all$SEQN)

saveRDS(lmf_all, "data/raw/lmf/lmf_combined.rds")

message(sprintf("\nLMF combined: %d rows", nrow(lmf_all)))
message(sprintf("  Eligible for follow-up (ELIGSTAT=1): %d",
                sum(lmf_all$ELIGSTAT == 1, na.rm = TRUE)))
message(sprintf("  Deceased (MORTSTAT=1): %d",
                sum(lmf_all$MORTSTAT == 1, na.rm = TRUE)))
