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

# LMF column positions per NCHS 2019 PUBLIC-USE release.
# Verified against actual file structure (files are 46-48 chars wide):
#   1-14: PUBLICID (SEQN, right-justified, space-padded)
#   15:   ELIGSTAT  (1=eligible, 2=under 18, 3=ineligible)
#   16:   MORTSTAT  (0=alive, 1=deceased, .=missing)
#   17-19: UCOD_LEADING (NCHS 113-cause recode, 3-char code; '.' = missing)
#   20:   DIABETES  (0/1, '.'=missing)
#   21:   HYPERTEN  (0/1, '.'=missing)
#   22:   DODQTR    (1-4, ' '=missing)
#   23-26: DODYEAR  (4-digit year, suppressed for privacy in public file)
#   27-34: WGT_NEW  (8-char, suppressed in public file)
#   35-42: SA_WGT_NEW (8-char, suppressed)
#   43-45: PERMTH_INT (3-char, person-months from interview to censoring)
#   46-48: PERMTH_EXM (3-char, person-months from MEC exam to censoring)
# Note: WGT_NEW, SA_WGT_NEW, DODQTR, DODYEAR are suppressed (spaces) in the
#       public-use file to protect respondent privacy.
lmf_widths <- fwf_positions(
  start = c(1,  15, 16, 17, 20, 21, 22, 23, 27, 35, 43, 46),
  end   = c(14, 15, 16, 19, 20, 21, 22, 26, 34, 42, 45, 48),
  col_names = c(
    "SEQN",            # NHANES Sequence ID (string, 14 char)
    "ELIGSTAT",        # Mortality eligibility
    "MORTSTAT",        # Vital status (0=alive, 1=deceased)
    "UCOD_LEADING",    # Underlying cause of death (NCHS 113 recode)
    "DIABETES",        # Diabetes flagged as contributing cause
    "HYPERTEN",        # Hypertension flagged as contributing cause
    "DODQTR",          # Quarter of death (suppressed in public file)
    "DODYEAR",         # Year of death (suppressed in public file)
    "WGT_NEW",         # MEC weight adjusted for mortality (suppressed)
    "SA_WGT_NEW",      # Subsample weight adjusted (suppressed)
    "PERMTH_INT",      # Person-months from interview to death/censor
    "PERMTH_EXM"       # Person-months from MEC exam to death/censor
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

# Coerce numeric columns (replace '.' and blank missing codes with NA)
for (col in c("ELIGSTAT", "MORTSTAT", "DIABETES", "HYPERTEN", "DODQTR",
              "DODYEAR", "PERMTH_INT", "PERMTH_EXM")) {
  v <- trimws(lmf_all[[col]])
  v[v == "." | v == ""] <- NA
  lmf_all[[col]] <- suppressWarnings(as.integer(v))
}
# UCOD_LEADING is a 3-char code string ("001"..."010"); coerce to integer.
# NCHS 113-cause recode summary (top 10 leading causes):
#   1=Heart disease, 2=Malignant neoplasms, 3=Chronic lower respiratory,
#   4=Accidents, 5=Cerebrovascular, 6=Alzheimer's, 7=Diabetes mellitus,
#   8=Influenza/pneumonia, 9=Nephritis, 10=All others
{
  v <- trimws(lmf_all[["UCOD_LEADING"]])
  v[v == "." | v == ""] <- NA
  lmf_all[["UCOD_LEADING"]] <- suppressWarnings(as.integer(v))
}
# WGT_NEW and SA_WGT_NEW are suppressed in the public file; keep as NA
lmf_all$WGT_NEW    <- NA_real_
lmf_all$SA_WGT_NEW <- NA_real_
lmf_all$SEQN <- as.integer(trimws(lmf_all$SEQN))

saveRDS(lmf_all, "data/raw/lmf/lmf_combined.rds")

message(sprintf("\nLMF combined: %d rows", nrow(lmf_all)))
message(sprintf("  Eligible for follow-up (ELIGSTAT=1): %d",
                sum(lmf_all$ELIGSTAT == 1, na.rm = TRUE)))
message(sprintf("  Deceased (MORTSTAT=1): %d",
                sum(lmf_all$MORTSTAT == 1, na.rm = TRUE)))
