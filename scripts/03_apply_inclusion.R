#!/usr/bin/env Rscript
# scripts/03_apply_inclusion.R
#
# Builds the analysis cohort from raw NHANES + LMF by:
#   1. Loading and stacking per-cycle NHANES tables
#   2. Standardizing column names across cycles
#   3. Applying inclusion criteria from the OSF pre-registration
#   4. Linking mortality outcomes
#   5. Capping follow-up at 10 years
#   6. Constructing competing-risk event indicators
#
# Output: data/processed/analysis_cohort.rds

.libPaths("/home/po/projects/work/longitudinal-mets-validation/renv/library/R-4.3/x86_64-pc-linux-gnu")

suppressMessages({
  library(dplyr)
  library(tidyr)
})

source("R/utils/load_nhanes.R")

message("Loading NHANES tables ...")
demo   <- load_nhanes_table("DEMO")
bmx    <- load_nhanes_table("BMX")
bpx    <- load_nhanes_table("BPX")
glu    <- load_nhanes_table("GLU")
trigly <- load_nhanes_table("TRIGLY")
hdl    <- load_nhanes_table("HDL")
tchol  <- load_nhanes_table("TCHOL")
mcq    <- load_nhanes_table("MCQ")
diq    <- load_nhanes_table("DIQ")
smq    <- load_nhanes_table("SMQ")
bpq    <- load_nhanes_table("BPQ")
paq    <- load_nhanes_table("PAQ")

message("Standardizing column names ...")

# Demographics: harmonize race coding across cycles.
# RIDRETH3 (categories: 1=Mex Am, 2=Other Hisp, 3=NH White, 4=NH Black,
# 6=NH Asian, 7=Other) is available 2011+ (cycles G-J).
# RIDRETH1 (no separate Asian category) is available all cycles.
# We prefer RIDRETH3 when present; for earlier cycles, NH Asian is collapsed
# into Other.
# nhanesA converts factor codes to their labels when the package returns
# data frames, and our load_nhanes_table strips factor structure to character.
# So RIAGENDR arrives as "Male"/"Female" and RIDRETH* arrive as text labels.
# Map these back to the canonical factor levels expected downstream.

# Normalize gender: handles both numeric (1/2) and text label forms.
demo <- demo %>% mutate(
  sex_raw = case_when(
    RIAGENDR %in% c("1", "Male")   ~ "male",
    RIAGENDR %in% c("2", "Female") ~ "female",
    TRUE ~ NA_character_
  ),
  # Normalize race: RIDRETH3 preferred (has Asian category); fall back to RIDRETH1.
  # Both may arrive as text labels or as character "1"-"7".
  race_src = if ("RIDRETH3" %in% names(demo)) {
    dplyr::coalesce(RIDRETH3, RIDRETH1)
  } else {
    RIDRETH1
  },
  race_raw = case_when(
    race_src %in% c("1", "Mexican American")                      ~ "MexicanAm",
    race_src %in% c("2", "Other Hispanic")                         ~ "OtherHisp",
    race_src %in% c("3", "Non-Hispanic White")                     ~ "NHWhite",
    race_src %in% c("4", "Non-Hispanic Black")                     ~ "NHBlack",
    race_src %in% c("6", "Non-Hispanic Asian")                     ~ "NHAsian",
    race_src %in% c("7", "Other Race - Including Multi-Racial", "Other") ~ "Other",
    TRUE ~ NA_character_
  )
)

demo_sub <- demo %>% transmute(
  SEQN,
  cycle,
  age      = RIDAGEYR,
  sex      = factor(sex_raw, levels = c("male", "female")),
  race     = factor(race_raw, levels = c("MexicanAm", "OtherHisp", "NHWhite",
                                          "NHBlack", "NHAsian", "Other")),
  wt_mec   = WTMEC2YR,
  sdmvpsu  = SDMVPSU,
  sdmvstra = SDMVSTRA
)

# Anthropometry
bmx_sub <- bmx %>% transmute(
  SEQN,
  bmi      = BMXBMI,
  waist_cm = BMXWAIST,
  height_cm = BMXHT,
  weight_kg = BMXWT
)

# Blood pressure: NHANES reports up to 4 readings (BPXSY1-4 / BPXDI1-4).
# Use the mean of available non-zero readings (zero indicates "not done").
bpx_sub <- bpx %>%
  mutate(across(starts_with("BPXSY"), ~ na_if(.x, 0))) %>%
  mutate(across(starts_with("BPXDI"), ~ na_if(.x, 0))) %>%
  rowwise() %>%
  mutate(
    sbp = mean(c_across(starts_with("BPXSY")), na.rm = TRUE),
    dbp = mean(c_across(starts_with("BPXDI")), na.rm = TRUE)
  ) %>%
  ungroup() %>%
  transmute(SEQN, sbp, dbp)

# Glucose (fasting subsample)
glu_sub <- glu %>% transmute(
  SEQN,
  fasting_glucose = LBXGLU,
  wt_fast = WTSAF2YR
)

# Triglycerides + HDL
trigly_sub <- trigly %>% transmute(SEQN, triglycerides = LBXTR)

# HDL column name changed across cycles:
#   Cycles A-B (1999-2002): LBDHDL  (from LAB13 / L13_B)
#   Cycle  C   (2003-2004): LBXHDD  (from L13_C)
#   Cycles D+  (2005+):     LBDHDD  (standard)
# Coalesce all three column names into a unified variable.
{
  hdl_cols <- list(
    LBDHDD = if ("LBDHDD" %in% names(hdl)) hdl[["LBDHDD"]] else rep(NA_real_, nrow(hdl)),
    LBXHDD = if ("LBXHDD" %in% names(hdl)) hdl[["LBXHDD"]] else rep(NA_real_, nrow(hdl)),
    LBDHDL = if ("LBDHDL" %in% names(hdl)) hdl[["LBDHDL"]] else rep(NA_real_, nrow(hdl))
  )
  hdl[["hdl_combined"]] <- dplyr::coalesce(hdl_cols$LBDHDD, hdl_cols$LBXHDD, hdl_cols$LBDHDL)
}
hdl_sub <- hdl %>% transmute(SEQN, hdl = hdl_combined)
tchol_sub  <- tchol  %>% transmute(SEQN, total_chol    = LBXTC)

# Prior CVD: MCQ160B (CHF), MCQ160C (CHD), MCQ160E (MI), MCQ160F (stroke).
# nhanesA converts to factor labels; after load_nhanes_table these are characters.
# Numeric 1 = yes; text label = "Yes". Handle both forms.
is_yes <- function(x) x %in% c("1", "Yes", "YES", 1L)

mcq_sub <- mcq %>% transmute(
  SEQN,
  prior_cvd = is_yes(MCQ160B) | is_yes(MCQ160C) |
              is_yes(MCQ160E) | is_yes(MCQ160F),
  prior_cvd = replace_na(prior_cvd, FALSE),
  # Family history of diabetes (FINDRISC item). MCQ300C asks whether a close
  # biological relative was ever told they had diabetes. Available cycles D+
  # (2005-2006 onward); absent in 1999-2004, which fall back to no-history.
  # MCQ300C does not separate first- from second-degree relatives, so any
  # positive is mapped to the first-degree FINDRISC tier.
  fam_dm_yes = is_yes(get_or_na(mcq, "MCQ300C"))
)

# Prior T2D: DIQ010 == 1 means "doctor told you have diabetes"
# Prediabetes history: DIQ160 == 1 means "ever told you have prediabetes"
# (FINDRISC high-blood-glucose item). Available cycles D+ (2005-2006 onward).
diq_sub <- diq %>% transmute(
  SEQN,
  prior_t2d = is_yes(DIQ010),
  prior_t2d = replace_na(prior_t2d, FALSE),
  prediabetes_hx = is_yes(get_or_na(diq, "DIQ160")),
  prediabetes_hx = replace_na(prediabetes_hx, FALSE)
)

# Smoking: SMQ040 (current smoking) — 1="Every day", 2="Some days", 3="Not at all"
smq_sub <- smq %>% transmute(
  SEQN,
  current_smoker = SMQ040 %in% c("1", "2", "Every day", "Some days",
                                   "Every day,", "Some days,")
)

# BP medication: BPQ050A == 1 means "now taking medicine for high BP"
bpq_sub <- bpq %>% transmute(
  SEQN,
  bp_treatment = is_yes(BPQ050A),
  bp_treatment = replace_na(bp_treatment, FALSE)
)

# Physical activity (FINDRISC item). The PAQ instrument was redesigned in 2007,
# so the variable names differ by era:
#   1999-2006 (A-D): PAD200 (vigorous) / PAD320 (moderate) over the past 30 days
#   2007-2018 (E-J): PAQ650/PAQ665 (vigorous/moderate recreational),
#                    PAQ605/PAQ620 (vigorous/moderate work)
# A "Yes" to any available item marks the person as physically active. is_yes()
# returns FALSE for the items absent in a given cycle, so this one expression
# harmonizes both eras.
paq_sub <- paq %>% transmute(
  SEQN,
  physical_active =
    is_yes(get_or_na(paq, "PAD200")) | is_yes(get_or_na(paq, "PAD320")) |
    is_yes(get_or_na(paq, "PAQ650")) | is_yes(get_or_na(paq, "PAQ665")) |
    is_yes(get_or_na(paq, "PAQ605")) | is_yes(get_or_na(paq, "PAQ620"))
)

message("Merging tables ...")
df <- demo_sub %>%
  left_join(bmx_sub,   by = "SEQN") %>%
  left_join(bpx_sub,   by = "SEQN") %>%
  left_join(glu_sub,   by = "SEQN") %>%
  left_join(trigly_sub, by = "SEQN") %>%
  left_join(hdl_sub,   by = "SEQN") %>%
  left_join(tchol_sub, by = "SEQN") %>%
  left_join(mcq_sub,   by = "SEQN") %>%
  left_join(diq_sub,   by = "SEQN") %>%
  left_join(smq_sub,   by = "SEQN") %>%
  left_join(bpq_sub,   by = "SEQN") %>%
  left_join(paq_sub,   by = "SEQN")

# FINDRISC questionnaire-derived inputs.
#   family_history_dm: MCQ300C positive -> first-degree tier, else none.
#   prior_high_glucose: self-reported prediabetes (DIQ160) OR measured impaired
#     fasting glucose (100-125 mg/dL). The measured component is available in
#     every cycle, so this input varies across the full cohort even where the
#     self-report item is absent (1999-2004).
df <- df %>% mutate(
  family_history_dm = if_else(fam_dm_yes, "first_degree", "none"),
  prior_high_glucose = prediabetes_hx |
    (!is.na(fasting_glucose) & fasting_glucose >= 100 & fasting_glucose < 126),
  prior_high_glucose = replace_na(prior_high_glucose, FALSE),
  physical_active = replace_na(physical_active, FALSE),
  # Public-use LMF releases only leading-cause codes 1, 2, 10 for the
  # 2015-2016 and 2017-2018 cycles; stroke (5), diabetes (7) and nephritis (9)
  # are folded into 10. Cause-specific (CV, diabetes) analyses are restricted
  # to cycles with full cause coding via this flag; all-cause uses every cycle.
  cause_coded = !(cycle %in% c("2015-2016", "2017-2018"))
)

message(sprintf("After merge: %d rows", nrow(df)))

# Inclusion criteria per OSF pre-registration
message("Applying inclusion criteria ...")
n_start <- nrow(df)

df <- df %>%
  filter(age >= 20, age <= 79) %>%
  filter(!is.na(waist_cm),
         !is.na(fasting_glucose),
         !is.na(sbp), !is.na(dbp),
         !is.na(triglycerides),
         !is.na(hdl)) %>%
  filter(!prior_cvd) %>%
  filter(!prior_t2d)

n_after_inclusion <- nrow(df)
message(sprintf("After inclusion: %d -> %d (dropped %d)",
                n_start, n_after_inclusion, n_start - n_after_inclusion))

# Link mortality
message("Linking mortality ...")
lmf <- readRDS("data/raw/lmf/lmf_combined.rds")

df <- df %>%
  inner_join(lmf, by = "SEQN") %>%
  filter(ELIGSTAT == 1)

n_with_mortality <- nrow(df)
message(sprintf("After mortality linkage (eligibility): %d", n_with_mortality))

# Construct outcome variables.
# All-cause and CV mortality cap at 10 years (120 person-months).
# Diabetes mortality cap extends to 15 years (180 person-months) because the
# narrow 10y diabetes-as-underlying-cause definition yielded only 6 events at
# N=17k, statistically unreliable for Fine-Gray. The diabetes-specific
# definition is broadened in parallel (see event_dm below).
message("Building outcomes ...")
df <- df %>% mutate(
  # 10y horizon for all-cause and CV
  followup_months = pmin(PERMTH_EXM, 120, na.rm = TRUE),
  followup_years  = followup_months / 12,

  # 15y horizon for diabetes only
  followup_months_dm = pmin(PERMTH_EXM, 180, na.rm = TRUE),
  followup_years_dm  = followup_months_dm / 12,

  # All-cause within 10 years
  event_allcause = as.integer(MORTSTAT == 1 & PERMTH_EXM <= 120),
  event_allcause = replace_na(event_allcause, 0L),

  # All-cause within 15 years (for diabetes competing-risk frame)
  event_allcause_dm = as.integer(MORTSTAT == 1 & PERMTH_EXM <= 180),
  event_allcause_dm = replace_na(event_allcause_dm, 0L),

  # CV mortality (NCHS UCOD113 recodes 1 = heart disease, 5 = stroke)
  # The UCOD_LEADING field uses these collapsed codes for cause classification.
  event_cv = as.integer(
    MORTSTAT == 1 &
    UCOD_LEADING %in% c(1, 5) &
    PERMTH_EXM <= 120
  ),
  event_cv = replace_na(event_cv, 0L),

  # Diabetes-related mortality, broadened definition within 15-year follow-up.
  # Triggers on any of:
  #   UCOD_LEADING == 7 (diabetes mellitus, NCHS 113-cause recode)
  #   UCOD_LEADING == 9 (nephritis/nephrotic syndrome/nephrosis, frequent
  #                      diabetic-kidney complication)
  #   DIABETES == 1     (LMF flag for diabetes as a contributing cause,
  #                      even when not the underlying cause)
  event_dm = as.integer(
    MORTSTAT == 1 &
    (UCOD_LEADING %in% c(7, 9) | DIABETES == 1) &
    PERMTH_EXM <= 180
  ),
  event_dm = replace_na(event_dm, 0L),

  # Competing-risks indicators
  # competing_cv: 1=CV death, 2=non-CV death, 0=censored (10y frame)
  competing_cv = case_when(
    event_cv == 1                     ~ 1L,
    event_allcause == 1 & event_cv == 0 ~ 2L,
    TRUE                              ~ 0L
  ),
  # competing_dm: 1=diabetes-related death, 2=non-diabetes death, 0=censored
  # (15y frame, broadened diabetes definition above)
  competing_dm = case_when(
    event_dm == 1                            ~ 1L,
    event_allcause_dm == 1 & event_dm == 0   ~ 2L,
    TRUE                                     ~ 0L
  )
)

dir.create("data/processed", recursive = TRUE, showWarnings = FALSE)
saveRDS(df, "data/processed/analysis_cohort.rds")

message(sprintf("\n=== Final analysis cohort: %d subjects ===", nrow(df)))
message(sprintf("  CV deaths within 10y:                %d (%.2f%%)",
                sum(df$event_cv), 100 * mean(df$event_cv)))
message(sprintf("  Diabetes deaths within 15y (broad):  %d (%.2f%%)",
                sum(df$event_dm), 100 * mean(df$event_dm)))
message(sprintf("  All-cause deaths within 10y:         %d (%.2f%%)",
                sum(df$event_allcause), 100 * mean(df$event_allcause)))
message(sprintf("  All-cause deaths within 15y:         %d (%.2f%%)",
                sum(df$event_allcause_dm), 100 * mean(df$event_allcause_dm)))
message(sprintf("  Median follow-up (10y frame):        %.2f years",
                median(df$followup_years)))
message(sprintf("  Mean follow-up (10y frame):          %.2f years",
                mean(df$followup_years)))
message(sprintf("  Median follow-up (15y frame, dm):    %.2f years",
                median(df$followup_years_dm)))
message(sprintf("  Mean follow-up (15y frame, dm):      %.2f years",
                mean(df$followup_years_dm)))
message("  Note: diabetes outcome uses extended 15y follow-up with broadened cause")
message("        definition (UCOD_LEADING in {7, 9} OR DIABETES contributing-cause flag).")
