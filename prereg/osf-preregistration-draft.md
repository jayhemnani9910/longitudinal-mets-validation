# OSF Pre-Registration: Longitudinal Validation of MetS Risk Scores

**Status:** DRAFT — to be submitted before any survival models are fit (after Task 1.10 in plan).
**Submitter:** Jay Hemnani
**Anticipated submission date:** Phase 1 completion (~end of week 3)

---

## 1. Title

Longitudinal validation of metabolic syndrome risk scores (RMRS and decision-tree) against the ACC/AHA Pooled Cohort Equations, Framingham, and FINDRISC on NHANES 1999-2018 with Linked Mortality File follow-up to 2019

## 2. Description

This pre-registration covers an external validation study of two metabolic syndrome (MetS) risk scoring methods developed by Shin, Shim, and Oh (PLoS One 2023; PeerJ Computer Science 2024) on the public NHANES + NHANES Linked Mortality File (LMF, 2019 release). The validation tests whether these MetS-derived scores predict three hard long-term outcomes (cardiovascular mortality, diabetes-related mortality, all-cause mortality) over a 10-year follow-up window, benchmarked against three established clinical scores: ACC/AHA Pooled Cohort Equations (PCE, 2013/2018), Framingham 2008 general CVD risk equation, and FINDRISC.

The original Shin/Shim/Oh papers reported cross-sectional discrimination of the MetS diagnosis label only. This pre-registration commits to testing longitudinal predictive validity using survey-weighted Fine-Gray competing-risks models, in line with TRIPOD+AI (Collins et al. BMJ 2024) and STROBE reporting standards.

**Positioning vs prior work.** Park et al. 2025 (JMIR 27:e67525) recently developed and validated a different noninvasive MetS predictive model with CVD risk assessments on multicohort data. This project differs in that it externally validates the *Shin/Shim/Oh-specific scores* (RMRS and the B9 decision tree) on NHANES + Linked Mortality File, rather than developing a new model. To our knowledge no prior work tests these specific scores against PCE / Framingham / FINDRISC on US longitudinal mortality outcomes. See `prereg/literature-scan.md` for full context.

## 3. Hypotheses

### Primary hypotheses (Bonferroni correction at α=0.01)

- **H1** (RMRS vs PCE on CV mortality): The Robust Metabolic Syndrome Risk Score (RMRS) will achieve a 10-year time-dependent AUC for cardiovascular mortality that is statistically non-inferior to the ACC/AHA Pooled Cohort Equations (PCE), defined as the lower 95% CI of the AUC difference being above -0.05.
- **H2** (RMRS vs Framingham on CV mortality): RMRS will achieve a 10-year time-dependent AUC for cardiovascular mortality that is statistically non-inferior to Framingham 2008 by the same non-inferiority margin.
- **H3** (B9 tree vs FINDRISC on diabetes-related mortality): The B9 decision tree (leaf-probability output) will achieve a 10-year time-dependent AUC for diabetes-related mortality that is statistically non-inferior to FINDRISC.
- **H4** (RMRS for all-cause mortality): RMRS will achieve a 10-year time-dependent AUC > 0.60 for all-cause mortality (null: AUC = 0.50).
- **H5** (B9 tree for all-cause mortality): The B9 decision tree will achieve a 10-year time-dependent AUC > 0.60 for all-cause mortality.

### Secondary hypotheses (Benjamini-Hochberg FDR within this set, α=0.05)

The remaining 10 of the 15 score × outcome pairs:

- B9 tree vs PCE on CV mortality
- B9 tree vs Framingham on CV mortality
- RMRS on diabetes-related mortality (no clinical-baseline comparator other than B9)
- RMRS vs FINDRISC on diabetes-related mortality
- B9 tree vs RMRS on CV mortality
- B9 tree vs RMRS on diabetes-related mortality
- PCE on all-cause mortality
- Framingham on all-cause mortality
- FINDRISC on all-cause mortality
- FINDRISC on CV mortality

### Exploratory analyses

- Subgroup analyses (sex, race/ethnicity, age band) for all primary outcomes
- Sensitivity analyses: without follow-up cap, complete-cases-only, single-cycle stratification
- US-refit B9 tree comparison (transportability)
- XGBoost ML upper-bound baseline (bounds achievable discrimination from MetS factors alone)

## 4. Design plan

- **Study type:** Observational cohort, secondary analysis of public data.
- **Blinding:** Not applicable. Outcomes are objective (death records).
- **Randomization:** Not applicable.
- **Study design:** Retrospective prediction-model external validation.

## 5. Sampling plan

### Data source
- NHANES continuous cycles 1999-2000 through 2017-2018 (10 cycles, ~80,000 baseline adults)
- NHANES Linked Mortality File 2019 release (public-use, NCHS probabilistic match)

### Inclusion criteria
- Age 20-79 at NHANES exam
- Complete fasting biomarkers: triglycerides, HDL cholesterol, fasting glucose (restricts to morning fasting subsample, ~50% of adults)
- Complete anthropometric: waist circumference, blood pressure
- Eligible for mortality follow-up (LMF ELIGSTAT = 1)
- No prior cardiovascular disease at baseline (MCQ160B/C/E/F)
- No prior type-2 diabetes at baseline (DIQ010)

### Expected sample size
- After all inclusion criteria: 20,000-30,000 adults (estimate based on prior NHANES analyses)
- Actual N to be filled in this section after Task 1.3 of the implementation plan
- **Actual N: [TO BE FILLED before pre-registration submission]**

## 6. Variables

### Outcomes
1. **All-cause mortality**: LMF MORTSTAT = 1 within 10 years of exam.
2. **CV mortality**: NCHS UCOD113 recodes 54-68 (heart disease) and 70 (cerebrovascular). Competing risk: non-CV death.
3. **Diabetes-related mortality**: NCHS UCOD113 recode 46. Competing risk: non-diabetes death.

### Risk scores under test (5)
1. **RMRS** (Shin et al. 2024 PeerJ CS): continuous, [0, 1]. Inputs: waist, sBP, dBP, HDL, fasting glucose, triglycerides, sex.
2. **B9 decision tree** (Shin et al. 2023 PLoS One): leaf-probability output, [0, 1]. Inputs: waist, sBP, dBP, sex.
3. **PCE** (Goff et al. 2014, Yadlowsky et al. 2018): 10-year ASCVD risk, [0, 1]. Inputs: age, sex, race, total cholesterol, HDL, sBP, BP treatment, smoking, diabetes status.
4. **Framingham 2008** (D'Agostino et al.): 10-year CVD risk, [0, 1]. Same inputs as PCE minus race.
5. **FINDRISC** (Lindström & Tuomilehto 2003): 0-26 point ordinal score. Inputs: age, BMI, waist, physical activity, vegetable intake, BP medication, prior high glucose, family history of diabetes.

### Other variables collected for covariates and stratification
- Sex (male, female)
- Race/ethnicity (NHANES RIDRETH3 or RIDRETH1): Mexican American, Other Hispanic, NH White, NH Black, NH Asian (2011 cycles forward only), Other
- Age band (20-39, 40-59, 60-79)
- NHANES survey design variables (SDMVPSU, SDMVSTRA, WTMEC2YR, fasting subsample WTSAF2YR)

## 7. Analysis plan

### Outcome models
- **All-cause mortality**: Cox proportional hazards survey-weighted (`survey::svycoxph`), assuming independent censoring.
- **CV mortality** and **diabetes-related mortality**: Fine-Gray subdistribution hazards (`riskRegression::FGR`) with competing-risks framework. Survey weights applied where the software supports it; sensitivity analysis without weights to assess weight impact.

### Discrimination
- Time-dependent ROC AUC at 5 years and 10 years using inverse probability of censoring weighting (IPCW), `timeROC::timeROC` with `iid = TRUE`.
- 500 bootstrap iterations for confidence intervals, accounting for NHANES PSU structure (cluster bootstrap).
- Pairwise DeLong tests via `timeROC::compare` for primary comparisons.

### Calibration
- Calibration-in-the-large, calibration slope, calibration plot at 10 years.
- IPCW-corrected Brier score via `riskRegression::Score`.
- Spiegelhalter z-test for significance.

### Decision analysis
- Decision Curve Analysis (DCA) at 10 years over threshold range 0.01-0.30, using `dcurves::dca`.
- Compared against "treat all" and "treat none" reference strategies.

### Reclassification metrics
- Integrated Discrimination Improvement (IDI) and Continuous Net Reclassification Index (NRI) for primary pairwise comparisons (`survIDINRI::IDI.INF`).

### Follow-up handling
- Primary: cap follow-up at 10 years for consistent prediction-window framing across cohorts.
- Sensitivity: rerun without cap.

### Missing data
- Multiple imputation by chained equations (`mice` package), 20 imputations.
- Pooling via Rubin's rules.
- Imputation accounts for survey weights.
- Sensitivity: complete-cases analysis.

### Multiplicity correction
- Primary tests (5): Bonferroni at α=0.01 (family-wise α=0.05).
- Secondary tests (10): Benjamini-Hochberg FDR at α=0.05.
- Exploratory tests: descriptive only, no correction.

### Reporting standards
- TRIPOD+AI checklist (Collins et al. BMJ 2024) completed and included as supplementary material.
- STROBE checklist for cohort study aspects.

## 8. Other

### Deviations log
Any deviation from this pre-registration will be documented in `prereg/deviations-log.md` in the project GitHub repository, with timestamp, reason, and impact assessment. Pre-registered analyses will be reported as planned even if deviations are also reported.

### Software
All analyses in R 4.3.3+ with `renv` lockfile pinning packages. ML sensitivity analysis in Python 3.11 with `uv` lockfile. Repository: https://github.com/jayhemnani9910/longitudinal-mets-validation

### Code availability
Repository public from project initialization. Final code release at submission time, with cleaned README and reproducibility instructions.

### Ethics
NHANES public-use files and the public-use NHANES Linked Mortality File are de-identified and exempt from additional IRB review. Manuscript will include this statement.

### Author
Jay Hemnani (independent researcher). Possible co-authors per outreach to Shin / Shim / Oh in week 7-8 of the implementation timeline.

---

## Submission notes

When submitting to OSF:
1. Create OSF project at osf.io/new
2. Title: "Longitudinal validation of metabolic syndrome risk scores"
3. Public visibility
4. Add this document as the registration content (use "Open-Ended Registration" template; copy sections into matching OSF fields)
5. After submission, capture the permanent OSF URL (osf.io/XXXXX) and update `README.md`

Do **not** submit until Task 1.3 is complete and the actual analysis-cohort N is filled in section 5.
