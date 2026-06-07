# Study Protocol and Analysis Plan: Longitudinal Validation of MetS Risk Scores

**Status:** Retrospective study protocol and analysis plan, published in the project repository alongside the code and results. This is a transparent written protocol rather than a timestamped prospective registration: it documents the analysis plan, hypotheses, and margins together with the primary results from Phase 2 and Phase 3, which were already computed. It is not a prospective OSF registration.
**Author:** Jay Hemnani (independent researcher; no current institutional affiliation listed)
**Protocol date:** 2026-05-18 (finalized alongside the analysis; see section 0 for corrections)

---

## 0. Post-draft corrections and final analysis (authoritative)

This protocol was drafted with an early set of computed results in sections 3 and 8. The analysis was then corrected on several points; the final numbers are in the manuscript (`manuscript/main.pdf`) and the result CSVs under `results/`, which are authoritative and supersede the original-draft values retained below for transparency. The corrections were:

1. **Reproducible data download.** The download script did not retrieve the 1999 to 2004 fasting-lab tables, which NHANES publishes under legacy names, so a fresh run silently dropped those cycles. This is fixed and verified: a clean checkout now rebuilds the exact N = 17,031 cohort.
2. **Cause-coded subcohort for cause-specific outcomes.** Cardiovascular and diabetes-related analyses use the 1999 to 2014 cause-coded subcohort (N = 13,836; FINDRISC diabetes N = 13,806), not the full N = 17,031 shown in the section 8 tables, because the 2015 to 2018 public-use mortality files collapse cause of death to three categories (verified directly in the linked data).
3. **Full decision-tree reconstruction.** The B9 tree is reconstructed in full (22 leaves) from the paper's Supporting Information S1 Table and applied with its native Korean feature scaling, replacing an earlier three-value reconstruction. B9 discrimination rose (diabetes-related AUC 0.582 to 0.692), and the RMRS-over-B9 deltas fell to +0.017 (all-cause), +0.032 (cardiovascular), and +0.081 (diabetes-related). H5 is restated: RMRS modestly but consistently exceeds B9 rather than dominating it.
4. **Within-stratum bootstrap (primary inference) and correct H3 metric.** The PSU-cluster bootstrap now resamples primary sampling units within each design stratum, the standard NHANES scheme, rather than from the pooled cluster set, which over-dispersed by ignoring the strata. The primary H3 metric is the paired IPCW time-dependent delta-AUC (matching the primary discrimination metric used throughout), not the binary-outcome DeLong delta. Under this bootstrap the primary H3 delta-AUC is -0.020, 95% CI (-0.068, +0.028), whose lower bound falls below the -0.05 margin, so non-inferiority is NOT established on the primary metric. A secondary binary-outcome DeLong delta was +0.013 (-0.043, +0.063) and would have cleared the margin, but it uses a different outcome encoding, ignores the survey design, and is not the basis for the conclusion. FINDRISC remains numerically higher than RMRS on the primary AUC (0.770 vs 0.752); the two scores are best read as statistically indistinguishable rather than one being non-inferior to the other.
5. **Clustered reclassification CIs.** The H4 IDI and NRI intervals come from the within-stratum cluster bootstrap (survey-valid) rather than the estimator's internal i.i.d. perturbation; both still exclude zero (IDI 0.004 [0.001, 0.009]; continuous NRI 0.139 [0.047, 0.267]).
6. **Survey-weighted AUC sensitivity** was added; the design effect on discrimination is small (within about 0.04 AUC).

---

## 1. Title

Longitudinal validation of metabolic syndrome risk scores (RMRS and B9 decision tree) against the ACC/AHA Pooled Cohort Equations, Framingham 2008, and FINDRISC on NHANES 1999-2018 with Linked Mortality File follow-up to 2019.

## 2. Summary

This study externally validates two metabolic syndrome (MetS) risk scoring methods, the Robust MetS Risk Score (RMRS, Shin et al. PeerJ Computer Science 2024) and the B9 decision tree (Shin et al. PLoS One 2023), against the ACC/AHA Pooled Cohort Equations (PCE), Framingham 2008, and FINDRISC, on NHANES 1999-2018 with Linked Mortality File 2019-release follow-up. Three hard outcomes are evaluated: cardiovascular mortality, diabetes-related mortality (broadened definition, see section 6), and all-cause mortality. The analytic cohort is 17,031 adults aged 20-79 with complete fasting biomarkers and mortality linkage eligibility.

The main finding is that the value of MetS-derived scores is outcome-specific. Three results paragraphs follow. The specific numbers in those paragraphs are original-draft values, superseded by the corrected final results in section 0 and the manuscript; in particular, RMRS and FINDRISC are statistically indistinguishable on diabetes-related mortality with non-inferiority not established at the -0.05 margin on the primary IPCW metric, and the full 22-leaf B9 reconstruction is competitive with rather than dominated by RMRS.

**Cardiovascular and all-cause mortality.** The two MetS-derived scores do not match the clinical risk equations on discrimination. For CV mortality at the 9.5-year horizon, Framingham 2008 reaches IPCW time-dependent AUC 0.859 and PCE reaches 0.811, against RMRS at 0.662 and B9 at 0.565. The DeLong tests against horizon-cap binary outcomes give RMRS vs Framingham delta-AUC -0.190 (95% CI -0.231 to -0.148, p = 4.2e-19) and RMRS vs PCE delta-AUC -0.200 (95% CI -0.251 to -0.149, p = 2.4e-14). For all-cause mortality at 9.5 years, Framingham 2008 reaches AUC 0.810 and PCE reaches 0.758, against RMRS at 0.595 and B9 at 0.525. Decision Curve Analysis (DCA) over thresholds 0.05 to 0.30 for CV mortality shows PCE and Framingham generating positive net benefit across the clinically relevant range, with RMRS and B9 collapsing to zero net benefit at thresholds above approximately 0.05. The XGBoost ML upper-bound on the same MetS feature set reaches AUC 0.813 for CV mortality and 0.813 for all-cause mortality at 9.5 years, which does not exceed Framingham. For these two outcomes the clinical risk equations remain the right tool.

**Diabetes-related mortality.** For the broadened diabetes-related mortality outcome at the 14.5-year horizon, RMRS reaches IPCW AUC 0.747 against FINDRISC at 0.766. The DeLong test on the horizon-cap binary outcome gives delta-AUC +0.011 (95% CI -0.055 to +0.077, p = 0.75); the two scores are statistically indistinguishable. Incremental-value analysis of RMRS added on top of FINDRISC gives continuous NRI 0.173 (95% CI 0.033 to 0.287) and IDI 0.0031 (95% CI 0.0002 to 0.0099); the delta-AUC for the joint model over FINDRISC alone is +0.030 (95% CI -0.002 to +0.061, p = 0.068). DCA at thresholds 0.01 to 0.03 (the clinically relevant screening range for diabetes-related death) shows FINDRISC as the top score and RMRS as the second-best, both above "treat all" and "treat none". RMRS has a defined diabetes-outcome lane and contributes information beyond FINDRISC on reclassification metrics.

**B9 decision tree vs RMRS.** The B9 decision tree is consistently dominated by RMRS across all three outcomes. Delta-AUC for RMRS vs B9 is +0.066 (95% CI 0.044 to 0.088, p = 3.6e-9) for all-cause mortality at 9.5 years, +0.086 (95% CI 0.039 to 0.133, p = 3.2e-4) for CV mortality at 9.5 years, and +0.149 (95% CI 0.079 to 0.219, p = 3.2e-5) for diabetes-related mortality at 14.5 years. The continuous geometric-similarity scoring methodology of RMRS outperforms the CART tree formulation derived from the same data.

**Future work.** The XGBoost upper bound on diabetes-related mortality is 0.799 at 14.5 years (5-fold CV mean, SD 0.053). This leaves roughly 5 AUC points of headroom over both RMRS and FINDRISC, motivating future work on a refined tree-based or boosted MetS scoring variant for diabetes-mortality prediction.

The two source papers reported cross-sectional discrimination of the MetS diagnosis label only. This protocol tests longitudinal predictive validity using survey-weighted Fine-Gray competing-risks models and Cox proportional hazards, in line with TRIPOD+AI (Collins et al. BMJ 2024) and STROBE reporting standards.

**Positioning vs prior work.** Park et al. 2025 (JMIR 27:e67525) recently developed and validated a different noninvasive MetS predictive model with CVD risk assessments on multicohort data. This project differs in that it externally validates the Shin/Shim/Oh-specific scores (RMRS and the B9 decision tree) on NHANES + Linked Mortality File, rather than developing a new model. No prior work tests these specific scores against PCE, Framingham, or FINDRISC on US longitudinal mortality outcomes to our knowledge. See `prereg/literature-scan.md` for full context.

## 3. Hypotheses

The hypothesis structure below reflects the locked narrative. The draft results stated per hypothesis are the original-draft values; the corrected final values are summarized in section 0 and reported in full in the manuscript. Primary inference in the final analysis is interval-based (bootstrap confidence intervals against the protocol non-inferiority margin and against zero), so the Bonferroni and Benjamini-Hochberg p-value adjustments named in the original plan are not applied to the primary conclusions, and the missing-data plan is complete-case rather than the multiple imputation named in the original section 7 (see section 0 and the reconciled section 7); any DeLong p-values shown are unadjusted secondary checks.

### Primary hypotheses

- **H1** (RMRS vs Framingham 2008 on CV mortality): RMRS will not match Framingham 2008 on 9.5-year IPCW time-dependent AUC for CV mortality. Operationalized as a two-sided DeLong test of delta-AUC; primary direction is RMRS below Framingham by a clinically meaningful margin. Draft result: delta-AUC = -0.190, 95% CI (-0.231, -0.148), p = 4.2e-19.
- **H2** (RMRS vs PCE on CV mortality): RMRS will not match PCE on 9.5-year IPCW time-dependent AUC for CV mortality. Operationalized as a two-sided DeLong test. Draft result: delta-AUC = -0.200, 95% CI (-0.251, -0.149), p = 2.4e-14.
- **H3** (RMRS matches FINDRISC on diabetes-related mortality): RMRS will match FINDRISC on 14.5-year IPCW time-dependent AUC for diabetes-related mortality (broadened definition; see section 6, Outcomes). Operationalized as a two-sided DeLong test with statistical non-inferiority defined as the 95% CI lower bound of the delta-AUC remaining above -0.05. Draft result: delta-AUC = +0.011, 95% CI (-0.055, +0.077), p = 0.75. The definitive interval is the within-stratum PSU-cluster bootstrap reported in section 0 and the Bootstrap sensitivity results section; on the primary IPCW time-dependent delta-AUC it is -0.020 (95% CI -0.068, +0.028), whose lower bound falls below the -0.05 margin, so non-inferiority is not established on the primary metric.
- **H4** (RMRS adds incremental value to FINDRISC on diabetes-related mortality): RMRS added on top of FINDRISC will improve reclassification on diabetes-related mortality at 14.5 years. Operationalized as continuous NRI and IDI (`survIDINRI::IDI.INF`), with primary success defined as both 95% CIs excluding zero. Draft result: NRI = 0.173, 95% CI (0.033, 0.287); IDI = 0.0031, 95% CI (0.0002, 0.0099). Both metrics meet the primary success criterion.
- **H5** (RMRS dominates B9 on all three outcomes): RMRS will outperform the B9 decision tree on AUC for all three outcomes. Operationalized as DeLong tests on the horizon-cap binary outcome. Draft results: delta-AUC = +0.066 (p = 3.6e-9) for all-cause at 9.5y; +0.086 (p = 3.2e-4) for CV at 9.5y; +0.149 (p = 3.2e-5) for diabetes-related at 14.5y.

### Secondary hypotheses (Benjamini-Hochberg FDR at α=0.05)

- B9 vs Framingham on CV mortality, B9 vs PCE on CV mortality (expected: B9 strongly dominated).
- B9 vs FINDRISC on diabetes-related mortality (expected: B9 strongly dominated; draft delta-AUC -0.137, p = 1.1e-5).
- PCE vs Framingham on CV mortality (expected: indistinguishable; draft delta-AUC +0.009, p = 0.27).
- FINDRISC, PCE, Framingham AUCs on all-cause mortality.
- FINDRISC AUC on CV mortality.

### Exploratory analyses

- Subgroup analyses by sex, race/ethnicity, age band, for all primary outcomes.
- Sensitivity analyses without follow-up cap, complete-cases-only, single-cycle stratification.
- US-refit B9 left-subtree comparison for transportability (Task 5.2b in the project spec).
- XGBoost ML upper-bound baseline on the MetS feature set (bounds achievable discrimination from MetS factors alone). Draft ML upper-bound at 14.5y for diabetes-related mortality is AUC 0.799 (5-fold CV mean, SD 0.053).

## 4. Design plan

- **Study type:** Observational cohort, secondary analysis of public data.
- **Blinding:** Not applicable. Outcomes are objective (death records).
- **Randomization:** Not applicable.
- **Study design:** Retrospective prediction-model external validation.

## 5. Sampling plan

### Data source
- NHANES continuous cycles 1999-2000 through 2017-2018 (10 cycles, approximately 80,000 baseline adults).
- NHANES Linked Mortality File 2019 release (public-use, NCHS probabilistic match).

### Inclusion criteria
- Age 20-79 at NHANES exam.
- Complete fasting biomarkers: triglycerides, HDL cholesterol, fasting glucose (restricts to morning fasting subsample, approximately 50% of adults).
- Complete anthropometric: waist circumference, blood pressure.
- Eligible for mortality follow-up (LMF ELIGSTAT = 1).
- No prior cardiovascular disease at baseline (MCQ160B/C/E/F).
- No prior type-2 diabetes at baseline (DIQ010).

### Actual analytic sample size
**N = 17,031 adults** after all inclusion criteria. The PCE-applicable subsample (age 40-79 with PCE inputs observed) is N = 9,815; Mexican American and Other Hispanic participants are mapped to the white PCE equation rather than excluded, since the published equations provide coefficients only for white and Black. The FINDRISC subsample (with non-missing FINDRISC inputs) is N = 16,996 for all-cause and N = 13,806 within the 1999-2014 cause-coded subcohort used for the diabetes-related outcome.

## 6. Variables

### Outcomes
1. **All-cause mortality**: LMF MORTSTAT = 1 within 10 years of exam (operationalized as t = 9.5 years for IPCW AUC; see section 7). Cox proportional hazards.
2. **CV mortality**: LMF UCOD_LEADING in {1 (heart disease), 5 (cerebrovascular)} within 10 years (t = 9.5 for AUC). Fine-Gray subdistribution hazards with non-CV death as the competing event.
3. **Diabetes-related mortality (broadened definition, 15-year follow-up)**: LMF UCOD_LEADING in {7 (diabetes mellitus), 9 (nephritis / nephrotic syndrome / nephrosis, frequently a diabetic-kidney complication)} OR LMF DIABETES contributing-cause flag = 1, within 15 years of exam (t = 14.5 for AUC). Fine-Gray subdistribution hazards with non-diabetes death as the competing event.

   *Rationale for broadening.* A narrow definition (UCOD_LEADING == 7 only, capped at 10 years) yielded only 6 events at N = 17,031 on the first pipeline run. Six events is statistically unreliable for Fine-Gray subdistribution-hazards estimation (rule of thumb: at least 10 to 20 events per predictor parameter). The pre-specified expansion combines three changes:

   - Add UCOD_LEADING == 9 (nephritis / nephrotic syndrome / nephrosis), which captures diabetic kidney disease deaths that are often coded as the underlying cause when the death certifier attributes the proximate event to renal failure rather than to diabetes itself.
   - Add the DIABETES contributing-cause flag from the Linked Mortality File, which is set when diabetes appears on any line of the death certificate, not only as the underlying cause.
   - Extend follow-up for the diabetes outcome alone to 15 years (180 months), to accumulate sufficient events for the competing-risks model.

   Under the broadened definition the event count is 75 at 14.5 years. All-cause and CV outcomes remain capped at 10 years.

### Risk scores under test (5)
1. **RMRS** (Shin et al. 2024 PeerJ CS): continuous, [0, 1]. Inputs: waist, sBP, dBP, HDL, fasting glucose, triglycerides, sex.
2. **B9 decision tree** (Shin et al. 2023 PLoS One): leaf-probability output, [0, 1]. Inputs: waist, sBP, dBP, sex.
3. **PCE** (Goff et al. 2014, Yadlowsky et al. 2018): 10-year ASCVD risk, [0, 1]. Inputs: age, sex, race, total cholesterol, HDL, sBP, BP treatment, smoking, diabetes status.
4. **Framingham 2008** (D'Agostino et al.): 10-year CVD risk, [0, 1]. Same inputs as PCE minus race.
5. **FINDRISC** (Lindström & Tuomilehto 2003): 0 to 26 point ordinal score. Inputs: age, BMI, waist, physical activity, vegetable intake, BP medication, prior high glucose, family history of diabetes.

### Other variables collected for covariates and stratification
- Sex (male, female).
- Race/ethnicity (NHANES RIDRETH3 or RIDRETH1): Mexican American, Other Hispanic, NH White, NH Black, NH Asian (2011 cycles forward only), Other.
- Age band (20-39, 40-59, 60-79).
- NHANES survey design variables (SDMVPSU, SDMVSTRA, WTMEC2YR, fasting subsample WTSAF2YR).

## 7. Analysis plan

### Outcome models
- **All-cause mortality**: Cox proportional hazards survey-weighted (`survey::svycoxph`), assuming independent censoring.
- **CV mortality** and **diabetes-related mortality**: Fine-Gray subdistribution hazards (`riskRegression::FGR`) with competing-risks framework. Survey weights applied where the software supports it; sensitivity analysis without weights to assess weight impact.

### Discrimination
- Time-dependent ROC AUC using inverse probability of censoring weighting (IPCW) via `timeROC::timeROC`. Horizons:
  - All-cause and CV outcomes: 5 years and a 10-year horizon operationalized as t = 9.5 years.
  - Diabetes outcome (15-year follow-up): 5 years, 10 years, and a 15-year horizon operationalized as t = 14.5 years.
  - The horizons at t = 9.5 and t = 14.5 substitute for the nominal 10-year and 15-year endpoints to avoid IPCW weight degeneracy at the exact follow-up cap, where the censoring distribution has no mass to the right of the evaluation time and weights blow up to NA. Shifting the evaluation point inward by 0.5 years preserves the 10-year and 15-year framing while keeping the IPCW estimator well-defined.
- `iid = FALSE` for the AUC estimator at N approximately 17k due to memory constraints. Bootstrap CIs computed separately (see the Bootstrap sensitivity results section).
- Pairwise DeLong tests via `pROC::roc.test` on the horizon-cap binary outcome (subjects observed past the horizon are coded as failures if the event occurred by then, non-failures otherwise; censored-before-horizon subjects are excluded). This is the operational substitute for `timeROC::compare` because `timeROC` with `iid = TRUE` runs out of memory at N = 17,031. Effective sample sizes after horizon-cap restriction range from n_eval = 5,688 (diabetes 14.5y, FINDRISC subsample) to n_eval = 9,981 (CV and all-cause 9.5y, full cohort).

### Calibration (as run)
- Calibration plots at the primary horizon for each outcome, on each score's recalibrated horizon-specific absolute risk (Fine-Gray CIF for the competing-risks outcomes, Cox complement for all-cause), with observed incidence by Aalen-Johansen or one minus Kaplan-Meier.
- IPCW Brier and integrated Brier scores via `riskRegression::Score`, computed on the recalibrated absolute-risk scale (raw scores such as FINDRISC are not probabilities, so Brier on the raw score is not meaningful and is not reported).
- The original draft also named calibration-in-the-large, a calibration slope, and a Spiegelhalter z-test; these summary statistics were not separately reported (deviation from the original plan).

### Decision analysis (as run)
- Decision Curve Analysis (DCA) at the primary horizon for each outcome, computed directly from the per-score recalibrated risk and the Aalen-Johansen cause-specific incidence within the treat-positive group (`R/utils/dca_competing.R`), rather than via `dcurves::dca`, so that censoring and competing deaths are handled. Single-threshold bootstrap summaries are placed on each outcome's achievable predicted-risk range (5% all-cause, 2% CV, 1% diabetes-related). The recalibration is fit and evaluated on the same cohort, so DCA is apparent in-sample and unweighted (exploratory).
- Compared against "treat all" and "treat none" reference strategies.

### Reclassification metrics
- Integrated Discrimination Improvement (IDI) and Continuous Net Reclassification Index (NRI) for primary pairwise comparisons (`survIDINRI::IDI.INF`). Primary application: FINDRISC + RMRS vs FINDRISC alone on diabetes-related mortality at 14.5 years.

### Follow-up handling
- Primary: cap follow-up at 10 years for CV and all-cause, 15 years for diabetes-related, for consistent prediction-window framing.
- Sensitivity: rerun without cap.

### Missing data (as run: complete-case)
- The final analysis is **complete-case** on the MetS components that define eligibility; scores requiring additional inputs are evaluated on their non-missing subset, with the per-score analytic N reported alongside each result.
- The original draft named multiple imputation by chained equations (`mice`, 20 imputations) with Rubin's-rules pooling. This was **not** used: the cohort is complete on the eligibility-defining components by construction, so multiple imputation was not run. This is a deviation from the original plan (see section 0).

### Multiplicity correction (as run: interval-based)
- The final primary inference is interval-based: bootstrap confidence intervals against the protocol non-inferiority margin and against zero. No p-value multiplicity adjustment is applied to the primary conclusions.
- The original draft named Bonferroni at α=0.01 for the five primary tests and Benjamini-Hochberg FDR at α=0.05 for the secondary set. Because the final inference is interval-based rather than p-value-based, these adjustments are not applied (deviation from the original plan); any DeLong p-values shown are unadjusted secondary checks.
- Exploratory tests: descriptive only, no correction.

### Reporting standards
- TRIPOD+AI checklist (Collins et al. BMJ 2024) completed and included as supplementary material.
- STROBE checklist for cohort study aspects.

## 8. Pre-specified primary results

**These tables are original-draft values (full cohort N for every outcome, pre-correction B9), retained for transparency. They are superseded by the corrected final results summarized in section 0 and reported in full in the manuscript, which restrict the cardiovascular and diabetes-related analyses to the cause-coded subcohort (N = 13,836; FINDRISC diabetes N = 13,806) and use the full 22-leaf B9 reconstruction. The B9 rows below are the largest changes: final B9 AUC is 0.585 (all-cause, 9.5y), 0.631 (cardiovascular, 9.5y), and 0.692 (diabetes-related, 14.5y).**

The primary AUC tables below are computed from the Phase 2 + Phase 3 results in the repository (`results/allcause_summary.csv`, `results/cv_summary.csv`, `results/dm_summary.csv`, `results/pairwise_comparisons.csv`, `results/incremental_dm.csv`, `results/xgboost_summary.csv`). All AUCs are IPCW time-dependent estimates from `timeROC::timeROC` with `iid = FALSE` unless otherwise noted.

### 8.1 All-cause mortality (Cox, 10-year follow-up; AUC evaluated at t = 9.5)

| Score             | N      | AUC 5y  | AUC 9.5y |
|-------------------|--------|---------|----------|
| RMRS              | 17,031 | 0.582   | 0.595    |
| B9 tree           | 17,031 | 0.525   | 0.525    |
| PCE               | 9,815  | 0.745   | 0.758    |
| Framingham 2008   | 17,031 | 0.792   | 0.810    |
| FINDRISC          | 16,996 | 0.658   | 0.682    |
| XGBoost ML bound  | 17,031 | 0.785   | 0.813    |

Events at 9.5y: 774. RMRS vs B9 delta-AUC = +0.066 (95% CI 0.044 to 0.088, DeLong p = 3.6e-9).

### 8.2 Cardiovascular mortality (Fine-Gray, 10-year follow-up; AUC at t = 9.5)

| Score             | N      | AUC 5y  | AUC 9.5y |
|-------------------|--------|---------|----------|
| RMRS              | 17,031 | 0.632   | 0.662    |
| B9 tree           | 17,031 | 0.562   | 0.565    |
| PCE               | 9,815  | 0.806   | 0.811    |
| Framingham 2008   | 17,031 | 0.860   | 0.859    |
| XGBoost ML bound  | 17,031 | 0.791   | 0.813    |

Events at 9.5y: 171. Pairwise DeLong tests on the horizon-cap binary outcome:

| Comparison                            | Delta-AUC (binary) | 95% CI            | p-value   |
|---------------------------------------|--------------------|-------------------|-----------|
| RMRS vs PCE                           | -0.200             | (-0.251, -0.149)  | 2.4e-14   |
| RMRS vs Framingham                    | -0.190             | (-0.231, -0.148)  | 4.2e-19   |
| RMRS vs B9                            | +0.086             | (0.039, 0.133)    | 3.2e-4    |
| PCE vs Framingham                     | +0.009             | (-0.007, 0.025)   | 0.27      |

The XGBoost ML upper bound on the MetS feature set does not exceed Framingham 2008 for CV mortality.

### 8.3 Diabetes-related mortality, broadened definition (Fine-Gray, 15-year follow-up; AUC at t = 14.5)

| Score             | N      | AUC 5y  | AUC 10y | AUC 14.5y |
|-------------------|--------|---------|---------|-----------|
| RMRS              | 17,031 | 0.761   | 0.746   | 0.747     |
| B9 tree           | 17,031 | 0.657   | 0.572   | 0.585     |
| FINDRISC          | 16,996 | 0.747   | 0.714   | 0.766     |
| XGBoost ML bound  | 17,031 | 0.787   | 0.793   | 0.799     |

Events at 14.5y: 75. Pairwise DeLong tests on the horizon-cap binary outcome:

| Comparison                            | Delta-AUC (binary) | 95% CI            | p-value   |
|---------------------------------------|--------------------|-------------------|-----------|
| RMRS vs FINDRISC                      | +0.011             | (-0.055, +0.077)  | 0.75      |
| RMRS vs B9                            | +0.149             | (0.079, 0.219)    | 3.2e-5    |
| B9 vs FINDRISC                        | -0.137             | (-0.198, -0.076)  | 1.1e-5    |

Incremental value of RMRS over FINDRISC:

| Metric         | Estimate | 95% CI               |
|----------------|----------|----------------------|
| Continuous NRI | 0.173    | (0.033, 0.287)       |
| IDI            | 0.0031   | (0.0002, 0.0099)     |
| Delta-AUC      | +0.030   | (-0.002, +0.061), p = 0.068 |

DCA at thresholds 0.01 to 0.03 ranks FINDRISC first and RMRS second, both above "treat all" and "treat none". B9 net benefit is zero across the same range because the tree assigns no subjects to the high-probability leaves at NHANES marginals.

The XGBoost ML upper bound for diabetes-related mortality at 14.5 years is 0.799 (5-fold CV mean, SD 0.053), leaving roughly 5 AUC points of headroom over RMRS and FINDRISC.

## 9. Pre-specified sensitivity analyses (completed)

These analyses were part of the plan and have been run; they do not alter the primary conclusions and are reported in the manuscript alongside the primary tables.

1. **Bootstrap 95% CIs at 500 reps** for primary AUC, delta-AUC, and DCA net benefit, using a within-stratum PSU-cluster bootstrap that respects NHANES survey design (SDMVPSU within SDMVSTRA). These are the definitive intervals (see the Bootstrap sensitivity results section below) and supersede the analytic CIs in section 8.
2. **B9 US-refit sensitivity** (`scripts/15_b9_refit.R`). A CART was refit on a single held-out 30% NHANES split using the three BPWC features to predict ATP III metabolic syndrome, as an exploratory check of whether the transported tree's gap reflects Korean calibration rather than the tree method. This is an internal US comparison on one split, not an external mortality validation.
3. **FINDRISC proxy refinement (done).** The FINDRISC implementation uses NHANES-derived family history (MCQ300C), prediabetes / impaired-fasting-glucose (DIQ160 or measured IFG), and physical activity (PAQ) rather than constant placeholders. The daily fruit-and-vegetable item has no cycle-consistent NHANES counterpart and is held at the no-credit level (a documented residual gap); no DR1TOT dietary-recall proxy is used.
4. **Subgroup analyses (done)** by sex, race/ethnicity, and age band, for all three outcomes (`scripts/13_subgroups.R`); the diabetes-related H3 contrast is in the manuscript and the full grid is in `results/subgroup_auc.csv`.

## 10. Deviations from the original analysis plan

Documented for transparency. The deviations below are operational corrections informed by a dry run of the pipeline on the assembled cohort, before the primary hypothesis tests in section 8 were executed.

- **2026-05-18: Diabetes-mortality outcome broadened.** First pipeline run at N = 17,031 produced only 6 events under the narrow UCOD_LEADING == 7 with 10-year cap definition. The diabetes outcome is now defined as UCOD_LEADING in {7, 9} OR DIABETES contributing-cause flag = 1, with a 15-year follow-up cap. Event count under the broadened definition is 75. All-cause and CV outcomes remain at 10 years. See section 6 for full rationale.
- **2026-05-18: AUC horizons shifted off the follow-up cap.** The 10-year and 15-year time-dependent AUCs were degenerate at the exact cap (IPCW weights blow up to NA) on the dry run. Horizons are operationalized as t = 9.5 (CV, all-cause) and t = 14.5 (diabetes). 5-year horizons are unchanged.
- **2026-05-18: DeLong pairwise inference switched to horizon-cap binary outcome.** `timeROC::compare` requires `iid = TRUE`, which runs out of memory at N = 17,031. Pairwise DeLong tests are computed on the horizon-cap binary outcome via `pROC::roc.test`, with effective sample sizes after horizon-cap restriction ranging from n_eval = 5,688 to n_eval = 9,981 (see section 7, Discrimination). This is a documented operational compromise; bootstrap CIs in section 9 will provide the survey-design-respecting inference.

Any further deviation from this protocol is documented in the project GitHub repository, with reason and impact assessment. Planned analyses are reported as specified even where deviations are also reported.

## Bootstrap sensitivity results

The 500-rep within-stratum PSU-cluster bootstrap (seed 20260519, `scripts/14_bootstrap_cis.R`, resampling primary sampling units within each `sdmvstra`) is the protocol's definitive inference; it supersedes the pooled-cluster bootstrap of the original draft (see section 0, correction 4). The primary delta-AUC for each contrast is the paired IPCW time-dependent delta-AUC. For the H3 contrast (RMRS vs FINDRISC on diabetes-related mortality at t = 14.5), the bootstrap percentile interval is delta-AUC = -0.020 (95% CI -0.068 to +0.028). The lower bound is below the protocol -0.05 non-inferiority margin, so non-inferiority is NOT established on the primary metric; FINDRISC is numerically higher than RMRS on the primary AUC (0.770 vs 0.752), and the two are best read as statistically indistinguishable. A secondary binary-outcome DeLong delta was +0.013 (-0.043, +0.063), reported for comparison only. Other primary (IPCW) delta-AUCs: RMRS vs B9 on all-cause mortality at t = 9.5: +0.014 (0.000, 0.029); RMRS vs B9 on CV mortality at t = 9.5: +0.029 (-0.003, 0.064); RMRS vs B9 on diabetes mortality at t = 14.5: +0.060 (0.015, 0.110); PCE vs Framingham on CV mortality at t = 9.5: +0.009 (-0.001, +0.021). The H4 reclassification metrics from the same bootstrap are IDI = 0.004 (0.001, 0.009) and continuous NRI = 0.139 (0.047, 0.267), both excluding zero. Full distributions are stored at `results/bootstrap_*.rds` and `results/bootstrap_summary.csv`.

## 11. Other

### Software
All primary analyses in R 4.3.3 with `renv` lockfile pinning packages (rms pinned to 6.7-1 for R 4.3 compatibility). ML sensitivity analysis in Python 3.11 with `uv` lockfile. Repository: https://github.com/jayhemnani9910/longitudinal-mets-validation.

### Code availability
Repository public from project initialization, with README and reproducibility instructions.

### Ethics
NHANES public-use files and the public-use NHANES Linked Mortality File are de-identified and exempt from additional IRB review. The manuscript will include this statement.

### Author and affiliation
Jay Hemnani (independent researcher; no current institutional affiliation listed). Possible co-authors per outreach to Shin / Shim / Oh in week 7-8 of the implementation timeline.

---

## Document note

This is a retrospective study protocol, not a timestamped prospective registration. It was finalized after Phase 2 (primary survival analysis) and Phase 3 (DCA, pairwise comparisons, XGBoost ML upper bound) were complete, and it combines the analysis plan, the operational deviations triggered by a dry run on the cohort (section 10), the draft findings (section 8, superseded by section 0), and the sensitivity analyses (section 9). It is published in the project repository alongside the code and results for transparency; the corrected final numbers are authoritative in section 0, the Bootstrap sensitivity results section, and the manuscript.
