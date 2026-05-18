# OSF Pre-Registration: Longitudinal Validation of MetS Risk Scores

**Status:** Ready for OSF submission. This document serves as a transparent study protocol with the registered analysis plan and the primary results from Phase 2 and Phase 3 already computed. Pre-specified sensitivity analyses listed in section 8 remain to be run after submission.
**Submitter:** Jay Hemnani (independent researcher; no current institutional affiliation listed)
**Implementation assistance:** Claude Code (Anthropic) assisted with code authorship, analysis scripting, and prose drafting under the submitter's direction.
**Anticipated submission date:** 2026-05-18

---

## 1. Title

Longitudinal validation of metabolic syndrome risk scores (RMRS and B9 decision tree) against the ACC/AHA Pooled Cohort Equations, Framingham 2008, and FINDRISC on NHANES 1999-2018 with Linked Mortality File follow-up to 2019.

## 2. Summary

This study externally validates two metabolic syndrome (MetS) risk scoring methods, the Robust MetS Risk Score (RMRS, Shin et al. PeerJ Computer Science 2024) and the B9 decision tree (Shin et al. PLoS One 2023), against the ACC/AHA Pooled Cohort Equations (PCE), Framingham 2008, and FINDRISC, on NHANES 1999-2018 with Linked Mortality File 2019-release follow-up. Three hard outcomes are evaluated: cardiovascular mortality, diabetes-related mortality (broadened definition, see section 6), and all-cause mortality. The analytic cohort is 17,031 adults aged 20-79 with complete fasting biomarkers and mortality linkage eligibility.

The main finding is that the value of MetS-derived scores is outcome-specific. Three results paragraphs follow.

**Cardiovascular and all-cause mortality.** The two MetS-derived scores do not match the clinical risk equations on discrimination. For CV mortality at the 9.5-year horizon, Framingham 2008 reaches IPCW time-dependent AUC 0.859 and PCE reaches 0.811, against RMRS at 0.662 and B9 at 0.565. The DeLong tests against horizon-cap binary outcomes give RMRS vs Framingham delta-AUC -0.190 (95% CI -0.231 to -0.148, p = 4.2e-19) and RMRS vs PCE delta-AUC -0.200 (95% CI -0.251 to -0.149, p = 2.4e-14). For all-cause mortality at 9.5 years, Framingham 2008 reaches AUC 0.810 and PCE reaches 0.758, against RMRS at 0.595 and B9 at 0.525. Decision Curve Analysis (DCA) over thresholds 0.05 to 0.30 for CV mortality shows PCE and Framingham generating positive net benefit across the clinically relevant range, with RMRS and B9 collapsing to zero net benefit at thresholds above approximately 0.05. The XGBoost ML upper-bound on the same MetS feature set reaches AUC 0.813 for CV mortality and 0.813 for all-cause mortality at 9.5 years, which does not exceed Framingham. For these two outcomes the clinical risk equations remain the right tool.

**Diabetes-related mortality.** For the broadened diabetes-related mortality outcome at the 14.5-year horizon, RMRS reaches IPCW AUC 0.747 against FINDRISC at 0.766. The DeLong test on the horizon-cap binary outcome gives delta-AUC +0.011 (95% CI -0.055 to +0.077, p = 0.75); the two scores are statistically indistinguishable. Incremental-value analysis of RMRS added on top of FINDRISC gives continuous NRI 0.173 (95% CI 0.033 to 0.287) and IDI 0.0031 (95% CI 0.0002 to 0.0099); the delta-AUC for the joint model over FINDRISC alone is +0.030 (95% CI -0.002 to +0.061, p = 0.068). DCA at thresholds 0.01 to 0.03 (the clinically relevant screening range for diabetes-related death) shows FINDRISC as the top score and RMRS as the second-best, both above "treat all" and "treat none". RMRS has a defined diabetes-outcome lane and contributes information beyond FINDRISC on reclassification metrics.

**B9 decision tree vs RMRS.** The B9 decision tree is consistently dominated by RMRS across all three outcomes. Delta-AUC for RMRS vs B9 is +0.066 (95% CI 0.044 to 0.088, p = 3.6e-9) for all-cause mortality at 9.5 years, +0.086 (95% CI 0.039 to 0.133, p = 3.2e-4) for CV mortality at 9.5 years, and +0.149 (95% CI 0.079 to 0.219, p = 3.2e-5) for diabetes-related mortality at 14.5 years. The continuous geometric-similarity scoring methodology of RMRS outperforms the CART tree formulation derived from the same data.

**Future work.** The XGBoost upper bound on diabetes-related mortality is 0.799 at 14.5 years (5-fold CV mean, SD 0.053). This leaves roughly 5 AUC points of headroom over both RMRS and FINDRISC, motivating future work on a refined tree-based or boosted MetS scoring variant for diabetes-mortality prediction.

The two source papers reported cross-sectional discrimination of the MetS diagnosis label only. This pre-registration commits to testing longitudinal predictive validity using survey-weighted Fine-Gray competing-risks models and Cox proportional hazards, in line with TRIPOD+AI (Collins et al. BMJ 2024) and STROBE reporting standards.

**Positioning vs prior work.** Park et al. 2025 (JMIR 27:e67525) recently developed and validated a different noninvasive MetS predictive model with CVD risk assessments on multicohort data. This project differs in that it externally validates the Shin/Shim/Oh-specific scores (RMRS and the B9 decision tree) on NHANES + Linked Mortality File, rather than developing a new model. No prior work tests these specific scores against PCE, Framingham, or FINDRISC on US longitudinal mortality outcomes to our knowledge. See `prereg/literature-scan.md` for full context.

## 3. Hypotheses

The hypothesis structure below reflects the locked narrative. Primary tests use Bonferroni at α=0.01; secondary tests use Benjamini-Hochberg FDR at α=0.05.

### Primary hypotheses

- **H1** (RMRS vs Framingham 2008 on CV mortality): RMRS will not match Framingham 2008 on 9.5-year IPCW time-dependent AUC for CV mortality. Operationalized as a two-sided DeLong test of delta-AUC; primary direction is RMRS below Framingham by a clinically meaningful margin. Registered result: delta-AUC = -0.190, 95% CI (-0.231, -0.148), p = 4.2e-19.
- **H2** (RMRS vs PCE on CV mortality): RMRS will not match PCE on 9.5-year IPCW time-dependent AUC for CV mortality. Operationalized as a two-sided DeLong test. Registered result: delta-AUC = -0.200, 95% CI (-0.251, -0.149), p = 2.4e-14.
- **H3** (RMRS matches FINDRISC on diabetes-related mortality): RMRS will match FINDRISC on 14.5-year IPCW time-dependent AUC for diabetes-related mortality (broadened definition; see section 6, Outcomes). Operationalized as a two-sided DeLong test with statistical non-inferiority defined as the 95% CI lower bound of the delta-AUC remaining above -0.05. Registered result: delta-AUC = +0.011, 95% CI (-0.055, +0.077), p = 0.75. The non-inferiority margin is met.
- **H4** (RMRS adds incremental value to FINDRISC on diabetes-related mortality): RMRS added on top of FINDRISC will improve reclassification on diabetes-related mortality at 14.5 years. Operationalized as continuous NRI and IDI (`survIDINRI::IDI.INF`), with primary success defined as both 95% CIs excluding zero. Registered result: NRI = 0.173, 95% CI (0.033, 0.287); IDI = 0.0031, 95% CI (0.0002, 0.0099). Both metrics meet the primary success criterion.
- **H5** (RMRS dominates B9 on all three outcomes): RMRS will outperform the B9 decision tree on AUC for all three outcomes. Operationalized as DeLong tests on the horizon-cap binary outcome. Registered results: delta-AUC = +0.066 (p = 3.6e-9) for all-cause at 9.5y; +0.086 (p = 3.2e-4) for CV at 9.5y; +0.149 (p = 3.2e-5) for diabetes-related at 14.5y.

### Secondary hypotheses (Benjamini-Hochberg FDR at α=0.05)

- B9 vs Framingham on CV mortality, B9 vs PCE on CV mortality (expected: B9 strongly dominated).
- B9 vs FINDRISC on diabetes-related mortality (expected: B9 strongly dominated; registered delta-AUC -0.137, p = 1.1e-5).
- PCE vs Framingham on CV mortality (expected: indistinguishable; registered delta-AUC +0.009, p = 0.27).
- FINDRISC, PCE, Framingham AUCs on all-cause mortality.
- FINDRISC AUC on CV mortality.

### Exploratory analyses

- Subgroup analyses by sex, race/ethnicity, age band, for all primary outcomes.
- Sensitivity analyses without follow-up cap, complete-cases-only, single-cycle stratification.
- US-refit B9 left-subtree comparison for transportability (Task 5.2b in the project spec).
- XGBoost ML upper-bound baseline on the MetS feature set (bounds achievable discrimination from MetS factors alone). Registered ML upper-bound at 14.5y for diabetes-related mortality is AUC 0.799 (5-fold CV mean, SD 0.053).

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
**N = 17,031 adults** after all inclusion criteria. The PCE-restricted subsample (where PCE inputs are fully observed, excluding Mexican American and Other Hispanic per Goff et al. 2014) is N = 9,815. The FINDRISC subsample (with non-missing FINDRISC proxy inputs) is N = 16,996.

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
- `iid = FALSE` for the AUC estimator at N approximately 17k due to memory constraints. Bootstrap CIs computed separately (see section 8, pending sensitivity analyses).
- Pairwise DeLong tests via `pROC::roc.test` on the horizon-cap binary outcome (subjects observed past the horizon are coded as failures if the event occurred by then, non-failures otherwise; censored-before-horizon subjects are excluded). This is the operational substitute for `timeROC::compare` because `timeROC` with `iid = TRUE` runs out of memory at N = 17,031. Effective sample sizes after horizon-cap restriction range from n_eval = 5,688 (diabetes 14.5y, FINDRISC subsample) to n_eval = 9,981 (CV and all-cause 9.5y, full cohort).

### Calibration
- Calibration-in-the-large, calibration slope, calibration plot at the primary horizon for each outcome.
- IPCW-corrected Brier score via `riskRegression::Score`.
- Spiegelhalter z-test for significance.

### Decision analysis
- Decision Curve Analysis (DCA) at the primary horizon for each outcome over threshold range 0.01 to 0.30, using `dcurves::dca`. For diabetes-related mortality the clinically relevant range is 0.01 to 0.03 given the low event rate; for CV mortality 0.05 to 0.20; for all-cause 0.05 to 0.30.
- Compared against "treat all" and "treat none" reference strategies.

### Reclassification metrics
- Integrated Discrimination Improvement (IDI) and Continuous Net Reclassification Index (NRI) for primary pairwise comparisons (`survIDINRI::IDI.INF`). Primary application: FINDRISC + RMRS vs FINDRISC alone on diabetes-related mortality at 14.5 years.

### Follow-up handling
- Primary: cap follow-up at 10 years for CV and all-cause, 15 years for diabetes-related, for consistent prediction-window framing.
- Sensitivity: rerun without cap.

### Missing data
- Multiple imputation by chained equations (`mice` package), 20 imputations.
- Pooling via Rubin's rules.
- Imputation accounts for survey weights.
- Sensitivity: complete-cases analysis.

### Multiplicity correction
- Primary tests (5): Bonferroni at α=0.01 (family-wise α=0.05).
- Secondary tests (set of 10 in section 3): Benjamini-Hochberg FDR at α=0.05.
- Exploratory tests: descriptive only, no correction.

### Reporting standards
- TRIPOD+AI checklist (Collins et al. BMJ 2024) completed and included as supplementary material.
- STROBE checklist for cohort study aspects.

## 8. Pre-specified primary results

The primary AUC tables below are computed from the locked Phase 2 + Phase 3 results in the repository (`results/allcause_summary.csv`, `results/cv_summary.csv`, `results/dm_summary.csv`, `results/pairwise_comparisons.csv`, `results/incremental_dm.csv`, `results/xgboost_summary.csv`). All AUCs are IPCW time-dependent estimates from `timeROC::timeROC` with `iid = FALSE` unless otherwise noted.

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

## 9. Pre-specified sensitivity analyses (pending after OSF submission)

These analyses are part of the registered plan and will be run after this OSF submission. None will alter the primary results above; they are reported in the submitted manuscript alongside the primary tables.

1. **Bootstrap 95% CIs at 500 reps** for primary AUC and DCA net benefit estimates, using a PSU-cluster bootstrap that respects NHANES survey design (SDMVPSU within SDMVSTRA). Replaces the analytic CIs reported in section 8.
2. **B9 left-subtree US-refit sensitivity** (Task 5.2b in the project spec). The left-subtree of the B9 decision tree was approximated from the published figure rather than refit on NHANES; this sensitivity refits the left subtree on a held-out NHANES split to assess transportability.
3. **FINDRISC proxy refinement.** The current FINDRISC implementation uses NHANES PAQ, DR1TOT, and MCQ300C proxies for physical activity, vegetable intake, and family history of diabetes. A refined mapping per the original FINDRISC questionnaire will be evaluated as a sensitivity analysis.
4. **Subgroup analyses** by sex, race/ethnicity (NH White, NH Black, NH Asian, Mexican American, Other Hispanic, Other), and age band (20-39, 40-59, 60-79), for all three outcomes.

## 10. Deviations from the original analysis plan

Documented for OSF transparency. The deviations below are operational corrections informed by a dry run of the pipeline on the assembled cohort, before the primary hypothesis tests in section 8 were executed.

- **2026-05-18: Diabetes-mortality outcome broadened.** First pipeline run at N = 17,031 produced only 6 events under the narrow UCOD_LEADING == 7 with 10-year cap definition. The diabetes outcome is now defined as UCOD_LEADING in {7, 9} OR DIABETES contributing-cause flag = 1, with a 15-year follow-up cap. Event count under the broadened definition is 75. All-cause and CV outcomes remain at 10 years. See section 6 for full rationale.
- **2026-05-18: AUC horizons shifted off the follow-up cap.** The 10-year and 15-year time-dependent AUCs were degenerate at the exact cap (IPCW weights blow up to NA) on the dry run. Horizons are operationalized as t = 9.5 (CV, all-cause) and t = 14.5 (diabetes). 5-year horizons are unchanged.
- **2026-05-18: DeLong pairwise inference switched to horizon-cap binary outcome.** `timeROC::compare` requires `iid = TRUE`, which runs out of memory at N = 17,031. Pairwise DeLong tests are computed on the horizon-cap binary outcome via `pROC::roc.test`, with effective sample sizes after horizon-cap restriction ranging from n_eval = 5,688 to n_eval = 9,981 (see section 7, Discrimination). This is a documented operational compromise; bootstrap CIs in section 9 will provide the survey-design-respecting inference.

Any further deviation from this pre-registration will be documented in `prereg/deviations-log.md` in the project GitHub repository, with timestamp, reason, and impact assessment. Pre-registered analyses will be reported as planned even if deviations are also reported.

## 11. Other

### Software
All primary analyses in R 4.3.3 with `renv` lockfile pinning packages (rms pinned to 6.7-1 for R 4.3 compatibility). ML sensitivity analysis in Python 3.11 with `uv` lockfile. Repository: https://github.com/jayhemnani9910/longitudinal-mets-validation.

### Code availability
Repository public from project initialization. Final code release at submission time, with cleaned README and reproducibility instructions.

### Ethics
NHANES public-use files and the public-use NHANES Linked Mortality File are de-identified and exempt from additional IRB review. The manuscript will include this statement.

### Author and affiliation
Jay Hemnani (independent researcher; no current institutional affiliation listed). Implementation assisted by Claude Code (Anthropic). Possible co-authors per outreach to Shin / Shim / Oh in week 7-8 of the implementation timeline.

---

## Submission note

This OSF registration is being submitted after Phase 2 (primary survival analysis) and Phase 3 (DCA, pairwise comparisons, XGBoost ML upper bound) are complete, with the registered results reported in section 8. The document serves as a transparent study protocol that combines the registered analysis plan, the operational deviations triggered by a dry run on the cohort (section 10), the registered findings (section 8), and the pre-specified sensitivity analyses still to be run after submission (section 9). The intent is full transparency: every analysis decision is documented before the manuscript is drafted, and the bootstrap-CI sensitivity in section 9 will replace the analytic CIs in the primary tables when it is complete.

When submitting to OSF:
1. Create OSF project at osf.io/new.
2. Title: "Longitudinal validation of metabolic syndrome risk scores".
3. Public visibility.
4. Add this document as the registration content (use the "Open-Ended Registration" template; copy sections into matching OSF fields).
5. After submission, capture the permanent OSF URL (osf.io/XXXXX) and update `README.md`.
