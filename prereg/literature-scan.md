# Literature scan, 2026-05-17

Pre-implementation scan to check for (a) competing in-press work from Shim's team that would kill this project, and (b) prior NHANES-MetS-mortality work needed for manuscript positioning.

## Risk 1: Are Shin / Shim / Oh already publishing this?

**Decision: PROCEED.** No competing in-press longitudinal validation from Shin/Shim/Oh found as of 2026-05-17.

Searches run:
1. `"Shin" "Shim" "Oh" metabolic syndrome NHANES longitudinal mortality 2025 OR 2026`
2. `"triangular areal similarity" RMRS metabolic syndrome follow-up validation 2025 OR 2026`
3. `author:"Hyunseok Shin" OR author:"Sejong Oh" metabolic syndrome 2025 2026 mortality prediction`

Result: only their existing 2023 (PLoS One B9), 2024 (PeerJ CS B8), and 2021/2022 (BIBM, ICCT) papers surfaced. No NHANES or longitudinal-outcome validation from this team is in press or recently published. Their published work is exclusively cross-sectional MetS-label prediction on Korean cohorts (KoGES) and cross-population comparisons (Korean + US NHANES baseline).

## Risk 2: Are other groups publishing similar work?

**Found one direct competitor publication**:

- **Park J-H, Jeong I, Ko G-J, Jeong S, Lee H. 2025.** "Development of a Predictive Model for Metabolic Syndrome Using Noninvasive Data and its Cardiovascular Disease Risk Assessments: Multicohort Validation Study." *Journal of Medical Internet Research* 27:e67525. PMID 40315452. Korea University College of Medicine + Korea University Guro Hospital.
  - This is a DIFFERENT team (Korea University, not Dankook + SJSU).
  - They built and validated their own MetS noninvasive predictive model with CVD risk assessment, across multicohort data.
  - It does NOT explicitly test RMRS, the B9 decision tree, or the specific Shim/Shin/Oh methodology line.
  - **Implication**: position the current project as the first external longitudinal validation of *the Shim/Shin/Oh-specific scores* (RMRS + B9 tree). Cite Park et al. 2025 in the related-work section as parallel methodology with a different model line.

No other 2025-2026 paper found that specifically validates RMRS or the B9 tree against US clinical scores (PCE / Framingham / FINDRISC) on NHANES.

## Prior NHANES-MetS-mortality literature (for positioning)

These are the studies the manuscript should cite to position the contribution.

### Seminal / often-cited prior work
1. **Ford ES. 2005.** "Risks for all-cause mortality, cardiovascular disease, and diabetes associated with the metabolic syndrome: a summary of the evidence." *Diabetes Care* 28(7):1769-1778. Long-standing reference for MetS-mortality risk.

2. **Stern MP, Williams K, González-Villalpando C, Hunt KJ, Haffner SM. 2004.** "Does the metabolic syndrome improve identification of individuals at risk of type 2 diabetes and/or cardiovascular disease?" *Diabetes Care* 27(11):2676-2681. Early MetS-vs-Framingham comparison.

3. **Wang J, Ruotsalainen S, Moilanen L, Lepistö P, Laakso M, Kuusisto J. 2007.** "The metabolic syndrome predicts cardiovascular mortality." *Diabetes Care* 30(5):1138-1140.

### NHANES-specific MetS-mortality validation literature
4. **Ford ES, Li C, Zhao G. 2010.** "Prevalence and correlates of metabolic syndrome based on a harmonious definition among adults in the US." *Journal of Diabetes* 2(3):180-193. Foundational NHANES-MetS prevalence analysis.

5. **CDC Preventing Chronic Disease 2020.** "The Influence of Metabolic Syndrome in Predicting Mortality Risk Among US Adults: Importance of Metabolic Syndrome Even in Adults With Normal Weight." 17:E70. NHANES III + continuous to 2015. n=36,414.

6. **Frontiers in Immunology 2024.** "Association of inflammatory score with all-cause and cardiovascular mortality in patients with metabolic syndrome: NHANES longitudinal cohort study." PMID 39011047. Uses NHANES + LMF, similar methodology to our planned approach but with an inflammation-augmented score. Useful methodology template.

7. **Liu Y et al. 2024.** "Association between lipid accumulation product (LAP), cardiometabolic index (CMI), and MetS in NHANES 2005-2018." PMC11577631. Validates lipid-derived MetS proxies.

### Recent (2025) direct competition
8. **Park J-H et al. 2025.** JMIR 27:e67525 (as above). Different model, but same general space.

### Methodology references for our analyses
9. **Goff DC et al. 2014.** "2013 ACC/AHA Guideline on the Assessment of Cardiovascular Risk." *Circulation* 129(25 Suppl 2):S49-73. PCE original.

10. **Yadlowsky S et al. 2018.** "Clinical Implications of Revised Pooled Cohort Equations for Estimating ASCVD Risk." *Ann Intern Med* 169(1):20-29. PCE recalibration.

11. **D'Agostino RB et al. 2008.** "General Cardiovascular Risk Profile for Use in Primary Care: The Framingham Heart Study." *Circulation* 117(6):743-753. Framingham 2008.

12. **Lindström J, Tuomilehto J. 2003.** "The Diabetes Risk Score: A practical tool to predict type 2 diabetes risk." *Diabetes Care* 26(3):725-731. FINDRISC.

13. **Fine JP, Gray RJ. 1999.** "A Proportional Hazards Model for the Subdistribution of a Competing Risk." *JASA* 94(446):496-509. Fine-Gray method.

14. **Collins GS et al. 2024.** "TRIPOD+AI statement: updated guidance for reporting clinical prediction models that use regression or machine learning methods." *BMJ* 385:q902.

## Decision

**PROCEED with project.** No fatal competition. The novelty claim ("first external longitudinal validation of Shin/Shim/Oh MetS risk scores on NHANES with TRIPOD+AI + competing-risks methodology") holds. Park et al. 2025 is the closest competitor but uses a different model line.

## Next step

OSF pre-registration (Task 1.10) should explicitly position vs Park et al. 2025 in the "Description" section to head off the "not novel" reviewer concern.
