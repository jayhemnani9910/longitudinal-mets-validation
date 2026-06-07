# Reporting checklists: TRIPOD+AI and STROBE

Supplement to "External longitudinal validation of metabolic syndrome risk scores for cardiovascular, diabetes, and all-cause mortality in the US NHANES Linked Mortality File." Each item lists where it is addressed in the manuscript (Methods, Results, Discussion, Tables, the Data and code availability section, or this protocol) or is marked not applicable with a reason.

This is an external evaluation (validation) of existing published prediction models, so TRIPOD+AI items specific to model development are marked accordingly.

## TRIPOD+AI (2024), 27 items

| # | Item | Where addressed |
|---|------|-----------------|
| 1 | Title identifies model evaluation, population, outcome | Title |
| 2 | Abstract (structured) | Abstract |
| 3a | Healthcare context and rationale | Introduction, paragraphs 1 and 4 |
| 3b | Target population, intended use, intended users | Introduction; the scores are screening/risk tools for adults |
| 3c | Known health inequalities across sociodemographic groups | Subgroup analyses by sex, race/ethnicity, age band (Results; Table on diabetes subgroups) |
| 4 | Objectives (this study evaluates, does not develop) | Introduction, final paragraph |
| 5a | Data sources, separately for development and evaluation | Methods, Study design and data sources; development cohorts are the original Korean studies (cited) |
| 5b | Dates of data collection | Methods (1999 to 2018 NHANES; 2019 LMF release) |
| 6a | Study setting | Methods (US civilian non-institutionalized population) |
| 6b | Eligibility criteria | Methods, Participants and inclusion criteria |
| 6c | Treatments received | Not a treatment study; baseline medication use (antihypertensive) is a model input, noted in Methods |
| 7 | Data pre-processing and quality checks | Methods, Participants (HDL/race coalescing, fasting-lab retrieval); cause-coding verified directly |
| 8a | Outcome definition and time horizon | Methods, Outcomes |
| 8b/8c | Outcome assessment, blinding | Outcomes are objective death-record linkage (Methods); blinding not applicable |
| 9a/9b/9c | Predictor choice, definitions, assessment | Methods, Risk scores (each score's inputs and source) |
| 10 | How study size was determined | Methods (all eligible NHANES fasting-subsample participants); event counts in Results |
| 11 | Missing-data handling | Methods, Statistical analysis (complete-case on the defining components; per-score subsets) |
| 12a-12g | Analysis details | Methods, Statistical analysis (Fine-Gray, Cox, IPCW AUC, recalibration, bootstrap) |
| 13 | Class-imbalance methods | Low event rates handled via competing-risks models and IPCW; the diabetes outcome was broadened to reach adequate events (Methods, Outcomes) |
| 14 | Fairness approaches | Subgroup discrimination by sex, race/ethnicity, age (Results); limitations note NH Asian identifiability |
| 15 | Model output and thresholds | Methods, Risk scores; decision thresholds in the net-benefit analysis |
| 16 | Differences between development and evaluation data | Discussion (Korean development cohorts vs US evaluation; transportability of the tree's calibration) |
| 17 | Ethics committee and consent | Methods, Study design (public de-identified data, exempt) |
| 18a | Funding and funder role | Conflicts of interest / no specific funding (title-page declarations) |
| 18b | Conflicts of interest | Conflicts of interest section |
| 18c | Protocol access | OSF pre-registration (Data and code availability) |
| 18d | Registration | OSF pre-registration (Data and code availability) |
| 18e | Data availability | Data and code availability (public NHANES and LMF) |
| 18f | Analytical code availability | Data and code availability (public repository, pinned environment, reproduces end to end) |
| 19 | Patient and public involvement | None; secondary analysis of existing public data |
| 20a | Participant flow | Results, Cohort (counts at each stage) |
| 20b | Characteristics | Table 1 |
| 20c | Comparison with development data | Discussion (transportability) |
| 21 | Participants and events per analysis | Results, Cohort and table captions (per-outcome N and events) |
| 22 | Full model details so others can apply it | Methods, Risk scores; code repository implements each score exactly |
| 23a | Performance with CIs | Tables (AUC, weighted AUC, delta-AUC, net benefit) with bootstrap CIs |
| 23b | Heterogeneity | Subgroup analyses (Results) |
| 24 | Model-updating results | Decision-tree refit sensitivity (Results, Table on the refit) |
| 25 | Overall interpretation, including fairness | Discussion |
| 26 | Limitations | Discussion, limitations paragraph |
| 27a-27c | Handling of poor/unavailable inputs; user interaction; next steps | Discussion (proxy inputs; future joint models and non-US validation) |

## STROBE (cohort), 22 items

| # | Item | Where addressed |
|---|------|-----------------|
| 1 | Design in title/abstract | Title and Abstract |
| 2 | Background and rationale | Introduction |
| 3 | Objectives with prespecified hypotheses | Introduction (H3 margin); pre-registration |
| 4 | Study design | Methods, Study design |
| 5 | Setting, locations, dates, follow-up | Methods, Study design and Outcomes |
| 6 | Participants, eligibility, follow-up | Methods, Participants; Results, Cohort |
| 7 | Variables (outcomes, predictors, confounders) | Methods, Risk scores and Outcomes |
| 8 | Data sources and measurement | Methods, Study design and Risk scores |
| 9 | Bias | Discussion, limitations (fasting-subsample selection, proxy inputs, unweighted discrimination) |
| 10 | Study size | Results, Cohort; event counts |
| 11 | Quantitative variables | Methods, Risk scores (scaling and thresholds) |
| 12 | Statistical methods | Methods, Statistical analysis |
| 13 | Participants at each stage | Results, Cohort |
| 14 | Descriptive data | Table 1 |
| 15 | Outcome events | Results, Cohort (events per outcome) |
| 16 | Main results with CIs | Tables (AUC, delta-AUC) |
| 17 | Other analyses (subgroups, sensitivity) | Results, Sensitivity analyses; weighted AUC; refit; subgroups |
| 18 | Key results | Discussion, opening paragraph |
| 19 | Limitations | Discussion, limitations paragraph |
| 20 | Interpretation | Discussion |
| 21 | Generalisability | Discussion (US representative sample; transportability of Korean-developed scores) |
| 22 | Funding | Declarations (no specific funding) |
