# Changelog

All notable changes to this project are listed here. The format follows [Keep a Changelog](https://keepachangelog.com/en/1.1.0/). There are no tagged releases yet.

## [Unreleased]

### Added
- Demo GIF, CONTRIBUTING.md, issue templates and a PR template.

### Changed
- README rewritten.
- The score tests now run on every pull request, not only when score files change.

## 2026-10-08

### Fixed
- `scripts/02` stops when a Linked Mortality File is missing or fails to parse.
- `scripts/10` checks that the cohort RDS exists before reading its timestamp.
- `scripts/15` and the manuscript match the 22-leaf B9 tree.
- The manuscript no longer claims an integrated Brier score. The code writes a point Brier only.
- Script headers for `05`, `13` and `14` match what the code does.

## 2026-06-08

### Changed
- The H3 test and the registered delta-AUCs use the paired IPCW time-dependent delta as the primary metric. The binary DeLong delta is kept as a secondary check.

### Fixed
- Brier scores use each score's recalibrated risk, not the raw score.
- Audit corrections to reproducibility and protocol consistency.

## 2026-05-21

### Added
- Related-work positioning in the manuscript, from a literature scan.

### Fixed
- Second correctness audit: fasting subsample weights (WTSAF), a dual XGBoost ceiling, and FINDRISC cut points.

## 2026-05-20

### Added
- Table 1, calibration curves, subgroup analyses and the B9 US refit.

### Fixed
- Correctness audit and a full re-run of the pipeline.

## 2026-05-19

### Added
- 500-rep PSU-cluster bootstrap CIs for the primary AUCs and delta-AUCs.
- Manuscript in LaTeX: introduction, methods, results, discussion and abstract.
- D2 diagrams of the pipeline, repository layout and scoring system.
- GitHub Pages site in `docs/`.

## 2026-05-18

### Added
- Phase 3: decision curve analysis, pairwise AUC comparisons with IDI and NRI, and the XGBoost ceiling.
- Study protocol draft with Phase 2 and 3 results.

### Changed
- Broadened the diabetes-related death definition and moved the AUC horizons off the follow-up cap.
- Softened the H3 non-inferiority claim to match the CI lower bound.

### Fixed
- `timeROC` with `iid = FALSE` to avoid running out of memory at N = 17k.
- Categorical variable handling and a missing diastolic pressure filter in cohort construction.

## 2026-05-17

### Added
- First version: the five risk scores (RMRS, B9 tree, PCE, Framingham 2008, FINDRISC) with testthat suites.
- Data pipeline scripts for NHANES download, mortality linkage and inclusion criteria.
- Survival analysis scripts for all-cause, cardiovascular and diabetes-related death.
- `renv` and `uv` environments, Makefile, protocol draft and literature scan.

### Fixed
- `renv` library path, Linked Mortality File column positions and HDL column names across cycles.
