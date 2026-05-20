# Longitudinal validation of metabolic syndrome risk scores

External longitudinal validation of two metabolic syndrome risk scoring methods (RMRS from Shin et al. 2024 PeerJ CS; decision tree from Shin et al. 2023 PLoS ONE) on the US NHANES Linked Mortality File. Benchmarked against the ACC/AHA Pooled Cohort Equations, Framingham 2008, and FINDRISC for cardiovascular, diabetes-related, and all-cause mortality.

[![License: MIT](https://img.shields.io/badge/License-MIT-blue.svg)](LICENSE)
[![Manuscript: draft](https://img.shields.io/badge/Manuscript-draft-orange.svg)](manuscript/main.pdf)
[![Pre-registration: drafted](https://img.shields.io/badge/Pre--registration-drafted-yellow.svg)](prereg/osf-preregistration-draft.md)
[![Language: R](https://img.shields.io/badge/Language-R%204.3.3-276DC3.svg)](https://www.r-project.org/)

## Headline finding

Outcome-specific value.

- For **cardiovascular and all-cause mortality**, the established clinical risk equations dominate. Framingham 2008 reaches AUC 0.858 at 9.5 years for cardiovascular mortality and 0.810 for all-cause; the metabolic-syndrome-derived RMRS lands at 0.660 and 0.600 respectively.
- For **diabetes-related mortality**, RMRS competes head-to-head with FINDRISC at 14.5 years (AUC 0.752 vs 0.770; delta-AUC +0.014 favoring RMRS, 95% bootstrap CI -0.051 to +0.090) and gives modest incremental value when stacked on FINDRISC (continuous NRI 0.139, IDI 0.0042, both CIs exclude zero). The FINDRISC comparator uses NHANES-derived family history, prediabetes, and physical-activity items rather than placeholders.
- The **B9 decision tree** is consistently outperformed by the continuous RMRS on every outcome (delta-AUC +0.071 all-cause, +0.093 cardiovascular, +0.158 diabetes-related; all CIs exclude zero). An XGBoost ceiling on the same five MetS inputs exceeds both, indicating the underperformance is a tree-structure problem rather than a MetS-signal problem.

## At a glance

| Item | Value |
|------|-------|
| Data source | NHANES 1999 to 2018 + Linked Mortality File 2019 release |
| Cohort N | 17,031 adults aged 20 to 79 years (fasting subsample); 13,836 in the 1999 to 2014 cause-coded subcohort for CV and diabetes outcomes |
| Outcomes | all-cause (n=806, full cohort), cardiovascular (n=173), diabetes-related (n=77, broadened) mortality |
| Scores tested | RMRS, B9 tree, ACC/AHA PCE, Framingham 2008, FINDRISC |
| Inference | survey-weighted Fine-Gray competing risks + Cox; IPCW time-dependent AUC; competing-risks DCA with Fine-Gray CIF recalibration; 500-rep PSU-cluster bootstrap |
| Target venue | JAMA Network Open / npj Digital Medicine |

Full results in `results/*_summary.csv`. Rendered manuscript at [`manuscript/main.pdf`](manuscript/main.pdf). Pre-registered analysis plan at [`prereg/osf-preregistration-draft.md`](prereg/osf-preregistration-draft.md).

## Repository structure

```
.
├── scripts/        Numbered pipeline 01..14 (R + Python)
├── R/scores/       5 risk score implementations + tests
├── R/utils/        NHANES loader and survey-design helpers
├── tests/scores/   testthat unit tests
├── data/           Raw NHANES + LMF (gitignored, regeneratable via Makefile)
├── results/        Output tables, plots, cached models
├── manuscript/     LaTeX source + rendered PDF
├── prereg/         OSF pre-registration draft + literature scan
├── diagrams/       D2 source + SVG renders of the pipeline, structure, and scoring system
└── docs/           Static landing page for GitHub Pages
```

Visual overview of the pipeline, repository structure, and scoring system: see [`diagrams/`](diagrams/).

## Reproducing

```bash
make data         # download NHANES + LMF and assemble the cohort (10-15 min)
make scores       # compute the 5 risk scores per subject (1 min)
make analyses     # all-cause Cox + CV/diabetes Fine-Gray (10 min)
make dca          # decision curve analysis + pairwise comparisons + XGBoost (5 min)
make bootstrap    # 500-rep PSU-cluster bootstrap CIs (about 100 min)
```

Or `make all` for the full pipeline end to end.

The R environment is managed via `renv` with `rms` pinned to 6.7-1 (CRAN archive) for compatibility with R 4.3.3. Python uses `uv` with the lockfile checked in. System dependencies (libxml2, libcurl, libssl, libgsl, fontconfig, freetype, harfbuzz, fribidi, libpng, libtiff, libjpeg) are listed in `DESCRIPTION`.

## Pre-registration

The full analysis plan, including hypothesis structure (H1 through H5), inclusion criteria, score implementations, outcome definitions, IPCW time-dependent AUC horizons, decision curve analysis thresholds, and pre-specified sensitivity analyses, is documented in [`prereg/osf-preregistration-draft.md`](prereg/osf-preregistration-draft.md). OSF submission with timestamped registration is pending.

## Citation

If you build on this work, please cite both the underlying scoring papers and this validation.

The two metabolic syndrome scoring methods evaluated here are from:

- Shin HJ, Oh J, Shim SS. A robust metabolic syndrome risk score using triangular areal similarity. *PeerJ Computer Science*, 2024.
- Shin HJ, Oh J, Shim SS. Machine learning predictive model for metabolic syndrome prevention. *PLoS ONE*, 2023.

This validation work:

- Hemnani J. External longitudinal validation of metabolic syndrome risk scores for cardiovascular, diabetes, and all-cause mortality in the US NHANES Linked Mortality File. Manuscript in preparation, 2026.

A `CITATION.cff` file is provided for GitHub's "Cite this repository" button.

## License

MIT for code (see [`LICENSE`](LICENSE)). CC-BY for figures and manuscript text.
