# Pipeline Pause Status — 2026-05-17

## Where we stopped

Completed `make data` and `make scores`. Partway through `make analyses` (scripts 05–07). Script 05 (`05_survival_allcause.R`) ran and wrote `results/cache/allcause_survival.rds` and `results/allcause_summary.csv` with 5-year AUC values. Scripts 06 and 07 have fixes committed but were not fully executed before the pause.

## What is done

- **renv library**: Fixed broken symlinks for nhanesA, plyr, rvest, xml2, and all other packages that were either at wrong symlink depth or missing from the renv lib directory. All 5 key packages (nhanesA, rms, riskRegression, timeROC, purrr) confirmed loading.
- **Makefile**: All `Rscript` calls updated to `Rscript --no-init-file` to bypass `.Rprofile`/renv auto-activation.
- **Script 01**: Downloaded 97/110 NHANES tables. Early cycles A/B/C for GLU/TRIGLY/HDL/TCHOL retrieved from their period-specific table names (LAB10AM, LAB13AM, LAB13, L10AM_B, etc.). DEMO_I downloaded successfully on retry.
- **Script 02**: LMF column positions corrected (actual file is 48 chars, not 95; PERMTH_EXM at cols 46-48). UCOD_LEADING coerced to integer. All 10 cycles parsed.
- **Script 03**: HDL column harmonization (LBDHDL/LBXHDD/LBDHDD across cycles). Factor label handling for RIAGENDR, RIDRETH1/3, MCQ, DIQ, BPQ, SMQ columns (nhanesA returns text labels after stripping; added is_yes() helper and explicit case_when mapping). Cohort N=17,031 with correct exclusions (prior CVD, prior T2D, missing labs).
- **Script 04**: All 5 scores computed on 17,031 subjects. RMRS, B9, Framingham, FINDRISC complete; PCE on 9,815 (age 40-79 with total cholesterol).
- **Script 05**: All-cause 5-year AUC: RMRS=0.582, B9=0.525, PCE=0.745, Framingham=0.792, FINDRISC=0.658.

## What is blocked / needs attention next run

1. **AUC at 10 years is NA for all outcomes**. Root cause: `PERMTH_EXM` is capped at 120 months, so at t=10 years, `survProb = 0` which makes IPCW weights degenerate. The `timeROC` package cannot compute AUC at the boundary of follow-up. Options: (a) use t=9.5 years as the "10-year" horizon, or (b) use a complementary dataset not capped at 10 years for that estimate, or (c) use `pROC` / `survivalROC` with a different weighting scheme.

2. **Diabetes mortality outcome has only 6 events** within 10 years (out of 17k). Fine-Gray models will likely fail to converge or produce unreliable estimates. The pre-registration should note this as an exploratory analysis; AUC estimates will have very wide uncertainty.

3. **Scripts 06 and 07** have all bugs fixed (FGR formula via as.formula, AUC field names AUC_1[1]/AUC_1[2]) but were not run to completion before pause. Run `make analyses` to complete.

4. **Scripts 08-10** (DCA, pairwise comparisons, XGBoost) have not been touched. The XGBoost script needs feather export (check if arrow is installed or use pyreadstat fallback).

## Cohort summary

| Metric | Value |
|--------|-------|
| Final N | 17,031 |
| All-cause deaths (10y) | 806 (4.73%) |
| CV deaths (10y) | 181 (1.06%) |
| Diabetes deaths (10y) | 6 (0.04%) |
| Median follow-up | 10.00 years |
| Mean follow-up | 7.82 years |

## Next run entry point

```bash
cd /home/po/projects/work/longitudinal-mets-validation
make analyses   # scripts 05 (already done), 06, 07, 07b
make dca        # scripts 08, 09, 10
```

Check `results/allcause_summary.csv` for the completed 5-year AUC table before re-running 05.
