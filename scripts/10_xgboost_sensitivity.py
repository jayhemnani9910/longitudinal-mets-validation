#!/usr/bin/env python3
"""scripts/10_xgboost_sensitivity.py

XGBoost ML upper-bound baseline as a sensitivity analysis. Bounds the
achievable AUC and provides a ceiling against which the clinical risk scores
can be interpreted. Two feature sets are reported:
  - mets5:     the five metabolic-syndrome components only (waist, systolic and
               diastolic BP, HDL, fasting glucose, triglycerides). This is the
               ceiling that speaks to whether the MetS signal itself limits the
               clinical scores.
  - clinical8: the same five components plus age and sex, an upper reference
               that includes the demographic predictors the PCE and Framingham
               equations also use.

Inputs:  data/processed/cohort_with_scores.rds (canonical)
         data/processed/cohort_with_scores.csv (generated from RDS via Rscript
         when missing, since the arrow R package is intentionally not installed)
Outputs: results/xgboost_{mets5,clinical8}_{label}.txt, results/xgboost_summary.csv

NHANES PSU-grouped CV respects the sampling design and avoids optimism from
within-PSU resampling. Median imputation is fit on each training fold only.
Horizons match the survival scripts: censor at 9.5y for allcause / cv
(followup_years), 14.5y for dm (followup_years_dm), with events counted only
inside the horizon.

Caveats (also written into each output file): the AUC is unweighted, so it
describes the sample rather than the weighted US population; and the label is a
binary horizon-capped event (a subject censored before the horizon is treated
as event-free), so these values are an approximation to the IPCW time-dependent
AUC reported by the survival scripts, not an identical statistic.
"""

import os
import subprocess
import sys

import numpy as np
import pandas as pd
import xgboost as xgb
from sklearn.metrics import roc_auc_score
from sklearn.model_selection import GroupKFold


CSV_PATH = "data/processed/cohort_with_scores.csv"
RDS_PATH = "data/processed/cohort_with_scores.rds"


def ensure_csv():
    """Generate the CSV intermediate from the canonical RDS if missing.

    The R-side renv intentionally omits the arrow package (the Apache Arrow
    C++ compile is too slow on this machine), so we round-trip through CSV
    rather than feather. Rscript --no-init-file matches the Makefile.
    """
    if os.path.exists(CSV_PATH) and os.path.getmtime(CSV_PATH) >= os.path.getmtime(RDS_PATH):
        return
    if not os.path.exists(RDS_PATH):
        print(f"Missing {RDS_PATH}. Run `make scores` first.", file=sys.stderr)
        sys.exit(1)
    print(f"Generating {CSV_PATH} from {RDS_PATH} via Rscript ...")
    subprocess.run(
        [
            "Rscript", "--no-init-file", "-e",
            (
                '.libPaths("/home/po/projects/work/longitudinal-mets-validation/'
                'renv/library/R-4.3/x86_64-pc-linux-gnu"); '
                'df <- readRDS("data/processed/cohort_with_scores.rds"); '
                'write.csv(df, "data/processed/cohort_with_scores.csv", row.names = FALSE)'
            ),
        ],
        check=True,
    )


def horizon_labels(df):
    """Return per-outcome (event vector, follow-up vector, horizon, label).

    Horizon-capped events: any event past the horizon is treated as censored,
    follow-up is min(followup, horizon). This matches the t* evaluation
    convention used by the Phase 2 survival AUCs.
    """
    # cause_coded flag: cause-specific (cv, dm) outcomes are restricted to the
    # 1999-2014 cycles with full leading-cause coding, matching the survival
    # scripts. All-cause uses every cycle.
    return [
        ("event_allcause", "followup_years", 9.5, "allcause_9_5y", False),
        ("event_allcause", "followup_years", 5.0, "allcause_5y", False),
        ("event_cv",       "followup_years", 9.5, "cv_9_5y", True),
        ("event_cv",       "followup_years", 5.0, "cv_5y", True),
        ("event_dm",       "followup_years_dm", 14.5, "dm_14_5y", True),
        ("event_dm",       "followup_years_dm", 10.0, "dm_10y", True),
        ("event_dm",       "followup_years_dm", 5.0,  "dm_5y", True),
    ]


# The five metabolic-syndrome components (BP enters as systolic and diastolic).
METS5 = ["waist_cm", "sbp", "dbp", "hdl", "fasting_glucose", "triglycerides"]
# Same components plus the demographic predictors PCE and Framingham also use.
CLINICAL8 = METS5 + ["age", "sex_int"]
FEATURE_SETS = (("mets5", METS5), ("clinical8", CLINICAL8))

CAVEAT_UNWEIGHTED = "caveat=unweighted_AUC_describes_sample_not_weighted_population"
CAVEAT_BINARY = "caveat=binary_horizon_capped_label_approximates_IPCW_not_identical"


def evaluate(df, feature_cols, event_col, time_col, horizon, cc_only,
             cause_coded, label):
    """PSU-grouped 5-fold CV AUC for one feature set and outcome.

    Imputation is fit on each training fold only, so no test-fold information
    leaks into the model.
    """
    time_v = df[time_col].astype(float).values
    ev_v = df[event_col].astype(int).values
    # Horizon-cap the event: any event past the horizon is treated as censored.
    y = ((ev_v == 1) & (time_v <= horizon)).astype(int)
    mask = ~np.isnan(time_v)
    if cc_only:
        mask = mask & cause_coded
    if mask.sum() < 100 or y[mask].sum() < 5:
        print(f"XGBoost {label}: insufficient events ({y[mask].sum()})")
        return None

    Xm = df.loc[mask, feature_cols].reset_index(drop=True)
    ym = y[mask]
    groups = df.loc[mask, "nhanes_cluster"].values

    aucs = []
    for train_idx, test_idx in GroupKFold(n_splits=5).split(Xm, ym, groups):
        if ym[test_idx].sum() == 0:
            continue
        train_median = Xm.iloc[train_idx].median(numeric_only=True)
        x_train = Xm.iloc[train_idx].fillna(train_median)
        x_test = Xm.iloc[test_idx].fillna(train_median)
        model = xgb.XGBClassifier(
            max_depth=4,
            n_estimators=200,
            learning_rate=0.05,
            subsample=0.8,
            eval_metric="auc",
        )
        model.fit(x_train, ym[train_idx])
        proba = model.predict_proba(x_test)[:, 1]
        aucs.append(roc_auc_score(ym[test_idx], proba))

    if not aucs:
        print(f"XGBoost {label}: no fold had test events")
        return None

    return {
        "mean_auc": float(np.mean(aucs)),
        "std_auc": float(np.std(aucs)),
        "n_events": int(ym.sum()),
        "n": int(mask.sum()),
        "fold_aucs": aucs,
    }


def main():
    ensure_csv()
    df = pd.read_csv(CSV_PATH)
    df = df.dropna(subset=["sdmvpsu", "sdmvstra"])
    df["sex_int"] = (df["sex"].astype(str).str.lower() == "male").astype(int)
    # NHANES has only 3 PSUs nested within ~148 strata. PSU alone is too
    # coarse for 5-fold grouping; use stratum * PSU which gives ~301 groups
    # and respects the design's primary sampling clusters.
    df["nhanes_cluster"] = (
        df["sdmvstra"].astype(int).astype(str)
        + "_"
        + df["sdmvpsu"].astype(int).astype(str)
    )
    df = df.reset_index(drop=True)

    os.makedirs("results", exist_ok=True)
    summary_rows = []

    cause_coded = (
        df["cause_coded"].astype(str).str.upper() == "TRUE"
    ).values if "cause_coded" in df.columns else np.ones(len(df), dtype=bool)

    for fs_name, feature_cols in FEATURE_SETS:
        for event_col, time_col, horizon, label, cc_only in horizon_labels(df):
            if event_col not in df.columns or time_col not in df.columns:
                print(f"Skipping {fs_name}/{label}: column missing")
                continue

            res = evaluate(df, feature_cols, event_col, time_col, horizon,
                           cc_only, cause_coded, f"{fs_name}/{label}")
            if res is None:
                continue

            print(
                f"XGBoost {fs_name}/{label}: n={res['n']} "
                f"events={res['n_events']} "
                f"AUC mean={res['mean_auc']:.3f} std={res['std_auc']:.3f}"
            )
            with open(f"results/xgboost_{fs_name}_{label}.txt", "w") as f:
                f.write(f"feature_set={fs_name}\n")
                f.write(f"features={','.join(feature_cols)}\n")
                f.write(f"n={res['n']}\n")
                f.write(f"events={res['n_events']}\n")
                f.write(f"horizon_y={horizon}\n")
                f.write(f"mean_auc={res['mean_auc']:.4f}\n")
                f.write(f"std_auc={res['std_auc']:.4f}\n")
                f.write(f"fold_aucs={res['fold_aucs']}\n")
                f.write(f"{CAVEAT_UNWEIGHTED}\n")
                f.write(f"{CAVEAT_BINARY}\n")
            summary_rows.append({
                "feature_set": fs_name,
                "label": label,
                "horizon_y": horizon,
                "n": res["n"],
                "events": res["n_events"],
                "mean_auc": res["mean_auc"],
                "std_auc": res["std_auc"],
            })

    if summary_rows:
        pd.DataFrame(summary_rows).to_csv(
            "results/xgboost_summary.csv", index=False
        )


if __name__ == "__main__":
    main()
