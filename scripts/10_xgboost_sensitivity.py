#!/usr/bin/env python3
"""scripts/10_xgboost_sensitivity.py

XGBoost ML upper-bound baseline as a sensitivity analysis. Bounds the
achievable AUC on the 5 MetS factors + age + sex, providing a ceiling
against which the clinical risk scores can be interpreted.

Inputs:  data/processed/cohort_with_scores.rds (canonical)
         data/processed/cohort_with_scores.csv (generated from RDS via Rscript
         when missing, since the arrow R package is intentionally not installed)
Outputs: results/xgboost_{allcause,cv,dm}.txt

NHANES PSU-grouped CV respects the sampling design and avoids optimism from
within-PSU resampling. Horizons match the survival scripts: censor at 9.5y
for allcause / cv (followup_years), 14.5y for dm (followup_years_dm), with
events counted only inside the horizon.
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
    return [
        ("event_allcause", "followup_years", 9.5, "allcause_9_5y"),
        ("event_allcause", "followup_years", 5.0, "allcause_5y"),
        ("event_cv",       "followup_years", 9.5, "cv_9_5y"),
        ("event_cv",       "followup_years", 5.0, "cv_5y"),
        ("event_dm",       "followup_years_dm", 14.5, "dm_14_5y"),
        ("event_dm",       "followup_years_dm", 10.0, "dm_10y"),
        ("event_dm",       "followup_years_dm", 5.0,  "dm_5y"),
    ]


def main():
    ensure_csv()
    df = pd.read_csv(CSV_PATH)

    features = [
        "age", "waist_cm", "sbp", "dbp", "hdl",
        "fasting_glucose", "triglycerides", "sex"
    ]
    df = df.dropna(subset=features + ["sdmvpsu", "sdmvstra"])
    df["sex_int"] = (df["sex"].astype(str).str.lower() == "male").astype(int)
    # NHANES has only 3 PSUs nested within ~148 strata. PSU alone is too
    # coarse for 5-fold grouping; use stratum * PSU which gives ~301 groups
    # and respects the design's primary sampling clusters.
    df["nhanes_cluster"] = (
        df["sdmvstra"].astype(int).astype(str)
        + "_"
        + df["sdmvpsu"].astype(int).astype(str)
    )

    X = df[
        ["age", "waist_cm", "sbp", "dbp", "hdl",
         "fasting_glucose", "triglycerides", "sex_int"]
    ].copy()
    X = X.fillna(X.median(numeric_only=True))
    X.reset_index(drop=True, inplace=True)
    df = df.reset_index(drop=True)

    os.makedirs("results", exist_ok=True)
    summary_rows = []

    for event_col, time_col, horizon, label in horizon_labels(df):
        if event_col not in df.columns or time_col not in df.columns:
            print(f"Skipping {label}: column missing")
            continue

        # Horizon-cap the event: any event with follow-up > horizon is censored
        time_v = df[time_col].astype(float).values
        ev_v = df[event_col].astype(int).values
        in_window = time_v <= horizon
        y = (ev_v == 1) & in_window
        y = y.astype(int)

        # Restrict CV to subjects with at least some follow-up (drop NA times)
        mask = ~np.isnan(time_v)
        if mask.sum() < 100 or y[mask].sum() < 5:
            print(f"XGBoost {label}: insufficient events ({y[mask].sum()})")
            continue

        Xm = X.loc[mask].reset_index(drop=True)
        ym = y[mask]
        groups = df.loc[mask, "nhanes_cluster"].values

        # PSU-grouped 5-fold CV
        gkf = GroupKFold(n_splits=5)
        aucs = []
        for train_idx, test_idx in gkf.split(Xm, ym, groups):
            if ym[test_idx].sum() == 0:
                continue
            model = xgb.XGBClassifier(
                max_depth=4,
                n_estimators=200,
                learning_rate=0.05,
                subsample=0.8,
                eval_metric="auc",
            )
            model.fit(Xm.iloc[train_idx], ym[train_idx])
            proba = model.predict_proba(Xm.iloc[test_idx])[:, 1]
            aucs.append(roc_auc_score(ym[test_idx], proba))

        if not aucs:
            print(f"XGBoost {label}: no fold had test events")
            continue

        mean_auc = float(np.mean(aucs))
        std_auc = float(np.std(aucs))
        n_events = int(ym.sum())
        print(
            f"XGBoost {label}: n={mask.sum()} events={n_events} "
            f"AUC mean={mean_auc:.3f} std={std_auc:.3f}"
        )
        with open(f"results/xgboost_{label}.txt", "w") as f:
            f.write(f"n={int(mask.sum())}\n")
            f.write(f"events={n_events}\n")
            f.write(f"horizon_y={horizon}\n")
            f.write(f"mean_auc={mean_auc:.4f}\n")
            f.write(f"std_auc={std_auc:.4f}\n")
            f.write(f"fold_aucs={aucs}\n")
        summary_rows.append({
            "label": label,
            "horizon_y": horizon,
            "n": int(mask.sum()),
            "events": n_events,
            "mean_auc": mean_auc,
            "std_auc": std_auc,
        })

    if summary_rows:
        pd.DataFrame(summary_rows).to_csv(
            "results/xgboost_summary.csv", index=False
        )


if __name__ == "__main__":
    main()
