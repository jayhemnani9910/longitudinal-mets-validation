#!/usr/bin/env python3
"""scripts/10_xgboost_sensitivity.py

XGBoost ML upper-bound baseline as a sensitivity analysis. Bounds the
achievable AUC on the 5 MetS factors + age + sex, providing a ceiling
against which the clinical risk scores can be interpreted.

Inputs:  data/processed/cohort_with_scores.feather (exported from R)
Outputs: results/xgboost_{allcause,cv,dm}.txt

NHANES PSU-grouped CV to respect sampling design (avoid optimism from
within-PSU resampling).
"""

import os
import sys

import numpy as np
import pandas as pd
import xgboost as xgb
from sklearn.metrics import roc_auc_score
from sklearn.model_selection import GroupKFold


def main():
    feather_path = "data/processed/cohort_with_scores.feather"
    if not os.path.exists(feather_path):
        print(
            f"Missing {feather_path}. "
            "Export it from R first:\n"
            "  R -e 'arrow::write_feather("
            "readRDS(\"data/processed/cohort_with_scores.rds\"), "
            "\"data/processed/cohort_with_scores.feather\")'",
            file=sys.stderr,
        )
        # Alternative without arrow: use pyreadstat to convert
        # The R-side scripts also write an rds; if arrow isn't available we can
        # save as a CSV intermediate (slower but works).
        sys.exit(1)

    df = pd.read_feather(feather_path)

    features = [
        "age", "waist_cm", "sbp", "dbp", "hdl",
        "fasting_glucose", "triglycerides", "sex"
    ]
    df = df.dropna(subset=features + ["sdmvpsu", "event_allcause"])
    df["sex_int"] = (df["sex"] == "male").astype(int)

    X = df[
        ["age", "waist_cm", "sbp", "dbp", "hdl",
         "fasting_glucose", "triglycerides", "sex_int"]
    ].copy()
    X = X.fillna(X.median(numeric_only=True))

    os.makedirs("results", exist_ok=True)

    for outcome_col, outcome_name in [
        ("event_allcause", "allcause"),
        ("event_cv", "cv"),
        ("event_dm", "dm"),
    ]:
        y = df[outcome_col].astype(int).values
        groups = df["sdmvpsu"].astype(str).values

        # PSU-grouped 5-fold CV
        gkf = GroupKFold(n_splits=5)
        aucs = []
        for train_idx, test_idx in gkf.split(X, y, groups):
            model = xgb.XGBClassifier(
                max_depth=4,
                n_estimators=200,
                learning_rate=0.05,
                subsample=0.8,
                eval_metric="auc",
            )
            model.fit(X.iloc[train_idx], y[train_idx])
            proba = model.predict_proba(X.iloc[test_idx])[:, 1]
            # If a fold has no positive class, skip
            if y[test_idx].sum() == 0:
                continue
            aucs.append(roc_auc_score(y[test_idx], proba))

        if not aucs:
            print(f"XGBoost {outcome_name}: insufficient events")
            continue

        mean_auc = float(np.mean(aucs))
        std_auc = float(np.std(aucs))
        print(
            f"XGBoost AUC for {outcome_name}: "
            f"mean={mean_auc:.3f} std={std_auc:.3f}"
        )

        with open(f"results/xgboost_{outcome_name}.txt", "w") as f:
            f.write(f"mean_auc={mean_auc:.4f}\n")
            f.write(f"std_auc={std_auc:.4f}\n")
            f.write(f"fold_aucs={aucs}\n")


if __name__ == "__main__":
    main()
