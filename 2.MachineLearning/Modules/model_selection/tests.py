"""
Statistical Model Comparison Tests
====================================
Performs DeLong's test,
and McNemar's test for pairwise model comparison.

Input:
  - predictions.xlsx : multiple sheets, one per model
      Columns: SampleID, TrueLabel, PredProb
  - thresholds.csv   : columns: Model, Threshold  (used for McNemar's test)

Output:
  - delong_results.csv
  - wilcoxon_results.csv
  - mcnemar_results.csv
"""

import numpy as np
import pandas as pd
from scipy import stats
from scipy.stats import wilcoxon
from itertools import combinations
import warnings
warnings.filterwarnings("ignore")


# ─────────────────────────────────────────────
# 1.  DATA LOADING
# ─────────────────────────────────────────────

def load_data(excel_path: str, threshold_csv_path: str):
    """
    Returns
    -------
    models_data : dict  {model_name: {'y_true': arr, 'y_prob': arr}}
    thresholds  : dict  {model_name: float}
    """
    xls = pd.ExcelFile(excel_path)
    models_data = {}
    for sheet in xls.sheet_names:
        df = xls.parse(sheet)
        # Normalise column names (strip whitespace)
        df.columns = df.columns.str.strip()
        models_data[sheet] = {
            "y_true": df["TrueLabel"].to_numpy(dtype=int),
            "y_prob": df["PredProb"].to_numpy(dtype=float),
        }

    thresh_df = pd.read_csv(threshold_csv_path)
    thresh_df.columns = thresh_df.columns.str.strip()
    thresholds = dict(zip(thresh_df["Model"].str.strip(),
                          thresh_df["Threshold"].astype(float)))

    return models_data, thresholds


# ─────────────────────────────────────────────
# 2.  DELONG TEST
# ─────────────────────────────────────────────

def _compute_auc(y_true, y_pred):
    """Compute AUC via the Mann-Whitney U statistic (exact, no sklearn needed)."""
    pos = y_pred[y_true == 1]
    neg = y_pred[y_true == 0]
    n1, n0 = len(pos), len(neg)
    u = sum(
        (p > n) + 0.5 * (p == n)
        for p in pos
        for n in neg
    )
    return u / (n1 * n0)


def _delong_test(y_true, y_pred1, y_pred2):
    """
    DeLong et al. (1988) test for comparing two correlated AUCs.
    Returns a result dict.
    """
    pos_idx = np.where(y_true == 1)[0]
    neg_idx = np.where(y_true == 0)[0]
    n1, n0 = len(pos_idx), len(neg_idx)

    def structural_components(y_pred):
        pos_preds = y_pred[pos_idx]
        neg_preds = y_pred[neg_idx]
        V10 = np.array([
            np.mean(pos_preds[i] > neg_preds) + 0.5 * np.mean(pos_preds[i] == neg_preds)
            for i in range(n1)
        ])
        V01 = np.array([
            np.mean(pos_preds > neg_preds[j]) + 0.5 * np.mean(pos_preds == neg_preds[j])
            for j in range(n0)
        ])
        return V10, V01

    V10_1, V01_1 = structural_components(y_pred1)
    V10_2, V01_2 = structural_components(y_pred2)

    auc1 = np.mean(V10_1)   # AUC = mean of V10 (equivalent to Mann-Whitney)
    auc2 = np.mean(V10_2)
    auc_diff = auc1 - auc2

    # 2×2 covariance matrices
    S10 = np.cov(np.column_stack([V10_1, V10_2]), rowvar=False, ddof=1)
    S01 = np.cov(np.column_stack([V01_1, V01_2]), rowvar=False, ddof=1)

    var_diff = (
        (S10[0, 0] + S10[1, 1] - 2 * S10[0, 1]) / n1 +
        (S01[0, 0] + S01[1, 1] - 2 * S01[0, 1]) / n0
    )

    if var_diff <= 0:
        se, z_stat, p_value = 0.0, 0.0, 1.0
    else:
        se = np.sqrt(var_diff)
        z_stat = auc_diff / se
        p_value = 2 * (1 - stats.norm.cdf(abs(z_stat)))

    ci_lower = auc_diff - 1.96 * se
    ci_upper = auc_diff + 1.96 * se

    return {
        "AUC_Model1": round(auc1, 4),
        "AUC_Model2": round(auc2, 4),
        "AUC_Diff":   round(auc_diff, 4),
        "SE":         round(se, 4),
        "Z_Statistic": round(z_stat, 4),
        "P_Value":    round(p_value, 4),
        "CI_95_Lower": round(ci_lower, 4),
        "CI_95_Upper": round(ci_upper, 4),
        "Significant": p_value < 0.05,
    }


def run_delong_all_pairs(models_data):
    """Run DeLong test for every pair of models."""
    rows = []
    model_names = list(models_data.keys())

    # All models must share the same test set (same y_true)
    ref_true = models_data[model_names[0]]["y_true"]
    for m in model_names[1:]:
        if not np.array_equal(models_data[m]["y_true"], ref_true):
            raise ValueError(
                f"y_true mismatch between {model_names[0]} and {m}. "
                "All models must be evaluated on the same test set."
            )

    for m1, m2 in combinations(model_names, 2):
        y_true = models_data[m1]["y_true"]
        result = _delong_test(y_true,
                              models_data[m1]["y_prob"],
                              models_data[m2]["y_prob"])
        rows.append({"Model1": m1, "Model2": m2, **result})

    return pd.DataFrame(rows)


# ─────────────────────────────────────────────
# 3.  McNEMAR'S TEST  (unchanged – verified correct)
# ─────────────────────────────────────────────

def _mcnemar_test(y_true, y_prob1, y_prob2, threshold1, threshold2):
    """
    McNemar's test with continuity correction (Edwards, 1948).
    Compares correctness of predictions, not raw labels.
    """
    y_pred1 = (y_prob1 >= threshold1).astype(int)
    y_pred2 = (y_prob2 >= threshold2).astype(int)

    correct1 = (y_pred1 == y_true)
    correct2 = (y_pred2 == y_true)

    a = int(np.sum( correct1 &  correct2))   # both correct
    b = int(np.sum( correct1 & ~correct2))   # only model1 correct
    c = int(np.sum(~correct1 &  correct2))   # only model2 correct
    d = int(np.sum(~correct1 & ~correct2))   # both wrong

    if b + c > 0:
        chi2 = (abs(b - c) - 1) ** 2 / (b + c)
        p_value = 1 - stats.chi2.cdf(chi2, df=1)
    else:
        chi2, p_value = 0.0, 1.0

    odds_ratio = (b / c) if c > 0 else np.inf

    return {
        "a_both_correct":    a,
        "b_only_m1_correct": b,
        "c_only_m2_correct": c,
        "d_both_wrong":      d,
        "Chi2":       round(chi2, 4),
        "P_Value":    round(p_value, 4),
        "Odds_Ratio": round(odds_ratio, 4) if np.isfinite(odds_ratio) else "Inf",
        "Significant": p_value < 0.05,
    }


def run_mcnemar_all_pairs(models_data, thresholds):
    """Run McNemar's test for every pair of models."""
    rows = []
    model_names = list(models_data.keys())

    for m in model_names:
        if m not in thresholds:
            raise KeyError(
                f"No threshold found for model '{m}' in the threshold CSV. "
                f"Available: {list(thresholds.keys())}"
            )

    for m1, m2 in combinations(model_names, 2):
        y_true = models_data[m1]["y_true"]
        result = _mcnemar_test(
            y_true,
            models_data[m1]["y_prob"],
            models_data[m2]["y_prob"],
            thresholds[m1],
            thresholds[m2],
        )
        rows.append({"Model1": m1, "Model2": m2, **result})

    return pd.DataFrame(rows)


# ─────────────────────────────────────────────
# 4.  MAIN
# ─────────────────────────────────────────────

def main(
    excel_path: str = "predictions.xlsx",
    threshold_csv_path: str = "thresholds.csv",
    n_bootstrap: int = 2000,
    output_dir: str = ".",
):
    print("Loading data …")
    models_data, thresholds = load_data(excel_path, threshold_csv_path)
    print(f"  Models loaded: {list(models_data.keys())}")
    print(f"  Thresholds:    {thresholds}\n")

    # --- DeLong ---
    print("Running DeLong test …")
    delong_df = run_delong_all_pairs(models_data)
    delong_path = f"{output_dir}/delong_results.csv"
    delong_df.to_csv(delong_path, index=False)
    print(delong_df.to_string(index=False))
    print(f"  → Saved to {delong_path}\n")

    # --- McNemar ---
    print("Running McNemar's test …")
    mcnemar_df = run_mcnemar_all_pairs(models_data, thresholds)
    mcnemar_path = f"{output_dir}/mcnemar_results.csv"
    mcnemar_df.to_csv(mcnemar_path, index=False)
    print(mcnemar_df.to_string(index=False))
    print(f"  → Saved to {mcnemar_path}\n")

    print("All tests complete.")
    return delong_df, mcnemar_df


if __name__ == "__main__":
    import sys
    
    excel_path        = sys.argv[1] if len(sys.argv) > 1 else "predictions.xlsx"
    threshold_csv     = sys.argv[2] if len(sys.argv) > 2 else "thresholds.csv"
    target_dir        = sys.argv[3] if len(sys.argv) > 3 else "."

    delong_df, mcnemar_df = main(excel_path, threshold_csv, target_dir)
    delong_df.to_csv(f"{target_dir}/delong_results.csv", index=False)
    mcnemar_df.to_csv(f"{target_dir}/mcnemar_results.csv", index=False)
