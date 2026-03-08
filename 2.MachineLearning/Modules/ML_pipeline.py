#!/usr/bin/env python3
import os
import warnings
from collections import defaultdict
from math        import floor, log10
from itertools   import groupby

import numpy  as np
import pandas as pd
import joblib

from sklearn.calibration           import calibration_curve
from sklearn.discriminant_analysis import LinearDiscriminantAnalysis
from sklearn.ensemble              import RandomForestClassifier
from sklearn.linear_model          import LogisticRegression
from sklearn.metrics               import (brier_score_loss,
                                           roc_auc_score, roc_curve)
from sklearn.model_selection       import (RepeatedStratifiedKFold,
                                           StratifiedKFold)
from sklearn.pipeline              import Pipeline
from sklearn.preprocessing         import StandardScaler
from sklearn.svm                   import SVC

from skopt       import BayesSearchCV
from skopt.space import Integer, Real

warnings.filterwarnings("ignore", message=".*Parameters:.*are not used.*")
warnings.filterwarnings("ignore", message="The objective has been evaluated.*")
warnings.filterwarnings("ignore", message=".*did not converge.*",
                        category=UserWarning, module="skopt")


# =============================================================================
# 0.  GLOBAL CONFIGURATION
# =============================================================================
RANDOM_STATE      = 42
OUTER_FOLDS       = 5
N_REPEATS         = 10
INNER_FOLDS       = 3
N_BAYES_ITER      = 25
MIN_PARAM_SUPPORT = 1
N_DECIMALS        = 4

ROOT_DIR  = "./Outputs"
PRED_DIR  = os.path.join(ROOT_DIR, "predictions")
STAT_DIR  = os.path.join(ROOT_DIR, "statistics")
CALIB_DIR = os.path.join(ROOT_DIR, "calibration")
MODEL_DIR = os.path.join(ROOT_DIR, "models")

for _d in [PRED_DIR, STAT_DIR, CALIB_DIR, MODEL_DIR]:
    os.makedirs(_d, exist_ok=True)


# =============================================================================
# 1.  MODEL REGISTRY
# =============================================================================
MODEL_CONFIGS: dict = {

    "LDA_lsqr": {
        "estimator": LinearDiscriminantAnalysis(solver="lsqr"),
        "scaling"  : True,
        "space"    : {"model__shrinkage": [None, "auto"]},
    },

    "LogReg_L2": {
        "estimator": LogisticRegression(
            solver="lbfgs", penalty="l2",
            max_iter=5000, random_state=RANDOM_STATE
        ),
        "scaling"  : True,
        "space"    : {
            "model__C"           : Real(1e-3, 1e2, prior="log-uniform"),
            "model__class_weight": [None, "balanced"],
        },
    },

    "SVM_Linear": {
        "estimator": SVC(kernel="linear", probability=True,
                         random_state=RANDOM_STATE),
        "scaling"  : True,
        "space"    : {
            "model__C"           : Real(1e-2, 1e2, prior="log-uniform"),
            "model__class_weight": [None, "balanced"],
        },
    },

    "RandomForest": {
        "estimator": RandomForestClassifier(n_jobs=-1,
                                            random_state=RANDOM_STATE),
        "scaling"  : False,
        "space"    : {
            "model__n_estimators"    : Integer(50, 150),
            "model__max_depth"       : Integer(2, 4),
            "model__min_samples_leaf": Integer(5, 15),
            "model__max_features"    : ["sqrt", 0.5, 0.8],
            "model__class_weight"    : [None, "balanced"],
        },
    },
}


# =============================================================================
# 2.  UTILITIES
# =============================================================================

def _fmt(x: float) -> str:
    return f"{x:.{N_DECIMALS}f}"


def _build_pipeline(cfg: dict) -> Pipeline:
    steps = [("scaler", StandardScaler())] if cfg["scaling"] else []
    steps.append(("model", cfg["estimator"]))
    return Pipeline(steps)


def _youden_threshold(y_true: np.ndarray, y_prob: np.ndarray) -> tuple:
    """Return (threshold, sensitivity, specificity, Youden-J)."""
    fpr, tpr, thr = roc_curve(y_true, y_prob)
    idx = int(np.argmax(tpr - fpr))
    return thr[idx], tpr[idx], 1.0 - fpr[idx], tpr[idx] - fpr[idx]


def _round_sig(x: float, sig: int = 2) -> float:
    if x == 0:
        return 0.0
    magnitude = floor(log10(abs(x)))
    factor    = 10 ** (sig - 1 - magnitude)
    return round(x * factor) / factor


def _bin_params(params: dict) -> tuple:
    binned = {}
    for k, v in params.items():
        if isinstance(v, bool):
            binned[k] = str(v)
        elif isinstance(v, float):
            binned[k] = _round_sig(v, sig=2)
        elif isinstance(v, int):
            binned[k] = v
        else:
            binned[k] = str(v)
    return tuple(sorted(binned.items()))


def _select_optimal_params(records: list) -> dict:
    if not records or all(r["_params"] == {} for r in records):
        return {}

    grouped_aucs: dict = defaultdict(list)
    grouped_reps: dict = {}

    for r in records:
        key = _bin_params(r["_params"])
        grouped_aucs[key].append(r["_auc"])
        if key not in grouped_reps:
            grouped_reps[key] = r["_params"]

    eligible = {k: v for k, v in grouped_aucs.items()
                if len(v) >= MIN_PARAM_SUPPORT}

    if eligible:
        best_key = max(eligible, key=lambda k: float(np.mean(eligible[k])))
        strategy = f"highest-mean-AUC (support >= {MIN_PARAM_SUPPORT})"
    else:
        best_key = max(grouped_aucs,
                       key=lambda k: (len(grouped_aucs[k]),
                                      float(np.mean(grouped_aucs[k]))))
        strategy = (f"FALLBACK: no config met support >= {MIN_PARAM_SUPPORT}; "
                    f"selected most frequent")

    best_count    = len(grouped_aucs[best_key])
    best_mean_auc = float(np.mean(grouped_aucs[best_key]))
    n_groups      = len(grouped_aucs)

    print(f"      [param selection]  {n_groups} distinct configs / {len(records)} folds"
          f"  |  {strategy}"
          f"  |  n={best_count} folds, mean AUC={_fmt(best_mean_auc)}")

    return grouped_reps[best_key]


# =============================================================================
# 3.  PHASE 1 -- REPEATED NESTED CROSS-VALIDATION
# =============================================================================

def nested_cv(X: pd.DataFrame, y: pd.Series) -> dict:
    rskf = RepeatedStratifiedKFold(
        n_splits=OUTER_FOLDS, n_repeats=N_REPEATS, random_state=RANDOM_STATE
    )

    oof_records: list = []
    oof_preds  : dict = defaultdict(list)

    for model_name, cfg in MODEL_CONFIGS.items():
        print(f"\n{'='*65}\n  Model : {model_name}\n{'='*65}")

        for fold_idx, (tr_idx, va_idx) in enumerate(rskf.split(X, y)):
            repeat = fold_idx // OUTER_FOLDS + 1
            fold   = fold_idx  % OUTER_FOLDS + 1
            print(f"  Repeat {repeat:2d}/{N_REPEATS}  Fold {fold}/{OUTER_FOLDS}",
                  end="  ...  ", flush=True)

            X_tr, X_va = X.iloc[tr_idx], X.iloc[va_idx]
            y_tr, y_va = y.iloc[tr_idx], y.iloc[va_idx]

            inner_cv  = StratifiedKFold(INNER_FOLDS, shuffle=True,
                                        random_state=RANDOM_STATE + fold_idx)
            base_pipe = _build_pipeline(cfg)

            if cfg["space"]:
                searcher = BayesSearchCV(
                    base_pipe, cfg["space"],
                    n_iter=N_BAYES_ITER, scoring="roc_auc",
                    cv=inner_cv, n_jobs=-1,
                    random_state=RANDOM_STATE + fold_idx,
                    refit=True,
                )
                searcher.fit(X_tr, y_tr)
                best_params = dict(searcher.best_params_)
            else:
                base_pipe.fit(X_tr, y_tr)
                best_params = {}

            uncalib_pipe = _build_pipeline(cfg)
            if best_params:
                uncalib_pipe.set_params(**best_params)
            uncalib_pipe.fit(X_tr, y_tr)

            probs    = uncalib_pipe.predict_proba(X_va)[:, 1]
            fold_auc = round(roc_auc_score(y_va, probs), N_DECIMALS)
            print(f"AUC = {_fmt(fold_auc)}  |  "
                  f"params: {best_params if best_params else 'fixed'}")

            oof_records.append({
                "Model": model_name, "Repeat": repeat, "Fold": fold,
                "_params": best_params, "_auc": fold_auc,
            })
            oof_preds[model_name].append(pd.DataFrame({
                "SampleID" : X_va.index,
                "TrueLabel": y_va.values,
                "PredProb" : probs.round(N_DECIMALS),
                "Repeat"   : repeat,
                "Fold"     : fold,
            }))

    oof_pred_dfs = {m: pd.concat(dfs, ignore_index=True)
                    for m, dfs in oof_preds.items()}

    oof_path = os.path.join(PRED_DIR, "nested_cv_oof_predictions.xlsx")
    with pd.ExcelWriter(oof_path) as w:
        for m, df in oof_pred_dfs.items():
            df.to_excel(w, sheet_name=m[:31], index=False)

    pd.DataFrame([
        {"Model": r["Model"], "Repeat": r["Repeat"], "Fold": r["Fold"],
         "AUC": _fmt(r["_auc"]), "Params": str(r["_params"])}
        for r in oof_records
    ]).to_csv(os.path.join(STAT_DIR, "nested_cv_fold_performance.csv"), index=False)

    print(f"\n[v] OOF predictions  ->  {oof_path}")
    print(f"[v] Fold performance ->  {STAT_DIR}/nested_cv_fold_performance.csv")
    return {"oof_records": oof_records, "oof_pred_dfs": oof_pred_dfs}


# =============================================================================
# 4.  PHASE 2 -- HYPERPARAMETER SELECTION
# =============================================================================

def select_stable_hyperparameters(oof_records: list) -> dict:
    print(f"\n{'='*65}")
    print(f"  Phase 2 -- Hyperparameter Selection")
    print(f"  (binned highest-mean AUC, MIN_PARAM_SUPPORT = {MIN_PARAM_SUPPORT})")
    print(f"{'='*65}")

    stable_params: dict = {}
    summary_rows : list = []

    for model_name, group in groupby(
            sorted(oof_records, key=lambda r: r["Model"]),
            key=lambda r: r["Model"]):

        records = list(group)
        aucs    = [r["_auc"] for r in records]
        params  = _select_optimal_params(records)
        stable_params[model_name] = params

        summary_rows.append({
            "Model"         : model_name,
            "Mean_AUC"      : _fmt(np.mean(aucs)),
            "SD_AUC"        : _fmt(np.std(aucs)),
            "Median_AUC"    : _fmt(np.median(aucs)),
            "Min_AUC"       : _fmt(np.min(aucs)),
            "Max_AUC"       : _fmt(np.max(aucs)),
            "SelectedParams": str(params),
        })
        print(f"  {model_name:25s}  "
              f"AUC = {_fmt(np.mean(aucs))} +/- {_fmt(np.std(aucs))}"
              f"   params: {params if params else 'fixed'}")

    (pd.DataFrame(summary_rows)
       .sort_values("Mean_AUC", ascending=False)
       .to_csv(os.path.join(STAT_DIR, "model_performance_summary.csv"), index=False))
    print(f"\n[v] Performance summary ->  {STAT_DIR}/model_performance_summary.csv")
    return stable_params


# =============================================================================
# 5.  PHASE 3 -- FINAL MODEL TRAINING ON FULL TRAINING SET
# =============================================================================

def train_final_models(X            : pd.DataFrame,
                       y            : pd.Series,
                       stable_params: dict,
                       oof_pred_dfs : dict) -> dict:
    print(f"\n{'='*65}")
    print(f"  Phase 3 -- Final Model Training  (FULL cohort)")
    print(f"{'='*65}")

    final_artifacts: dict = {}
    calib_summary  : list = []

    for model_name, cfg in MODEL_CONFIGS.items():
        print(f"\n  {model_name}", flush=True)
        params    = stable_params.get(model_name, {})
        base_pipe = _build_pipeline(cfg)
        if params:
            base_pipe.set_params(**params)
        base_pipe.fit(X, y)

        oof_df   = oof_pred_dfs[model_name]
        oof_mean = (oof_df
                    .groupby("SampleID")
                    .agg(TrueLabel=("TrueLabel", "first"),
                         PredProb =("PredProb",  "mean"))
                    .reset_index())

        thr, sen, spe, yj = _youden_threshold(
            oof_mean["TrueLabel"].values, oof_mean["PredProb"].values
        )

        pt, pp = calibration_curve(oof_mean["TrueLabel"].values,
                                   oof_mean["PredProb"].values,
                                   n_bins=10, strategy="quantile")
        b_unc = round(brier_score_loss(oof_mean["TrueLabel"].values,
                                       oof_mean["PredProb"].values), N_DECIMALS)
        e_unc = round(float(np.mean(np.abs(pt - pp))), N_DECIMALS)

        pd.DataFrame([{"Condition": "Uncalibrated", "Brier": b_unc,
                        "ECE": e_unc, "Model": model_name}]).to_csv(
            os.path.join(CALIB_DIR, f"{model_name}_calibration.csv"), index=False
        )

        artifact = {
            "model"       : base_pipe,
            "threshold"   : round(float(thr), N_DECIMALS),
            "youden_stats": {"Threshold": _fmt(thr), "Sensitivity": _fmt(sen),
                             "Specificity": _fmt(spe), "Youden_J": _fmt(yj)},
            "stable_params": params,
        }
        final_artifacts[model_name] = artifact
        joblib.dump(artifact, os.path.join(MODEL_DIR, f"{model_name}_final.pkl"))

        calib_summary.append({
            "Model": model_name, "Threshold": _fmt(thr),
            "Sensitivity": _fmt(sen), "Specificity": _fmt(spe),
            "Youden_J": _fmt(yj), "Brier_Uncal": _fmt(b_unc),
            "ECE_Uncal": _fmt(e_unc), "SelectedParams": str(params),
        })
        print(f"    Threshold={_fmt(thr)}  Sen={_fmt(sen)}  "
              f"Spe={_fmt(spe)}  Youden={_fmt(yj)}")
        print(f"    Brier(uncal)={_fmt(b_unc)}  ECE(uncal)={_fmt(e_unc)}")

    pd.DataFrame(calib_summary).to_csv(
        os.path.join(CALIB_DIR, "calibration_summary.csv"), index=False
    )
    print(f"\n[v] Final models  ->  {MODEL_DIR}")
    print(f"[v] Uncal metrics ->  {CALIB_DIR}")
    return final_artifacts


# =============================================================================
# 6.  PHASE 4 -- TESTING / VALIDATION  (probabilities only)
# =============================================================================

def _predict_and_save(X: pd.DataFrame, y: pd.Series,
                      final_artifacts: dict,
                      label: str, out_csv: str, out_xlsx: str) -> None:
    """Shared logic: predict probabilities, save predictions, no metrics."""
    print(f"\n{'='*65}\n  Phase 4 -- {label}  (frozen models)\n{'='*65}")

    preds: dict = {}
    for model_name, artifact in final_artifacts.items():
        probs = artifact["model"].predict_proba(X)[:, 1]
        preds[model_name] = pd.DataFrame({
            "SampleID" : X.index,
            "TrueLabel": y.values,
            "PredProb" : probs.round(N_DECIMALS),
        })
        print(f"  {model_name:25s}  done")

    # Wide format: one row per sample, one column per model
    wide = pd.DataFrame({"SampleID": X.index, "TrueLabel": y.values})
    for model_name, df in preds.items():
        wide[model_name] = df["PredProb"].values
    wide.to_csv(out_csv, index=False)

    with pd.ExcelWriter(out_xlsx) as w:
        for model_name, df in preds.items():
            df.to_excel(w, sheet_name=model_name[:31], index=False)

    print(f"[v] Predictions (wide) ->  {out_csv}")
    print(f"[v] Predictions (long) ->  {out_xlsx}")


def internal_testing(X_int: pd.DataFrame, y_int: pd.Series,
                     final_artifacts: dict) -> None:
    _predict_and_save(
        X_int, y_int, final_artifacts,
        label    = "Internal Testing",
        out_csv  = os.path.join(STAT_DIR, "internal_testing_predictions.csv"),
        out_xlsx = os.path.join(PRED_DIR, "internal_testing_predictions.xlsx"),
    )


def external_validation(X_ext: pd.DataFrame, y_ext: pd.Series,
                        final_artifacts: dict) -> None:
    _predict_and_save(
        X_ext, y_ext, final_artifacts,
        label    = "External Validation",
        out_csv  = os.path.join(STAT_DIR, "external_validation_predictions.csv"),
        out_xlsx = os.path.join(PRED_DIR, "external_validation_predictions.xlsx"),
    )


# =============================================================================
# 7.  ENTRY POINT
# =============================================================================

if __name__ == "__main__":

    derivation = pd.read_csv(
        "./Data/training_set.csv",
        index_col=0)
    internal   = pd.read_csv(
        "./Data/internal_testing_set.csv",
        index_col=0)
    external   = pd.read_csv(
        "./Data/external_validation_cohort.csv",
        index_col=0)

    FEATURES = ["has-miR-223-3p", "rsRNA-28s-12", "tsRNA-Asp-GTC-1", "has-let-7f-5p"]
    LABEL    = "label"

    X_dev, y_dev = derivation[FEATURES], derivation[LABEL]
    X_int, y_int = internal[FEATURES],   internal[LABEL]
    X_ext, y_ext = external[FEATURES],   external[LABEL]

    print("\n" + "#"*65)
    print(f"  PHASE 1  --  Repeated Nested CV  "
          f"({N_REPEATS} x {OUTER_FOLDS}-fold = {N_REPEATS*OUTER_FOLDS} outer folds)")
    print("#"*65)
    cv_results = nested_cv(X_dev, y_dev)

    print("\n" + "#"*65 + "\n  PHASE 2  --  Hyperparameter Selection\n" + "#"*65)
    stable_params = select_stable_hyperparameters(cv_results["oof_records"])

    print("\n" + "#"*65 + "\n  PHASE 3  --  Final Model Training\n" + "#"*65)
    final_artifacts = train_final_models(
        X_dev, y_dev, stable_params, cv_results["oof_pred_dfs"]
    )

    print("\n" + "#"*65 + "\n  PHASE 4  --  Internal Testing\n" + "#"*65)
    internal_testing(X_int, y_int, final_artifacts)

    print("\n" + "#"*65 + "\n  PHASE 4  --  External Validation\n" + "#"*65)
    external_validation(X_ext, y_ext, final_artifacts)

    print("\n" + "="*65 + "\n  PIPELINE COMPLETE\n" + "="*65)