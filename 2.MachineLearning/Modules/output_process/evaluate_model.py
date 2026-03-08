import argparse
import pandas as pd
import numpy as np
import joblib

from sklearn.metrics import roc_auc_score, confusion_matrix, brier_score_loss
from collections import defaultdict

# --------------------------------------------------
# bootstrap CI
# --------------------------------------------------

from sklearn.utils import resample

def stratified_bootstrap_ci_sklearn(y_true, y_prob, threshold, n=2000, seed=42):
    """
    Stratified bootstrap using sklearn's resample with stratify option
    """
    rng = np.random.default_rng(seed)
    boot_stats = defaultdict(list)
    
    for _ in range(n):
        # Stratified bootstrap sampling
        idx = resample(
            range(len(y_true)), 
            replace=True, 
            n_samples=len(y_true),
            stratify=y_true,  # This ensures class balance is preserved
            random_state=rng.integers(0, 1000000)
        )
        
        yt = y_true[idx]
        yp = y_prob[idx]
        yb = (yp >= threshold).astype(int)
        
        # Calculate metrics
        tn, fp, fn, tp = confusion_matrix(yt, yb).ravel()
        
        boot_stats["AUC"].append(roc_auc_score(yt, yp))
        boot_stats["Sensitivity"].append(tp/(tp+fn) if (tp+fn)>0 else np.nan)
        boot_stats["Specificity"].append(tn/(tn+fp) if (tn+fp)>0 else np.nan)
        boot_stats["PPV"].append(tp/(tp+fp) if (tp+fp)>0 else np.nan)
        boot_stats["NPV"].append(tn/(tn+fn) if (tn+fn)>0 else np.nan)
    
    # Calculate confidence intervals
    ci = {}
    for k, v in boot_stats.items():
        arr = np.array(v)
        arr = arr[~np.isnan(arr)]
        ci[f"{k}_CI_lower"] = np.percentile(arr, 2.5)
        ci[f"{k}_CI_upper"] = np.percentile(arr, 97.5)
    
    return ci

# --------------------------------------------------
# evaluate dataset
# --------------------------------------------------
def evaluate_dataset(csv_file, pipeline, threshold):
    df = pd.read_csv(csv_file)
    
    y_true = df["pred_label"].values
    
    # Get expected feature names from the scaler
    if hasattr(pipeline.named_steps['scaler'], 'feature_names_in_'):
        expected_features = list(pipeline.named_steps['scaler'].feature_names_in_)
    else:
        expected_features = ['miR-223', 'rsRNA-28s', 'tsRNA-Asp-GTC', 'let-7f']
    
    # Select only the features the model expects
    X = df[expected_features]  # Keep as DataFrame, don't use .values
    
    probs = pipeline.predict_proba(X)[:,1]  # This will now use feature names
    preds = (probs >= threshold).astype(int)

    tn, fp, fn, tp = confusion_matrix(y_true, preds).ravel()
    auc = roc_auc_score(y_true, probs)
    sens = tp/(tp+fn)
    spec = tn/(tn+fp)
    ppv  = tp/(tp+fp) if (tp+fp)>0 else np.nan
    npv  = tn/(tn+fn) if (tn+fn)>0 else np.nan

    ci = stratified_bootstrap_ci_sklearn(y_true, probs, threshold)

    result = {
        "Task": csv_file.split("/")[-1],
        "AUC": auc,
        "Threshold": threshold,
        "Sensitivity": sens,
        "Specificity": spec,
        "PPV": ppv,
        "NPV": npv,
        **ci
    }

    return result


# --------------------------------------------------
# main
# --------------------------------------------------

def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--model", required=True)
    parser.add_argument("--csvs", nargs="+", required=True)
    parser.add_argument("--output", default="summary_results.csv")
    args = parser.parse_args()

    print("Loading model...")
    model_artifact = joblib.load(args.model)
    pipeline = model_artifact["model"]
    threshold = float(model_artifact["threshold"])
    
    print(f"Model type: {type(pipeline)}")
    print(f"Threshold: {threshold}")
    
    # Print pipeline steps for debugging
    if hasattr(pipeline, 'named_steps'):
        print("Pipeline steps:")
        for step_name, step_obj in pipeline.named_steps.items():
            print(f"  - {step_name}: {type(step_obj)}")
            
        # Check if scaler has feature names
        if 'scaler' in pipeline.named_steps:
            if hasattr(pipeline.named_steps['scaler'], 'feature_names_in_'):
                print(f"Scaler expected features: {list(pipeline.named_steps['scaler'].feature_names_in_)}")

    results = []
    for csv_file in args.csvs:
        try:
            res = evaluate_dataset(csv_file, pipeline, threshold)
            results.append(res)
        except Exception as e:
            print(f"Error processing {csv_file}: {e}")
            continue

    df = pd.DataFrame(results)
    df.to_excel(args.output, index=False)
    print("\nFinal Results:")
    print(df)

if __name__ == "__main__":
    main()