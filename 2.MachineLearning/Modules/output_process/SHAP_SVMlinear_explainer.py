import shap
import joblib
import numpy as np
import copy
import pandas as pd
import matplotlib.pyplot as plt
import sys

if __name__ == "__main__":
    # Take ONLY the first argument as base_dir
    if len(sys.argv) < 2:
        print("❌ Error: Please provide base directory path")
        print("Usage: python SHAP_SVMlinear_explainer.py /path/to/directory")
        sys.exit(1)

    base_dir = sys.argv[1]  # This is a string, not a list

    # =========================
    # Load trained model (dictionary)
    # =========================
    saved = joblib.load(
        f"{base_dir}/models/SVM_Linear_final.pkl"
    )

    pipeline = saved["model"]
    scaler = pipeline.named_steps["scaler"]
    svm_model = pipeline.named_steps["model"]

    # =========================
    # Feature definitions
    # =========================
    features = ["has-miR-223-3p", "rsRNA-28s-12", "tsRNA-Asp-GTC-1", "has-let-7f-5p"]

    # =========================
    # Load data
    # =========================
    X_train = pd.read_csv(
        f"{base_dir}/../Data/training_set.csv",
        usecols=features
    )

    X_test_internal = pd.read_csv(
        f"{base_dir}/../Data/internal_testing_set.csv",
        usecols=features
    )

    X_test_external = pd.read_csv(
        f"{base_dir}/../Data/external_validation_cohort.csv",
        usecols=features
    )

    # =========================
    # Scale data — wrap back into DataFrames to preserve column names
    # =========================
    X_train_scaled = pd.DataFrame(
        scaler.transform(X_train),
        columns=features
    )

    X_test_internal_scaled = pd.DataFrame(
        scaler.transform(X_test_internal),
        columns=features
    )

    X_test_external_scaled = pd.DataFrame(
        scaler.transform(X_test_external),
        columns=features
    )

    # =========================
    # SHAP explainer — built on scaled training data (DataFrame)
    # =========================
    explainer = shap.LinearExplainer(svm_model, X_train_scaled)

    # SHAP values
    shap_int = explainer(X_test_internal_scaled)
    shap_ext = explainer(X_test_external_scaled)

    # =========================
    # Rename feature names for plotting
    # =========================
    def rename_shap_features(shap_values, rename_dict):
        new_names = [rename_dict.get(name, name) for name in shap_values.feature_names]
        shap_values.feature_names = new_names
        return shap_values

    # =========================
    # Clip SHAP values (deepcopy to preserve originals)
    # =========================
    def clip_shap_values(shap_values, percentile=99):
        """Clip extreme SHAP values to improve plot readability. All samples retained."""
        vals = shap_values.values
        limit = np.percentile(np.abs(vals), percentile)
        shap_values.values = np.clip(vals, -limit, limit)
        return shap_values

    shap_int_clipped = clip_shap_values(copy.deepcopy(shap_int), percentile=99)
    shap_ext_clipped = clip_shap_values(copy.deepcopy(shap_ext), percentile=99)

    # =========================
    # Investigate outlier in internal test set
    # =========================
    shap_df = pd.DataFrame(shap_int.values, columns=shap_int.feature_names)
    print("Max absolute SHAP value per feature:")
    print(shap_df.abs().max())
    print("\nSample index with max absolute SHAP per feature:")
    print(shap_df.abs().idxmax())
    outlier_idx = shap_df.abs().sum(axis=1).idxmax()
    print(f"\nOverall outlier sample index: {outlier_idx}")
    print(f"SHAP values:\n{shap_df.loc[outlier_idx]}")
    print(f"Raw feature values:\n{X_test_internal.iloc[outlier_idx]}")

    # =========================
    # SHAP Swarm Plot – INTERNAL TEST (unclipped, for supplementary)
    # =========================
    plt.figure(figsize=(6, 4))
    shap.plots.beeswarm(shap_int, max_display=len(features), show=False)
    plt.title("SHAP Swarm Plot (Internal Testing Cohort)", fontsize=9)
    plt.tight_layout()
    plt.savefig(
        f"{base_dir}/vis/SHAP/SHAP_swarm_internal_testing_SVM_linear_unclipped.pdf",
        format="pdf", bbox_inches="tight"
    )
    plt.close()

    # =========================
    # SHAP Swarm Plot – INTERNAL TEST (clipped, for main figure)
    # =========================
    plt.figure(figsize=(6, 4))
    shap.plots.beeswarm(shap_int_clipped, max_display=len(features), show=False)
    plt.title("SHAP Swarm Plot (Internal Testing Cohort)", fontsize=9)
    plt.tight_layout()
    plt.savefig(
        f"{base_dir}/vis/SHAP/SHAP_swarm_internal_testing_SVM_linear_clipped.pdf",
        format="pdf", bbox_inches="tight"
    )
    plt.close()

    # =========================
    # SHAP Swarm Plot – EXTERNAL TEST (unclipped, for supplementary)
    # =========================
    plt.figure(figsize=(6, 4))
    shap.plots.beeswarm(shap_ext, max_display=len(features), show=False)
    plt.title("SHAP Swarm Plot (External Testing Cohort)", fontsize=9)
    plt.tight_layout()
    plt.savefig(
        f"{base_dir}/vis/SHAP/SHAP_swarm_external_testing_SVM_linear_unclipped.pdf",
        format="pdf", bbox_inches="tight"
    )
    plt.close()

    # =========================
    # SHAP Swarm Plot – EXTERNAL TEST (clipped, for main figure)
    # =========================
    plt.figure(figsize=(6, 4))
    shap.plots.beeswarm(shap_ext_clipped, max_display=len(features), show=False)
    plt.title("SHAP Swarm Plot (External Testing Cohort)", fontsize=9)
    plt.tight_layout()
    plt.savefig(
        f"{base_dir}/vis/SHAP/SHAP_swarm_external_testing_SVM_linear_clipped.pdf",
        format="pdf", bbox_inches="tight"
    )
    plt.close()

    print("SHAP analysis completed successfully.")