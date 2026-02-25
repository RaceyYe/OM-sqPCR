# OM-sqPCR Designer

**Primer Design and Machine Learning Framework for OM-sqPCR Applications**

🔗 Repository: [https://github.com/RaceyYe/OM-sqPCR](https://github.com/RaceyYe/OM-sqPCR)

---

## Overview

**OM-sqPCR Designer** is an integrated computational framework that combines:

1. A primer design software package for OM-sqPCR applications
2. Machine learning pipelines for predictive modeling and feature-based classification

The platform enables systematic primer design, optimization, and evaluation, alongside statistically rigorous machine learning workflows for predictive modeling based on OM-sqPCR-derived features.

The repository is structured to support reproducible bioinformatics workflows and modular expansion for future algorithmic development.

---

## Repository Structure

```
OM-sqPCR/
│
├── Primer_Design/
│   ├── App.py
│   ├── core_algorithms/
│   ├── GUI/
│   └── utilities/
│
├── Machine_Learning/
│   ├── preprocessing/
│   ├── feature_engineering/
│   ├── model_training/
│   └── evaluation/
│
└── README.md
```

---

# 1. Primer Design Software

The primer design module implements algorithms tailored for OM-sqPCR primer construction and evaluation.

### Functional Capabilities

* Design OM-sqPCR primers for target sequences
* Analyze thermodynamic properties:

  * Melting temperature (Tm)
  * GC content
  * Secondary structure formation
* Detect primer-dimer formation
* Evaluate cross-reactivity
* Score and rank candidate primer sets
* Batch design for multiple sRNAs
* Custom stem-loop sequence support

### Optimization Parameters

* Adjustable Tm thresholds
* GC content constraints
* Amplicon length filters
* Self-dimer and hetero-dimer penalties
* Customizable scoring weights

### Graphical User Interface

The software includes a user-friendly GUI to facilitate:

* Batch primer submission
* Parameter configuration
* Real-time result visualization
* Export of optimized primer sets

---

## Running the Primer Designer

```bash
python App.py
```

Ensure that all required dependencies are installed prior to execution (see Installation section below).

---

# 2. Machine Learning Components

The machine learning module provides a reproducible framework for model training, evaluation, and prediction using OM-sqPCR-derived features.

The pipeline is implemented using the scikit-learn framework to ensure modularity and reproducibility.

---

## Data Processing Workflow

* Raw feature extraction from OM-sqPCR output
* Feature filtering and selection
* Z-score normalization
* Pipeline-based model integration
* Strict separation of training and test data
* Nested cross-validation within the training set

Nested cross-validation ensures unbiased estimation of generalization performance.

---

## Implemented Models

The following classifiers are implemented:

* Linear Discriminant Analysis (SVD)
* Linear Discriminant Analysis (LSQR)
* SGD Logistic Regression
* Logistic Regression (ElasticNet)
* Support Vector Machine (RBF kernel)
* Random Forest
* XGBoost
* LightGBM
* Gaussian Naive Bayes
* K-Nearest Neighbors

External gradient boosting libraries:

* XGBoost
* LightGBM

---

## Model Outputs

* Probability scores for binary classification
* Performance metrics (AUC, accuracy, sensitivity, specificity)
* Independent test set evaluation

---

## Installation

```bash
git clone https://github.com/RaceyYe/OM-sqPCR.git
cd OM-sqPCR
pip install -r requirements.txt
```

Dependencies include:

* Python ≥ 3.8
* numpy
* pandas
* scikit-learn
* xgboost
* lightgbm
* matplotlib
* seaborn

---

## Reproducibility and Experimental Design

* Feature scaling is performed inside the pipeline to prevent data leakage.
* Hyperparameter tuning is conducted exclusively within nested cross-validation loops.
* Independent test sets are never used during model selection.

This ensures rigorous statistical validity and publication-grade reproducibility.

---

## Example Workflow

### Primer Design

```bash
python App.py
```

### Machine Learning Training (example)

```bash
python train_model.py --model svm --cv nested
```

---

## Citation

```
[manuscript citation]
```

---

## License

MIT licience

---

## Contact

Maintainer: Racey Ye
GitHub: [https://github.com/RaceyYe](https://github.com/RaceyYe)

---
