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
|
├── 1.OM-SqPCR_primerDesigner
│   ├── App.py
│   ├── integrated_execution.py
│   ├── Modules/
│   └── outputs_integration.py
├── 2.MachineLearning
│   ├── Modules/
│   └── pipeline.sh
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

* Pipeline-based model integration
* Z-score normalization
* Strict separation of training and test data
* Nested cross-validation within the training set

10 repeatations of nested cross-validation ensure unbiased estimation of generalization performance.

---

## Implemented Models

The following classifiers are implemented:

* Linear Discriminant Analysis (LSQR)
* Logistic Regression (L2)
* Support Vector Machine (Linear kernel)
* Random Forest

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
* matplotlib

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
bash pipeline.sh
```

---

## Citation

```
[manuscript citation]
```

---

## License

MIT license

---

## Contact

Maintainer: Racey Ye
GitHub: [https://github.com/RaceyYe](https://github.com/RaceyYe)

---
