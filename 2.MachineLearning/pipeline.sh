#!/bin/bash
# env: rapids-25.04

# basic settings
OUT_DIR="./Outputs"
SCRIPT_PATH="./Modules"
DATA_PATH="./Data"
VIS_PATH=$OUT_DIR/vis
features="has-miR-223-3p,rsRNA-28s-12,tsRNA-Asp-GTC-1,has-let-7f-5p"
PRED_DIR=$OUT_DIR/predictions

#############################################################
# 1. Train models and Predict on hold-out sets
#############################################################
# run the nested cross-validation part for training on the training set
python $SCRIPT_PATH/ML_pipeline.py
# average the repeated OOF into single format
python $SCRIPT_PATH/output_process/transform_repeated_OOF.py $PRED_DIR
# Youden threshold extration
python $SCRIPT_PATH/output_process/Youden_threshold_extraction.py $OUT_DIR/calibration/calibration_summary.csv $PRED_DIR/Youden_thresholds.csv


#############################################################
# 2. Characterize the models'performance, and use statistics for comparison
#############################################################
# auroc curve for nested cv oof predictions
python $SCRIPT_PATH/model_selection/roc_curve_generator.py -i $PRED_DIR/nested_cv_oof_predictions_averaged.xlsx \
 -o $VIS_PATH/auROC/nestedCV -l "TrueLabel" -s "PredProb" --combine --title-suffix "training set" --no-grid

# statistical model comparison tests
mkdir -p $OUT_DIR/statistical_comparison
python $SCRIPT_PATH/model_selection/tests.py $PRED_DIR/nested_cv_oof_predictions_averaged.xlsx $PRED_DIR/Youden_thresholds.csv $OUT_DIR/statistical_comparison

#############################################################
# 3. SHAP analysis for SVM_Linear model
#############################################################
mkdir -p $VIS_PATH/SHAP $VIS_PATH/confusionMatrix
python $SCRIPT_PATH/output_process/SHAP_SVMlinear_explainer.py $OUT_DIR

#############################################################
# 4. Confusion matrix for SVM_Linear model
#############################################################
python $SCRIPT_PATH/output_process/confusionMatrix.py -i $PRED_DIR/internal_testing_predictions.xlsx \
 -o $OUT_DIR/vis/confusionMatrix -t 0.4603 -s "SVM_Linear" -l "non-CRC" "CRC" --title "SVM_linear on internal testing set" --filename "confusion_matrix_internal_testing.pdf"
python $SCRIPT_PATH/output_process/confusionMatrix.py -i $PRED_DIR/external_validation_predictions.xlsx \
 -o $OUT_DIR/vis/confusionMatrix -t 0.4603 -s "SVM_Linear" -l "non-CRC" "CRC" --title "SVM_linear on external validation set" --filename "confusion_matrix_external_validation.pdf"
for cohort in "NJ" "LYG" "ZJ"
do
    python $SCRIPT_PATH/output_process/confusionMatrix.py -i $PRED_DIR/${cohort}_predictions.xlsx \
    -o $OUT_DIR/vis/confusionMatrix -t 0.4603 -s "SVM_Linear" -l "non-CRC" "CRC" --title "SVM_linear on ${cohort} set" --filename "confusion_matrix_${cohort}.pdf"
done

#############################################################
# 5. metrics for models
#############################################################
mkdir -p $PRED_DIR/Subgroup_Results
python $SCRIPT_PATH/output_process/evaluate_model.py \
--model $OUT_DIR/models/SVM_Linear_final.pkl \
--csvs $$DATA_PATH/EV_center/EV_NJ.csv $$DATA_PATH/EV_center/EV_LYG.csv $$DATA_PATH/EV_center/EV_ZJ.csv \
 $$DATA_PATH/Pathological/General/all_test_All.csv $$DATA_PATH/Age/all_test_Elder.csv $$DATA_PATH/Age/all_test_Medium.csv \
 $$DATA_PATH/Age/all_test_Youth.csv $$DATA_PATH/Gender/all_test_Female.csv $$DATA_PATH/Gender/all_test_Male.csv \
 $$DATA_PATH/Pathological/CRP_vs_CRC/all_test_CP.csv  $$DATA_PATH/Pathological/HC_vs_CRC/all_test_HC.csv $$DATA_PATH/Staging/all_test_Staging.csv \
--output $PRED_DIR/Subgroup_Results/results_summary.xlsx