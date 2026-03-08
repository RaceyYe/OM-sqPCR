#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Publication-Quality ROC Curve Generator (Enhanced Version)
============================================================

Features:
- Nature journal style (Arial font, 300 DPI)
- Multiple models comparison on single plot
- Excel input with multiple sheets (training, test, validation)
- **NEW: Feature comparison mode** (compare model vs individual features)
- Smart ranking by performance
- Customizable colors, annotations, grid, fonts
- Square ROC plot with perfect aspect ratio
- Separate PDFs per dataset
- No bootstrap analysis (already done in separate script)
"""

import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from matplotlib import rcParams
from sklearn.metrics import roc_curve, auc
import argparse
import warnings
import os
import sys

warnings.filterwarnings('ignore')

# Set publication-quality defaults (Nature journal style)
rcParams['font.family'] = 'Arial'
rcParams['font.size'] = 8
rcParams['axes.linewidth'] = 0.5
rcParams['xtick.major.width'] = 0.5
rcParams['ytick.major.width'] = 0.5
rcParams['xtick.major.size'] = 3
rcParams['ytick.major.size'] = 3
rcParams['figure.dpi'] = 300
rcParams['savefig.dpi'] = 300
rcParams['savefig.format'] = 'pdf'
rcParams['pdf.fonttype'] = 42  # TrueType fonts for editability


class ROCCurvePlotter:
    """Generate publication-quality ROC curves without bootstrap."""
    
    def __init__(self, colors=None, figsize=(3.5, 3.5),  # Nature single column width
                 show_grid=True, show_diagonal=True, 
                 line_width=1.5, legend_fontsize=6):
        """
        Initialize ROC curve plotter.
        
        Parameters:
        -----------
        colors : list, optional
            List of colors for different models
        figsize : tuple
            Figure size in inches (default: 3.5x3.5 for Nature single column)
        show_grid : bool
            Show grid on plot
        show_diagonal : bool
            Show chance diagonal line
        line_width : float
            Width of ROC curves
        legend_fontsize : float
            Font size for legend
        """
        self.colors = colors or [
            "#ffc400", '#ffa600', '#ff8531', '#ff6361', 
            '#de5a79', '#bc5090', '#58508d', '#2c4875', 
            '#003f5c', "#021720"
        ]
        self.figsize = figsize
        self.show_grid = show_grid
        self.show_diagonal = show_diagonal
        self.line_width = line_width
        self.legend_fontsize = legend_fontsize
    
    def calculate_roc_curve(self, y_true, y_scores):
        """
        Calculate ROC curve and AUC for a single model.
        
        Parameters:
        -----------
        y_true : array-like
            True binary labels
        y_scores : array-like
            Predicted probabilities or scores
            
        Returns:
        --------
        dict : Dictionary containing ROC curve points and AUC
        """
        fpr, tpr, thresholds = roc_curve(y_true, y_scores)
        roc_auc = auc(fpr, tpr)
        
        return {
            'fpr': fpr,
            'tpr': tpr,
            'thresholds': thresholds,
            'auc': roc_auc
        }
    
    def plot_roc_curves(self, data_dict, dataset_name, rank_by_auc=True,
                       output_path=None, title=None, 
                       include_ci_in_legend=False, ci_dict=None):
        """
        Plot ROC curves for multiple models.
        
        Parameters:
        -----------
        data_dict : dict
            Dictionary with model names as keys and tuples of (y_true, y_scores) as values
        dataset_name : str
            Name of the dataset (e.g., 'Training', 'Test', 'Validation')
        rank_by_auc : bool
            Whether to rank models by AUC in legend
        output_path : str, optional
            Path to save the figure
        title : str, optional
            Custom title for the plot
        include_ci_in_legend : bool
            Include 95% CI in legend (requires ci_dict)
        ci_dict : dict, optional
            Dictionary with model AUC confidence intervals
            
        Returns:
        --------
        fig, ax : matplotlib figure and axes objects
        results : dict with model statistics
        """
        fig, ax = plt.subplots(figsize=self.figsize)
        results = {}
        
        # Calculate ROC curves and AUCs for all models
        for model_name, (y_true, y_scores) in data_dict.items():
            roc_results = self.calculate_roc_curve(y_true, y_scores)
            results[model_name] = roc_results
        
        # Sort models by AUC if requested
        if rank_by_auc:
            sorted_models = sorted(results.items(), 
                                 key=lambda x: x[1]['auc'], 
                                 reverse=True)
        else:
            sorted_models = list(results.items())
        
        # Plot ROC curves
        for idx, (model_name, res) in enumerate(sorted_models):
            color = self.colors[idx % len(self.colors)]
            
            # Create legend label
            if include_ci_in_legend and ci_dict and model_name in ci_dict:
                ci_lower = ci_dict[model_name]['lower']
                ci_upper = ci_dict[model_name]['upper']
                label = f"{model_name} (AUC = {res['auc']:.3f}, 95% CI: {ci_lower:.3f}-{ci_upper:.3f})"
            else:
                label = f"{model_name} (AUC = {res['auc']:.3f})"
            
            # Plot ROC curve
            ax.plot(res['fpr'], res['tpr'], 
                   color=color, lw=self.line_width, label=label)
        
        # Plot chance diagonal
        if self.show_diagonal:
            # ax.plot([0, 1], [0, 1], 'k--', lw=0.5, label='Chance', alpha=0.5)
            ax.plot([0, 1], [0, 1], 'k--', lw=0.5) #, label='Chance', alpha=0.5)


        # Set labels and title
        ax.set_xlabel('False Positive Rate (1 - Specificity)', 
                     fontsize=8, fontweight='bold')
        ax.set_ylabel('True Positive Rate (Sensitivity)', 
                     fontsize=8, fontweight='bold')
        
        if title:
            ax.set_title(title, fontsize=9, fontweight='bold', pad=10)
        else:
            ax.set_title(f'ROC Curves - {dataset_name}', 
                        fontsize=9, fontweight='bold', pad=10)
        
        # Set axis limits and aspect ratio (square plot)
        ax.set_xlim([-0.02, 1.02])
        ax.set_ylim([-0.02, 1.02])
        ax.set_aspect('equal', adjustable='box')
        
        # Grid
        if self.show_grid:
            ax.grid(True, alpha=0.3, linestyle=':', linewidth=0.5)
        
        # Legend
        ax.legend(loc='lower right', frameon=True, fancybox=False,
                 shadow=False, fontsize=self.legend_fontsize,
                 edgecolor='none', framealpha=1)
                #  edgecolor='black', framealpha=1)
                 
        
        # Tick parameters
        ax.tick_params(labelsize=7)
        
        # Tight layout
        plt.tight_layout()
        
        # Save if path provided
        if output_path:
            plt.savefig(output_path, dpi=300, bbox_inches='tight', format='pdf')
            print(f"  ✓ Saved: {output_path}")
        
        return fig, ax, results


def load_excel_data(file_path, sheet_name=None, label_col='y_true', 
                   score_col='y_prob', model_name=None):
    """
    Load data from Excel file.
    
    Parameters:
    -----------
    file_path : str
        Path to Excel file
    sheet_name : str or list, optional
        Name of the sheet(s) to load. If None, load all sheets
    label_col : str or list
        Column name for true labels. If list, first valid column is used
    score_col : str or list
        Column name(s) for predicted scores/probabilities
    model_name : str, optional
        Override model name. If None, uses sheet name
        
    Returns:
    --------
    dict : Dictionary with model names as keys and (y_true, y_scores) tuples as values
    """
    try:
        # Read Excel file
        xls = pd.ExcelFile(file_path)
        
        # Determine which sheets to process
        if sheet_name is None:
            sheet_names = xls.sheet_names
        elif isinstance(sheet_name, str):
            sheet_names = [sheet_name]
        else:
            sheet_names = sheet_name
        
        data_dict = {}
        
        for sheet in sheet_names:
            if sheet not in xls.sheet_names:
                print(f"  ⚠️  Warning: Sheet '{sheet}' not found in Excel file")
                continue
            
            df = pd.read_excel(xls, sheet_name=sheet)
            
            # Determine label column
            if isinstance(label_col, list):
                label_col_found = None
                for col in label_col:
                    if col in df.columns:
                        label_col_found = col
                        break
                if label_col_found is None:
                    print(f"  ⚠️  Warning: No label column found in sheet '{sheet}'")
                    continue
                y_true = df[label_col_found].values
            else:
                if label_col not in df.columns:
                    print(f"  ⚠️  Warning: Label column '{label_col}' not found in sheet '{sheet}'")
                    continue
                y_true = df[label_col].values
            
            # Determine score column
            if isinstance(score_col, list):
                score_col_found = None
                for col in score_col:
                    if col in df.columns:
                        score_col_found = col
                        break
                if score_col_found is None:
                    print(f"  ⚠️  Warning: No score column found in sheet '{sheet}'")
                    continue
                y_scores = df[score_col_found].values
            else:
                if score_col not in df.columns:
                    print(f"  ⚠️  Warning: Score column '{score_col}' not found in sheet '{sheet}'")
                    continue
                y_scores = df[score_col].values
            
            # Clean data (remove NaNs)
            mask = ~(np.isnan(y_true) | np.isnan(y_scores))
            y_true = y_true[mask]
            y_scores = y_scores[mask]
            
            if len(y_true) < 2:
                print(f"  ⚠️  Warning: Insufficient data in sheet '{sheet}' after cleaning")
                continue
            
            # Determine model name
            if model_name:
                model_key = model_name
            else:
                model_key = sheet
            
            data_dict[model_key] = (y_true, y_scores)
            print(f"  ✓ Loaded '{model_key}' from sheet '{sheet}': {len(y_true)} samples")
        
        return data_dict
    
    except Exception as e:
        print(f"  ❌ Error loading Excel file: {str(e)}")
        return {}


def load_model_feature_data(file_path, label_col, score_col, features, feature_only=False):
    """
    **NEW FEATURE MODE**
    Load data for comparing a model against individual features, or features only.
    
    Expected Excel structure:
    - Exactly ONE sheet
    - Columns: [label_col, score_col (model predictions), feature1, feature2, ...]
    - For feature_only mode: [label_col, feature1, feature2, ...]
    
    Parameters:
    -----------
    file_path : str
        Path to Excel file
    label_col : str
        Column name for true labels
    score_col : str
        Column name for model predictions (ignored if feature_only=True)
    features : list of str
        List of feature column names to compare
    feature_only : bool, optional
        If True, only plot features (no model comparison)
        
    Returns:
    --------
    tuple : (data_dict, model_name)
        data_dict: Dictionary with model+features as keys and (y_true, y_scores) as values
        model_name: Name of the sheet (used as model name)
    """
    xls = pd.ExcelFile(file_path)
    
    if len(xls.sheet_names) != 1:
        raise ValueError(
            f"Feature comparison mode requires exactly ONE sheet. "
            f"Found {len(xls.sheet_names)} sheets: {xls.sheet_names}"
        )
    
    sheet_name = xls.sheet_names[0]
    df = pd.read_excel(xls, sheet_name=sheet_name)
    
    # Validate label column
    if label_col not in df.columns:
        raise ValueError(
            f"Label column '{label_col}' not found. "
            f"Available columns: {list(df.columns)}"
        )
    
    y_true = df[label_col].values
    
    # Build data dictionary
    data_dict = {}
    
    # Add model predictions (skip if feature_only mode)
    if not feature_only:
        # Validate model score column
        if score_col not in df.columns:
            raise ValueError(
                f"Model score column '{score_col}' not found. "
                f"Available columns: {list(df.columns)}"
            )
        
        data_dict[sheet_name] = (y_true, df[score_col].values)
        print(f"  ✓ Loaded model '{sheet_name}': {len(y_true)} samples")
    else:
        print(f"  ℹ️  Feature-only mode: skipping model column")
    
    # Add individual features
    for feature in features:
        if feature not in df.columns:
            raise ValueError(
                f"Feature column '{feature}' not found. "
                f"Available columns: {list(df.columns)}"
            )
        
        data_dict[feature] = (y_true, df[feature].values)
        print(f"  ✓ Loaded feature '{feature}': {len(y_true)} samples")
    
    # Clean NaNs from all data
    for key in list(data_dict.keys()):
        y_true_temp, y_scores_temp = data_dict[key]
        mask = ~(np.isnan(y_true_temp) | np.isnan(y_scores_temp))
        data_dict[key] = (y_true_temp[mask], y_scores_temp[mask])
    
    return data_dict, sheet_name


def load_ci_from_csv(ci_file_path):
    """
    Load AUC confidence intervals from CSV file.
    
    Supports two CSV formats:
    1. Format with CI_Lower and CI_Upper columns
    2. Format with formatted CI string in 'Format' column
    
    Parameters:
    -----------
    ci_file_path : str
        Path to CSV file with confidence intervals
        
    Returns:
    --------
    dict : Dictionary with model names as keys and CI dicts as values
           Returns None if file not found or no valid data
    """
    if not os.path.exists(ci_file_path):
        print(f"  ⚠️  Warning: CI file not found: {ci_file_path}")
        return None
    
    try:
        ci_df = pd.read_csv(ci_file_path)
        
        # Initialize empty dictionary
        ci_dict = {}
        
        # Filter for AUROC metrics only
        if 'Metric' in ci_df.columns:
            auroc_df = ci_df[ci_df['Metric'] == 'AUROC']
        else:
            auroc_df = ci_df
        
        if len(auroc_df) == 0:
            print(f"  ⚠️  No AUROC metrics found in CI file")
            return None
        
        # Process each AUROC row
        for _, row in auroc_df.iterrows():
            model_name = str(row['Model']).strip()
            
            # Check for different column name possibilities
            # Format 1: CI_Lower and CI_Upper columns
            if 'CI_Lower' in row and 'CI_Upper' in row:
                ci_dict[model_name] = {
                    'lower': float(row['CI_Lower']),
                    'upper': float(row['CI_Upper'])
                }
            elif 'CI_lower' in row and 'CI_upper' in row:
                ci_dict[model_name] = {
                    'lower': float(row['CI_lower']),
                    'upper': float(row['CI_upper'])
                }
            # Format 2: Formatted string in 'Format' column
            elif 'Format' in row and pd.notna(row['Format']):
                format_str = str(row['Format'])
                if '(' in format_str and ')' in format_str:
                    try:
                        ci_part = format_str.split('(')[1].split(')')[0]
                        lower, upper = map(float, ci_part.split('-'))
                        ci_dict[model_name] = {'lower': lower, 'upper': upper}
                    except:
                        continue
        
        if ci_dict:
            print(f"  ✓ Loaded CI data for {len(ci_dict)} models")
            return ci_dict
        else:
            print(f"  ⚠️  No valid CI data found")
            return None
    
    except Exception as e:
        print(f"  ❌ Error loading CI data: {str(e)}")
        return None


def generate_roc_for_all_sheets(input_file, output_dir, label_col='y_true',
                                score_col='y_prob', combine_all=False,
                                ci_file=None, colors=None, title_suffix=None,
                                no_grid=False, no_diagonal=False):
    """
    Generate ROC curves for all sheets in an Excel file.
    
    Parameters:
    -----------
    input_file : str
        Path to input Excel file
    output_dir : str
        Directory to save output PDFs
    label_col : str or list
        Column name for true labels
    score_col : str or list
        Column name for predicted scores
    combine_all : bool
        If True, combine all sheets into one ROC plot
    ci_file : str, optional
        Path to CSV file with confidence intervals for legend
    colors : list, optional
        Custom color palette
    title_suffix : str, optional
        Suffix to add to plot titles
    no_grid : bool
        Disable grid on plots
    no_diagonal : bool
        Disable chance diagonal line
        
    Returns:
    --------
    list : Paths to generated PDF files
    """
    # Create output directory
    os.makedirs(output_dir, exist_ok=True)
    
    print(f"\n📊 Processing file: {os.path.basename(input_file)}")
    print(f"📁 Output directory: {output_dir}")
    
    # Load CI data if provided
    ci_dict = load_ci_from_csv(ci_file) if ci_file else None
    
    # Load Excel file to get sheet names
    try:
        xls = pd.ExcelFile(input_file)
        sheet_names = xls.sheet_names
        print(f"📄 Found {len(sheet_names)} sheets: {', '.join(sheet_names)}")
    except Exception as e:
        print(f"❌ Error reading Excel file: {str(e)}")
        return []
    
    # Initialize plotter
    plotter = ROCCurvePlotter(
        colors=colors,
        show_grid=not no_grid,
        show_diagonal=not no_diagonal
    )
    
    output_paths = []
    
    if combine_all:
        # Combine all sheets into one ROC plot
        print(f"\n📈 Generating combined ROC plot for all sheets...")
        all_data = {}
        
        for sheet in sheet_names:
            data_dict = load_excel_data(input_file, sheet, label_col, score_col, 
                                       model_name=sheet)
            all_data.update(data_dict)
        
        if all_data:
            output_path = os.path.join(output_dir, "all_models_roc.pdf")
            title = "ROC Curves - All Models"
            if title_suffix:
                title += f" ({title_suffix})"
            
            fig, ax, results = plotter.plot_roc_curves(
                all_data,
                dataset_name="All Models",
                rank_by_auc=True,
                output_path=output_path,
                title=title,
                include_ci_in_legend=ci_dict is not None,
                ci_dict=ci_dict
            )
            plt.close(fig)
            output_paths.append(output_path)
            
            # Print summary
            print(f"\n📋 Combined ROC Summary:")
            print("-" * 40)
            for model_name, res in sorted(results.items(), 
                                         key=lambda x: x[1]['auc'], 
                                         reverse=True):
                print(f"  {model_name}: AUC = {res['auc']:.4f}")
    
    else:
        # Generate separate ROC plot for each sheet
        for sheet in sheet_names:
            print(f"\n📈 Processing sheet: '{sheet}'...")
            
            data_dict = load_excel_data(input_file, sheet, label_col, score_col)
            
            if not data_dict:
                print(f"  ⚠️  No valid data found in sheet '{sheet}'")
                continue
            
            # Create output filename
            safe_sheet_name = sheet.replace(' ', '_').replace('/', '_').replace('\\', '_')
            output_path = os.path.join(output_dir, f"{safe_sheet_name}_roc.pdf")
            
            title = f"ROC Curves - {sheet}"
            if title_suffix:
                title += f" ({title_suffix})"
            
            fig, ax, results = plotter.plot_roc_curves(
                data_dict,
                dataset_name=sheet,
                rank_by_auc=True, 
                output_path=output_path,
                title=title,
                include_ci_in_legend=ci_dict is not None,
                ci_dict=ci_dict
            )
            plt.close(fig)
            output_paths.append(output_path)
            
            # Print summary for this sheet
            print(f"  📋 ROC Summary for '{sheet}':")
            for model_name, res in sorted(results.items(), 
                                         key=lambda x: x[1]['auc'], 
                                         reverse=True):
                print(f"    {model_name}: AUC = {res['auc']:.4f}")
    
    return output_paths


def generate_feature_comparison_roc(input_file, output_dir, label_col, 
                                   score_col, features, ci_file=None, 
                                   colors=None, title_suffix=None,
                                   no_grid=False, no_diagonal=False,
                                   feature_only=False):
    """
    **NEW FEATURE MODE**
    Generate ROC curve comparing model vs individual features, or features only.
    
    Parameters:
    -----------
    input_file : str
        Path to input Excel file (must have exactly ONE sheet)
    output_dir : str
        Directory to save output PDF
    label_col : str
        Column name for true labels
    score_col : str
        Column name for model predictions (ignored if feature_only=True)
    features : list of str
        List of feature column names to compare
    ci_file : str, optional
        Path to CSV file with confidence intervals
    colors : list, optional
        Custom color palette (first color for model, rest for features)
    title_suffix : str, optional
        Suffix to add to plot title
    no_grid : bool
        Disable grid on plot
    no_diagonal : bool
        Disable chance diagonal line
    feature_only : bool, optional
        If True, only plot features (no model comparison)
        
    Returns:
    --------
    str : Path to generated PDF file
    """
    # Create output directory
    os.makedirs(output_dir, exist_ok=True)
    
    if feature_only:
        print(f"\n📊 Feature-Only Comparison Mode")
    else:
        print(f"\n📊 Feature Comparison Mode")
    
    print(f"📁 Input file: {os.path.basename(input_file)}")
    print(f"📁 Output directory: {output_dir}")
    print(f"🔍 Features to compare: {', '.join(features)}")
    
    # Load CI data if provided
    ci_dict = load_ci_from_csv(ci_file) if ci_file else None
    
    # Load data in feature comparison mode
    try:
        data_dict, model_name = load_model_feature_data(
            input_file, label_col, score_col, features, feature_only=feature_only
        )
    except Exception as e:
        print(f"❌ Error loading feature comparison data: {str(e)}")
        return None
    
    # Use specialized color palette for feature comparison
    if colors is None:
        if feature_only:
            # Feature-only mode: use diverse palette (no need for model gold)
            colors = [
                "#8C1A34",  # Feature 1 (dark red)
                "#295677",  # Feature 2 (dark blue)
                "#CE814C",  # Feature 3 (orange/brown)
                "#448A66",  # Feature 4 (teal/green)
                "#9B59B6",  # Feature 5 (purple)
                "#E74C3C",  # Feature 6 (red)
                "#3498DB",  # Feature 7 (blue)
                "#F39C12",  # Feature 8 (orange)
                "#16A085",  # Feature 9 (turquoise)
                "#C0392B",  # Feature 10 (dark red)
            ]
        else:
            # Model + features mode: gold for model, diverse for features
            colors = [
                "#FFC300",  # Model (bright yellow/gold)
                "#8C1A34",  # Feature 1 (dark red)
                "#295677",  # Feature 2 (dark blue)
                "#CE814C",  # Feature 3 (orange/brown)
                "#448A66",  # Feature 4 (teal/green)
                "#9B59B6",  # Feature 5 (purple)
                "#E74C3C",  # Feature 6 (red)
                "#3498DB",  # Feature 7 (blue)
            ]
    
    # Initialize plotter
    plotter = ROCCurvePlotter(
        colors=colors,
        show_grid=not no_grid,
        show_diagonal=not no_diagonal
    )
    
    # Create output filename
    if feature_only:
        output_path = os.path.join(output_dir, f"{model_name}_features_only_roc.pdf")
    else:
        output_path = os.path.join(output_dir, f"{model_name}_feature_roc.pdf")
    
    # Create title
    if feature_only:
        title = f"ROC Curves - Features Only"
    else:
        title = f"ROC Curves - {model_name}"
    
    if title_suffix:
        title += f" ({title_suffix})"
    
    # Plot
    if feature_only:
        print(f"\n📈 Generating feature-only ROC plot...")
    else:
        print(f"\n📈 Generating feature comparison ROC plot...")
    
    fig, ax, results = plotter.plot_roc_curves(
        data_dict,
        dataset_name=model_name if not feature_only else "Features",
        rank_by_auc=False, # rank in fixed order
        output_path=output_path,
        title=title,
        include_ci_in_legend=ci_dict is not None,
        ci_dict=ci_dict
    )
    plt.close(fig)
    
    # Print summary
    if feature_only:
        print(f"\n📋 Feature-Only ROC Summary:")
    else:
        print(f"\n📋 Feature ROC Summary:")
    
    print("-" * 40)
    for name, res in sorted(results.items(), key=lambda x: x[1]['auc'], reverse=True):
        auc_val = res['auc']
        if not feature_only and name == model_name:
            marker = "🏆 "
        else:
            marker = "   "
        print(f"{marker}{name}: AUC = {auc_val:.4f}")
    
    return output_path


def main():
    """Main function with argparse."""
    parser = argparse.ArgumentParser(
        description='Generate publication-quality ROC curves from Excel files',
        formatter_class=argparse.RawDescriptionHelpFormatter,
        epilog="""
Examples:
  
  STANDARD MODE - Multi-sheet comparison:
  ----------------------------------------
  # Generate separate ROC plots for each sheet
  python roc_plotter_enhanced.py -i predictions.xlsx -o ./output
  
  # Combine all sheets into one ROC plot
  python roc_plotter_enhanced.py -i predictions.xlsx -o ./output --combine
  
  # Use specific label and score columns
  python roc_plotter_enhanced.py -i data.xlsx -o ./plots -l TrueLabel -s PredProb
  
  # Add confidence intervals from CSV
  python roc_plotter_enhanced.py -i models.xlsx -o ./results --ci-file ./ci_results.csv
  
  FEATURE COMPARISON MODE - Model vs Features:
  ---------------------------------------------
  # Compare model against individual features
  python roc_plotter_enhanced.py -i model_data.xlsx -o ./output \\
      --model-features --features Age,BMI,Glucose
  
  # Feature comparison with CI from bootstrap analysis
  python roc_plotter_enhanced.py -i model_data.xlsx -o ./output \\
      --model-features --features Age,BMI,Glucose \\
      --ci-file ./bootstrap_ci.csv
  
  # Feature comparison with custom label/score columns
  python roc_plotter_enhanced.py -i data.xlsx -o ./output \\
      --model-features -l TrueLabel -s ModelProb \\
      --features Feature1,Feature2,Feature3
  
  FEATURE-ONLY MODE - Compare Features Only:
  ------------------------------------------
  # Plot only features without model comparison
  python roc_plotter_enhanced.py -i features.xlsx -o ./output \\
      --feature-only --features Age,BMI,Glucose,Cholesterol
  
  # Feature-only with confidence intervals
  python roc_plotter_enhanced.py -i features.xlsx -o ./output \\
      --feature-only --features Age,BMI,Glucose \\
      --ci-file ./feature_ci.csv
  
  # Feature-only with custom styling
  python roc_plotter_enhanced.py -i features.xlsx -o ./output \\
      --feature-only --features Age,BMI,Glucose \\
      --no-grid --colors "#8C1A34,#295677,#CE814C"
"""
    )
    
    # Required arguments
    parser.add_argument('-i', '--input', type=str, required=True,
                       help='Input Excel file path')
    parser.add_argument('-o', '--output', type=str, required=True,
                       help='Output directory for PDF files')
    
    # Optional arguments
    parser.add_argument('-l', '--label-col', type=str, default='y_true',
                       help='Column name for true labels (default: y_true). '
                            'Can also be a comma-separated list')
    parser.add_argument('-s', '--score-col', type=str, default='y_prob',
                       help='Column name for predicted scores (default: y_prob). '
                            'Can also be a comma-separated list')
    
    parser.add_argument('--sheets', type=str, nargs='+',
                       help='Specific sheet names to process (default: all sheets)')
    parser.add_argument('--combine', action='store_true',
                       help='Combine all sheets into one ROC plot')
    
    # **NEW: Feature comparison mode**
    parser.add_argument('--model-features', action='store_true',
                       help='**FEATURE MODE**: Compare model vs individual features. '
                            'Requires --features argument. Excel must have ONE sheet.')
    parser.add_argument('--feature-only', action='store_true',
                       help='**FEATURE-ONLY MODE**: Plot only features (no model). '
                            'Requires --features argument. Excel must have ONE sheet.')
    parser.add_argument('--features', type=str,
                       help='Comma-separated list of feature column names for comparison. '
                            'Required when --model-features or --feature-only is used.')
    
    # Styling arguments
    parser.add_argument('--ci-file', type=str,
                       help='CSV file with AUC confidence intervals for legend')
    parser.add_argument('--title-suffix', type=str,
                       help='Suffix to add to plot titles')
    parser.add_argument('--colors', type=str,
                       help='Comma-separated list of hex colors (e.g., #ff0000,#00ff00)')
    parser.add_argument('--no-grid', action='store_true',
                       help='Disable grid on plots')
    parser.add_argument('--no-diagonal', action='store_true',
                       help='Disable chance diagonal line')
    
    args = parser.parse_args()
    
    print("=" * 70)
    print("Publication-Quality ROC Curve Generator (Enhanced)")
    print("=" * 70)
    print()
    
    # Validate feature mode
    if args.model_features and args.feature_only:
        parser.error("Cannot use both --model-features and --feature-only. Choose one.")
    
    if (args.model_features or args.feature_only) and not args.features:
        parser.error("--features is required when using --model-features or --feature-only")
    
    # Process colors argument
    if args.colors:
        colors = [color.strip() for color in args.colors.split(',')]
    else:
        colors = None
    
    # **FEATURE COMPARISON MODE or FEATURE-ONLY MODE**
    if args.model_features or args.feature_only:
        if args.feature_only:
            print("🔬 MODE: Feature-Only Comparison")
        else:
            print("🔬 MODE: Feature Comparison")
        print("-" * 70)
        
        features = [f.strip() for f in args.features.split(',')]
        
        output_path = generate_feature_comparison_roc(
            input_file=args.input,
            output_dir=args.output,
            label_col=args.label_col,
            score_col=args.score_col,
            features=features,
            ci_file=args.ci_file,
            colors=colors,
            title_suffix=args.title_suffix,
            no_grid=args.no_grid,
            no_diagonal=args.no_diagonal,
            feature_only=args.feature_only
        )
        
        # Print summary
        print("\n" + "=" * 70)
        print("SUMMARY")
        print("=" * 70)
        if output_path:
            if args.feature_only:
                print(f"✓ Successfully generated feature-only ROC plot:")
            else:
                print(f"✓ Successfully generated feature comparison ROC plot:")
            print(f"  - {os.path.basename(output_path)}")
            print(f"\n📁 Plot saved to: {os.path.abspath(args.output)}")
        else:
            print("❌ ROC plot generation failed.")
    
    # **STANDARD MULTI-SHEET MODE**
    else:
        print("📊 MODE: Multi-Sheet Comparison")
        print("-" * 70)
        
        # Process label column argument (could be comma-separated list)
        if ',' in args.label_col:
            label_col = [col.strip() for col in args.label_col.split(',')]
        else:
            label_col = args.label_col
        
        # Process score column argument (could be comma-separated list)
        if ',' in args.score_col:
            score_col = [col.strip() for col in args.score_col.split(',')]
        else:
            score_col = args.score_col
        
        # Generate ROC curves
        output_paths = generate_roc_for_all_sheets(
            input_file=args.input,
            output_dir=args.output,
            label_col=label_col,
            score_col=score_col,
            combine_all=args.combine,
            ci_file=args.ci_file,
            colors=colors,
            title_suffix=args.title_suffix,
            no_grid=args.no_grid,
            no_diagonal=args.no_diagonal
        )
        
        # Print summary
        print("\n" + "=" * 70)
        print("SUMMARY")
        print("=" * 70)
        if output_paths:
            print(f"✓ Successfully generated {len(output_paths)} ROC plot(s):")
            for path in output_paths:
                print(f"  - {os.path.basename(path)}")
            print(f"\n📁 All plots saved to: {os.path.abspath(args.output)}")
        else:
            print("❌ No ROC plots were generated. Please check your input data.")
    
    print("\n" + "=" * 70)


if __name__ == "__main__":
    main()