#!/usr/bin/env python3

import os
import argparse
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import matplotlib as mpl
from sklearn.metrics import confusion_matrix, ConfusionMatrixDisplay


def plot_confusion_matrices(cm, cm_normalized, args):
    """
    Plot confusion matrices with PDF compatibility for Adobe Illustrator
    """
    # Set global matplotlib parameters for better PDF compatibility
    mpl.rcParams['pdf.fonttype'] = 42  # Use TrueType fonts (more compatible with Illustrator)
    mpl.rcParams['pdf.use14corefonts'] = True  # Use base 14 PDF fonts
    mpl.rcParams['font.family'] = 'sans-serif'
    mpl.rcParams['font.sans-serif'] = ['Arial', 'Helvetica', 'DejaVu Sans']
    mpl.rcParams['axes.unicode_minus'] = False  # Fix minus sign display
    mpl.rcParams['savefig.dpi'] = 300  # High resolution
    mpl.rcParams['savefig.bbox'] = 'tight'
    
    # Create figure with two subplots
    fig, axes = plt.subplots(1, 2, figsize=(12, 5))
    
    titles = [
        args.title,
        args.title + " (Normalized by True Label)"
    ]
    
    cms = [cm, cm_normalized]
    
    for ax, matrix, title in zip(axes, cms, titles):
        # Format values based on matrix type
        if matrix is cm_normalized:
            # Normalized matrix - show as percentages with 2 decimal places
            disp = ConfusionMatrixDisplay(
                confusion_matrix=matrix,
                display_labels=args.label
            )
            disp.plot(
                ax=ax, 
                cmap=plt.cm.Blues, 
                colorbar=False, 
                values_format='.2f',  # 2 decimal places for normalized
                text_kw={'fontsize': 11, 'fontfamily': 'sans-serif', 'fontweight': 'normal'}
            )
        else:
            # Count matrix - show as integers
            disp = ConfusionMatrixDisplay(
                confusion_matrix=matrix,
                display_labels=args.label
            )
            disp.plot(
                ax=ax, 
                cmap=plt.cm.Blues, 
                colorbar=False, 
                values_format='d',  # integer format for counts
                text_kw={'fontsize': 11, 'fontfamily': 'sans-serif', 'fontweight': 'normal'}
            )
        
        # Customize labels with proper font settings
        ax.set_title(title, fontfamily='sans-serif', fontsize=12, fontweight='bold')
        ax.set_xlabel('Predicted Label', fontfamily='sans-serif', fontsize=11)
        ax.set_ylabel('True Label', fontfamily='sans-serif', fontsize=11)
        
        # Ensure tick labels have proper font
        for label in ax.get_xticklabels() + ax.get_yticklabels():
            label.set_fontfamily('sans-serif')
            label.set_fontsize(10)
    
    plt.tight_layout()
    
    # Save with PDF/A-1b compatibility settings
    output_path = os.path.join(args.output_dir, args.filename)
    plt.savefig(
        output_path,
        format="pdf",
        bbox_inches="tight",
        dpi=300,
        metadata={
            'Creator': 'Python Matplotlib',
            'Producer': 'Matplotlib PDF Backend',
            'Title': args.title,
            'Subject': 'Confusion Matrix',
            'Keywords': 'confusion matrix, classification, machine learning'
        }
    )
    plt.close()
    
    # Also save as EPS as backup (often more compatible with Illustrator)
    eps_path = output_path.replace('.pdf', '.eps')
    try:
        plt.figure(figsize=(12, 5))
        axes_eps = [plt.subplot(1, 2, 1), plt.subplot(1, 2, 2)]
        
        for ax, matrix, title in zip(axes_eps, cms, titles):
            if matrix is cm_normalized:
                disp = ConfusionMatrixDisplay(confusion_matrix=matrix, display_labels=args.label)
                disp.plot(ax=ax, cmap=plt.cm.Blues, colorbar=False, values_format='.2f')
            else:
                disp = ConfusionMatrixDisplay(confusion_matrix=matrix, display_labels=args.label)
                disp.plot(ax=ax, cmap=plt.cm.Blues, colorbar=False, values_format='d')
            ax.set_title(title)
            ax.set_xlabel('Predicted Label')
            ax.set_ylabel('True Label')
        
        plt.tight_layout()
        plt.savefig(eps_path, format='eps', bbox_inches="tight", dpi=300)
        plt.close()
        print(f"EPS backup saved to: {eps_path}")
    except Exception as e:
        print(f"Note: Could not save EPS backup: {e}")
    
    return output_path


def main():
    # =========================
    # Argument parser
    # =========================
    parser = argparse.ArgumentParser(
        description="Draw confusion matrix from Excel predictions file."
    )

    parser.add_argument(
        '-i', "--input",
        type=str,
        required=True,
        help="Path to Excel file containing predictions"
    )
    parser.add_argument(
        '-s', "--sheet",
        type=str,
        default="Model",
        help="Sheet name in Excel file"
    )
    parser.add_argument(
        '-t', "--threshold",
        type=float,
        required=True,
        help="Probability threshold for classification"
    )
    parser.add_argument(
        '-o', "--output_dir",
        type=str,
        default=".",
        help="Output directory for confusion matrix plot"
    )
    parser.add_argument(
        '-l', "--label",
        type=str,
        nargs=2,
        default=["Negative", "Positive"],
        help="Labels for the two classes (e.g., 'Non-CRC' 'CRC')"
    )
    parser.add_argument(
        "--title",
        type=str,
        default="Confusion Matrix",
        help="Title for confusion matrix plot"
    )
    parser.add_argument(
        "--filename",
        type=str,
        default="confusion_matrix.pdf",
        help="Filename for confusion matrix plot (PDF)"
    )

    args = parser.parse_args()

    # =========================
    # Load Excel file
    # =========================
    print(f"Loading data from {args.input}, sheet: {args.sheet}")
    df = pd.read_excel(args.input, sheet_name=args.sheet)

    required_columns = {"SampleID", "TrueLabel", "PredProb"}
    if not required_columns.issubset(df.columns):
        raise ValueError(
            f"Input file must contain columns: {required_columns}"
        )

    # =========================
    # Generate binary predictions
    # =========================
    y_true = df["TrueLabel"].values
    y_prob = df["PredProb"].values

    y_pred = (y_prob >= args.threshold).astype(int)

    # =========================
    # Compute confusion matrix
    # =========================
    cm = confusion_matrix(y_true, y_pred)
    cm_normalized = confusion_matrix(y_true, y_pred, normalize="true")

    np.set_printoptions(precision=2)

    print("\n" + "="*50)
    print("Confusion Matrix Results")
    print("="*50)
    print("Confusion matrix (counts):")
    print(cm)

    print("\nNormalized confusion matrix (by true labels):")
    print(cm_normalized)

    # =========================
    # Calculate metrics
    # =========================
    tn = cm[0, 0]  # True negatives
    fp = cm[0, 1]  # False positives
    fn = cm[1, 0]  # False negatives
    tp = cm[1, 1]  # True positives

    total = tn + fp + fn + tp
    accuracy = (tp + tn) / total
    specificity = tn / (tn + fp) if (tn + fp) > 0 else 0  # True Negative Rate
    sensitivity = tp / (tp + fn) if (tp + fn) > 0 else 0  # True Positive Rate (Recall)
    precision = tp / (tp + fp) if (tp + fp) > 0 else 0
    f1_score = 2 * (precision * sensitivity) / (precision + sensitivity) if (precision + sensitivity) > 0 else 0

    print("\n" + "="*50)
    print("Performance Metrics")
    print("="*50)
    print(f"Threshold: {args.threshold:.3f}")
    print(f"Total samples: {total}")
    print(f"Accuracy: {accuracy:.4f}")
    print(f"Specificity (True Negative Rate): {specificity:.4f}")
    print(f"Sensitivity (Recall): {sensitivity:.4f}")
    print(f"Precision: {precision:.4f}")
    print(f"F1 Score: {f1_score:.4f}")
    print(f"True Negatives: {tn}")
    print(f"False Positives: {fp}")
    print(f"False Negatives: {fn}")
    print(f"True Positives: {tp}")

    # =========================
    # Plot confusion matrices
    # =========================
    print("\n" + "="*50)
    print("Generating plots...")
    print("="*50)
    
    # Create output directory if it doesn't exist
    os.makedirs(args.output_dir, exist_ok=True)
    
    output_path = plot_confusion_matrices(cm, cm_normalized, args)

    print(f"\n✓ Confusion matrix plot saved to: {output_path}")
    
    # Provide Illustrator compatibility tips
    print("\n" + "="*50)
    print("Adobe Illustrator Compatibility Tips")
    print("="*50)
    print("If the PDF still doesn't open correctly in Illustrator:")
    print("1. Try opening the EPS backup file instead")
    print("2. Open the PDF in Adobe Acrobat and save as 'Illustrator compatible PDF'")
    print("3. In Illustrator, try: File → Open and select the PDF, then choose 'Import PDF' options")
    print("4. Use Illustrator's 'Place' command instead of 'Open'")
    print("5. Update Illustrator to the latest version")


if __name__ == "__main__":
    main()