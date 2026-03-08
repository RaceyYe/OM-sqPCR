import pandas as pd
import os
import sys

if __name__ == "__main__":
    # Take ONLY the first argument as base_dir
    if len(sys.argv) < 2:
        print("❌ Error: Please provide base directory path")
        print("Usage: python transform_repeated_OOF.py /path/to/directory")
        sys.exit(1)
    
    base_dir = sys.argv[1]  # This is a string, not a list
    
    file_path = os.path.join(base_dir, "nested_cv_oof_predictions.xlsx")
    output_path = os.path.join(base_dir, 'nested_cv_oof_predictions_averaged.xlsx')

    # Rest of your code...
    sheet_names = pd.ExcelFile(file_path).sheet_names
    print(f"Processing {len(sheet_names)} sheets: {sheet_names}")

    with pd.ExcelWriter(output_path) as writer:
        for sheet in sheet_names:
            df = pd.read_excel(file_path, sheet_name=sheet)
            result_df = (df.groupby(['SampleID', 'TrueLabel'])['PredProb']
                        .mean()
                        .round(4)
                        .reset_index())
            result_df.to_excel(writer, sheet_name=sheet, index=False)
            print(f"✓ {sheet}: {len(df)} → {len(result_df)} rows")

    print(f"\n✅ All sheets saved to: {output_path}")