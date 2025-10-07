import pandas as pd
import numpy as np
import multiprocessing as mp
from functools import partial, lru_cache
import os
from modules.auxiliary_modules.SecStructure import analyze_primers_interaction
# from modules.auxiliary_modules.check import RatingCheck


@lru_cache(maxsize=1024)
def cached_primer_interaction(fp, rp, temperature):
    """Cached version of primer interaction analysis"""
    df = analyze_primers_interaction(fp, rp)
    return df[f'dG_HET_{temperature}C'].iloc[0]


def add_rating_column(df, temp, primer_type):
    """Add rating scores based on ΔG values"""
    dg_hd_col = f'dG_HD_{temp}C_{primer_type}'
    dg_hp_col = f'dG_HP_{temp}C_{primer_type}'

    if dg_hd_col in df.columns and dg_hp_col in df.columns:
        df[f'Rating_{primer_type}'] = 100 + 1.8 * df[dg_hd_col] + 1.4 * df[dg_hp_col]

    return df


# def process_single_file(file_path, temp, primer_type):
#     """Process a single primer file with batch processing"""
#     print(f"Processing {file_path}...")
#
#     if not os.path.exists(file_path):
#         print(f"File not found: {file_path}")
#         return False
#
#     # Read the file
#     df = pd.read_csv(file_path)
#
#     # Process in batches to reduce memory pressure
#     batch_size = 5000
#     updated_dfs = []
#
#     for i in range(0, len(df), batch_size):
#         batch = df.iloc[i:i + batch_size].copy()
#         checker = RatingCheck(batch, temp=temp)
#         updated_batch = checker.run_from_dataframe(batch)
#
#         # Add rating column if needed
#         if f'Rating_{primer_type}' not in updated_batch.columns:
#             updated_batch = add_rating_column(updated_batch, temp, primer_type)
#
#         updated_dfs.append(updated_batch)
#
#     # Combine batches and save
#     final_df = pd.concat(updated_dfs)
#     final_df.to_csv(file_path, index=False)
#
#     print(f"Updated {file_path} with rating scores")
#     return True


# def process_primer_files(directory, temp):
#     """Process primer files using multiprocessing"""
#     fp_file = os.path.join(directory, "designed_FP.csv")
#     rp_file = os.path.join(directory, "designed_RP_PC.csv")
#
#     # Verify files exist
#     files_to_process = []
#     if os.path.exists(fp_file):
#         files_to_process.append((fp_file, temp, "FP"))
#     else:
#         print(f"Warning: Forward primer file not found at {fp_file}")
#
#     if os.path.exists(rp_file):
#         files_to_process.append((rp_file, temp, "RP"))
#     else:
#         print(f"Warning: Reverse primer file not found at {rp_file}")
#
#     if not files_to_process:
#         print("Error: No primer files found for processing!")
#         return False
#
#     # Use multiprocessing for file processing
#     with mp.Pool(processes=min(2, mp.cpu_count())) as pool:
#         results = [pool.apply_async(process_single_file, args=args) for args in files_to_process]
#
#         # Get results
#         all_success = all(result.get() for result in results)
#
#     if all_success:
#         print("Preprocessing complete!")
#         return True
#     else:
#         print("Preprocessing encountered issues.")
#         return False

def process_primer_pair(args):
    """Process a single forward-reverse primer pair"""
    (forward_row, reverse_row,
     delta_G_min, delta_G_max,
     delta_G_FR_min, delta_G_FR_max,
     delta_G_PO_min, delta_G_PO_max,
     temperature) = args

    try:
        # Extract primer sequences and properties
        fp = forward_row["Forward_Primer"]
        rp = reverse_row["RP_Sequence"]
        pc = reverse_row["Probe_Sequence"]

        tm_fp = float(forward_row["Tm_FP"])
        tm_rp = float(reverse_row["Tm_RP"])
        tm_pc = float(reverse_row["Tm_Probe"])

        fp_woe = float(forward_row["FP_WOE"])
        gap = float(reverse_row["Gap"])
        rp_start = float(reverse_row["RP_Start"])

        # if fp_woe + (-1)*gap + (-1.5)*rp_start < 0:
        #     return None

        # Skip if Tm difference is too large
        if abs(tm_fp - tm_rp) > 1:
            return None

        # Helper function to safely convert any value to float
        def safe_float(val):
            if isinstance(val, tuple):
                if len(val) > 0:
                    try:
                        return float(val[0])
                    except (ValueError, TypeError):
                        return None
                return None
            try:
                return float(val)
            except (ValueError, TypeError):
                return None

        # Get and convert ΔG values for homodimers and hairpins
        dg_hd_fp = safe_float(forward_row[f'dG_HD_{temperature}C_FP'])
        dg_hp_fp = safe_float(forward_row[f'dG_HP_{temperature}C_FP'])
        dg_hd_rp = safe_float(reverse_row[f'dG_HD_{temperature}C_RP'])
        dg_hp_rp = safe_float(reverse_row[f'dG_HP_{temperature}C_RP'])
        dg_hd_pc = safe_float(reverse_row[f'dG_HD_{temperature}C_Probe'])
        dg_hp_pc = safe_float(reverse_row[f'dG_HP_{temperature}C_Probe'])

        # Check if any value is None after conversion
        self_values = [dg_hd_fp, dg_hp_fp, dg_hd_rp, dg_hp_rp, dg_hd_pc, dg_hp_pc]
        if None in self_values:
            return None

        # Check self-structure values
        if not all(delta_G_min <= val <= delta_G_max for val in self_values):
            return None

        # Get and convert heterodimer values
        dg_rp_pc = safe_float(reverse_row[f"dG_HET_{temperature}C_HET"])
        dg_fp_rp = safe_float(cached_primer_interaction(fp, rp, temperature))
        dg_fp_pc = safe_float(cached_primer_interaction(fp, pc, temperature))

        # Check if any heterodimer value is None after conversion
        if None in [dg_rp_pc, dg_fp_rp, dg_fp_pc]:
            return None

        # Check interaction values
        if not (delta_G_FR_min <= dg_fp_rp <= delta_G_FR_max and
                delta_G_PO_min <= dg_fp_pc <= delta_G_PO_max and
                delta_G_PO_min <= dg_rp_pc <= delta_G_PO_max):
            return None

        # Return valid primer pair data
        return {
            "Forward_Primer": fp,
            "Reverse_Primer": rp,
            "Probe_Carrier": pc,
            "Tm_FP": tm_fp,
            "Tm_RP": tm_rp,
            "Tm_PC": tm_pc,
            "FP_start_5T": fp_woe,
            "RP_gap": gap,
            "RP_start_5T": rp_start,
            "Tm_diff": abs(tm_rp - tm_fp),
            f'dG_HD_{temperature}C_FP': dg_hd_fp,
            f'dG_HD_{temperature}C_RP': dg_hd_rp,
            f'dG_HD_{temperature}C_PC': dg_hd_pc,
            f'dG_HP_{temperature}C_FP': dg_hp_fp,
            f'dG_HP_{temperature}C_RP': dg_hp_rp,
            f'dG_HP_{temperature}C_PC': dg_hp_pc,
            f"dG_CD_{temperature}C_FP_RP": dg_fp_rp,
            f"dG_CD_{temperature}C_FP_PC": dg_fp_pc,
            f"dG_CD_{temperature}C_RP_PC": dg_rp_pc,
        }
    except Exception as e:
        # Add debugging to see which primer pair is causing issues
        print(f"Error processing primer pair: {str(e)}")
        return None

def integrate_primers(forward_df, reverse_df, **kwargs):
    """Integrate primers with parallel processing"""
    # Set default parameters
    params = {
        'delta_G_min': -4.5,
        'delta_G_max': 4.5,
        'delta_G_FR_min': -5.5,
        'delta_G_FR_max': 5.5,
        'delta_G_PO_min': -7.5,
        'delta_G_PO_max': 7.5,
        'temperature': 37
    }
    params.update(kwargs)

    # Check for required columns
    required_forward_cols = ["Forward_Primer", "Tm_FP", f"dG_HD_{params['temperature']}C_FP",
                             f"dG_HP_{params['temperature']}C_FP"]
    required_reverse_cols = ["RP_Sequence", "Tm_RP", "Probe_Sequence", "Tm_Probe",
                             f"dG_HD_{params['temperature']}C_RP", f"dG_HP_{params['temperature']}C_RP",
                             f"dG_HD_{params['temperature']}C_Probe", f"dG_HP_{params['temperature']}C_Probe",
                             f"dG_HET_{params['temperature']}C_HET"]

    for df, cols, name in [(forward_df, required_forward_cols, "forward"),
                           (reverse_df, required_reverse_cols, "reverse")]:
        missing = [col for col in cols if col not in df.columns]
        if missing:
            raise ValueError(f"Missing required columns in {name} dataframe: {', '.join(missing)}")

    # Create tasks for multiprocessing
    tasks = []
    forbidden_check = "Forbidden_Region" in reverse_df.columns and "FP_WOE" in forward_df.columns

    # Apply forbidden region filter if columns exist
    if forbidden_check:
        print("Applying forbidden region filter...")
        for _, rev_row in reverse_df.iterrows():
            forbidden = rev_row['Forbidden_Region']
            filtered_forward = forward_df[forward_df['FP_WOE'] < forbidden]
            for _, fwd_row in filtered_forward.iterrows():
                tasks.append((fwd_row, rev_row,
                              params['delta_G_min'], params['delta_G_max'],
                              params['delta_G_FR_min'], params['delta_G_FR_max'],
                              params['delta_G_PO_min'], params['delta_G_PO_max'],
                              params['temperature']))
    else:
        print("No forbidden region columns found, processing all combinations...")
        for _, rev_row in reverse_df.iterrows():
            for _, fwd_row in forward_df.iterrows():
                tasks.append((fwd_row, rev_row,
                              params['delta_G_min'], params['delta_G_max'],
                              params['delta_G_FR_min'], params['delta_G_FR_max'],
                              params['delta_G_PO_min'], params['delta_G_PO_max'],
                              params['temperature']))

    print(f"Processing {len(tasks)} potential primer pairs...")

    # Process in parallel
    num_cores = max(1, mp.cpu_count() - 1)  # Leave one core free
    with mp.Pool(processes=num_cores) as pool:
        results = pool.map(process_primer_pair, tasks)

    # Filter None results and convert to DataFrame
    valid_results = [r for r in results if r is not None]
    print(f"Found {len(valid_results)} valid primer pairs.")

    if not valid_results:
        return pd.DataFrame()

    return pd.DataFrame(valid_results)

def sort_primers(integrated_df: pd.DataFrame,
                  optimization_sorting: list[str], task_type = "type2") -> pd.DataFrame:
    """
    Sorts primer pairs based on specified optimization criteria.

    Parameters
    ----------
    integrated_df : pd.DataFrame
        DataFrame containing primer metrics.
    optimization_sorting : list[str]
        Subset and order of ["TM", "DG", "LP"] to define sorting.

    Returns
    -------
    pd.DataFrame
        Sorted DataFrame with additional columns for metrics and rank.
    """
    # 1) TM: absolute distance of mean Tm to 59°C (smaller is better)
    df = integrated_df.copy()

    tm_base_1 = 71.0
    tm_base_2 = 59.0

    if task_type == "type2":
        # 1) TM: absolute distance of each Tm to 59°C (smaller is better)
        df['Key_TM'] = (
                (df['Tm_FP'] - tm_base_2).abs() +
                (df['Tm_RP'] - tm_base_2).abs() +
                (df['Tm_PC'] - tm_base_2).abs()
        )
    else:
        # 2) TM: absolute distance of each Tm to 71°C
        df['Key_TM'] = (
                (df['Tm_FP'] - tm_base_2).abs() +
                (df['Tm_RP'] - tm_base_2).abs() +
                (df['Tm_PC'] - tm_base_1).abs()
        )

    # 2) DG: average delta G across all specified columns (larger is better)
    dg_cols = [
        'dG_HD_37C_FP', 'dG_HD_37C_RP', 'dG_HD_37C_PC',
        'dG_HP_37C_FP', 'dG_HP_37C_RP', 'dG_HP_37C_PC',
        'dG_CD_37C_FP_RP', 'dG_CD_37C_FP_PC', 'dG_CD_37C_RP_PC'
    ]
    df['Avg_dG'] = df[dg_cols].mean(axis=1)
    # for sorting purposes, we want reverse order: higher Avg_dG => smaller key
    df['Key_DG'] = -df['Avg_dG']

    # 3) LP: sum of RP_Start and Gap (smaller is better)
    df['Key_LP'] = df['RP_start_5T'] + df['RP_gap'] - df['FP_start_5T']

    # Build sort keys based on user-specified order
    key_map = {'TM': 'Key_TM', 'DG': 'Key_DG', 'LP': 'Key_LP'}
    sort_keys = [key_map[opt] for opt in optimization_sorting]

    # Define a weight for each optimization option (must align 1:1 with optimization_sorting)
    # e.g. if you want TM×100 (smaller is better), DG×2.0 (smaller is better), LP×0.5 (smaller is better):
    weight_map = {'TM': 0.5, 'DG': 10.0, 'LP': 0.5}
    weights = [weight_map[opt] for opt in optimization_sorting]

    # Compute a weighted score instead of a plain sum:
    #   each column multiplied by its weight, then row-sum
    df['Score'] = (df[sort_keys] * weights).sum(axis=1)
    # alternatively you can do:
    # df['Score'] = df[sort_keys].mul(weights, axis=1).sum(axis=1)

    # Sort by that weighted score
    df_sorted = df.sort_values(by='Score', ascending=True)

    # # Build sort keys based on user-specified order
    # key_map = {'TM': 'Key_TM', 'DG': 'Key_DG', 'LP': 'Key_LP'}
    # sort_keys = [key_map[opt] for opt in optimization_sorting]
    # # Compute the score as the sum of the selected keys
    # df['Score'] = df[sort_keys].sum(axis=1)
    # # Sort the dataframe by the score (ascending by default)
    # df_sorted = df.sort_values(by='Score', ascending=True)

    # # Build sort keys based on user-specified order
    # key_map = {'TM': 'Key_TM', 'DG': 'Key_DG', 'LP': 'Key_LP'}
    # sort_keys = [key_map[opt] for opt in optimization_sorting]
    # # Ascending on Key_TM and Key_LP, but for Key_DG since we negated, also ascending
    # ascending = [True] * len(sort_keys)
    # # Sort
    # df_sorted = df.sort_values(by=sort_keys, ascending=ascending, ignore_index=True)

    # Assign ranks: identical key tuples get same rank (dense ranking)
    tuples = df_sorted[sort_keys].apply(tuple, axis=1)
    df_sorted['Rank'] = pd.factorize(tuples)[0] + 1

    return df_sorted

def interaction_check(directory, criteria=["TM","LP"], **kwargs):
    """Main function to process and integrate primer files"""
    # Extract parameters
    params = {
        'delta_G_min': -4.5,
        'delta_G_max': 4.5,
        'delta_G_FR_min': -5.5,
        'delta_G_FR_max': 5.5,
        'delta_G_PO_min': -7.5,
        'delta_G_PO_max': 7.5,
        'temperature': 37
    }
    params.update(kwargs)

    # Ensure directory exists
    if not os.path.isdir(directory):
        raise ValueError(f"Directory not found: {directory}")

    # # Process primer files
    # if not process_primer_files(directory, params['temperature']):
    #     print("Error preprocessing files. Stopping.")
    #     return

    # Load files
    fp_path = os.path.join(directory, "designed_FP.csv")
    rp_path = os.path.join(directory, "designed_RP_PC.csv")

    if not os.path.exists(fp_path) or not os.path.exists(rp_path):
        print(f"Required files not found in {directory}")
        return

    forward_df = pd.read_csv(fp_path)
    reverse_df = pd.read_csv(rp_path)

    print(f"Loaded {len(forward_df)} forward primers and {len(reverse_df)} reverse primers")

    # Integrate primers
    integrated_df = integrate_primers(forward_df, reverse_df, **params)

    # Save results
    if not integrated_df.empty:
        output_file = os.path.join(directory, "Integrated_df.csv")
        ########## add the checks for the criteria check
        integrated_df = sort_primers(integrated_df,optimization_sorting=criteria)

        integrated_df.to_csv(output_file, index=False)
        print(f"Saved {len(integrated_df)} integrated primer pairs to {output_file}")
    else:
        print("No valid primer pairs found. No output file generated.")


def main():
    """Entry point function with default parameters"""
    # Default parameters
    params = {
        'delta_G_min': -4.5,
        'delta_G_max': 4.5,
        'delta_G_FR_min': -5.5,
        'delta_G_FR_max': 5.5,
        'delta_G_PO_min': -7.5,
        'delta_G_PO_max': 7.5,
        'temperature': 37,
        # 'directory': "/home/racey/Desktop/primerDesign/Outputs/rsRNA"  # Default directory
        'directory': "/home/racey/Desktop/primerDesign/Outputs/Combination_execution/Task1"  # Default directory
    }

    # Override with command line arguments if provided
    import argparse
    parser = argparse.ArgumentParser(description='Primer interaction analysis tool')
    parser.add_argument('--dir', type=str, help='Output directory path')
    parser.add_argument('--temp', type=int, help='Temperature for analysis (°C)')
    parser.add_argument('--dg_min', type=float, help='Minimum ΔG for homodimers/hairpins')
    parser.add_argument('--dg_max', type=float, help='Maximum ΔG for homodimers/hairpins')
    parser.add_argument('--dg_fr_min', type=float, help='Minimum ΔG for forward-reverse interactions')
    parser.add_argument('--dg_fr_max', type=float, help='Maximum ΔG for forward-reverse interactions')
    parser.add_argument('--dg_po_min', type=float, help='Minimum ΔG for primer-probe interactions')
    parser.add_argument('--dg_po_max', type=float, help='Maximum ΔG for primer-probe interactions')

    args = parser.parse_args()

    # Update params with provided arguments
    if args.dir:
        params['directory'] = args.dir
    if args.temp:
        params['temperature'] = args.temp
    if args.dg_min:
        params['delta_G_min'] = args.dg_min
    if args.dg_max:
        params['delta_G_max'] = args.dg_max
    if args.dg_fr_min:
        params['delta_G_FR_min'] = args.dg_fr_min
    if args.dg_fr_max:
        params['delta_G_FR_max'] = args.dg_fr_max
    if args.dg_po_min:
        params['delta_G_PO_min'] = args.dg_po_min
    if args.dg_po_max:
        params['delta_G_PO_max'] = args.dg_po_max

    # Run the analysis
    try:
        directory = params.pop('directory')
        interaction_check(directory, **params)
    except Exception as e:
        print(f"Error: {e}")
        import traceback
        traceback.print_exc()


if __name__ == "__main__":
    # Set start method to 'spawn' for better compatibility
    mp.set_start_method('spawn', force=True)
    main()