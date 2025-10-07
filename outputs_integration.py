import os
import re
import argparse
import pandas as pd
from scipy.optimize import linear_sum_assignment


def load_head_dfs(base_dir):
    """
    Scan subdirectories of base_dir for 'Integrated_df.csv' files,
    load each into a DataFrame, take head(10), and store in a dict keyed by subdir name.
    Skips directories without valid CSV files.
    """
    dfs = {}
    missing_combinations = []

    # Check if base_dir exists
    if not os.path.exists(base_dir):
        print(f"Error: Input directory '{base_dir}' does not exist.")
        return dfs, missing_combinations

    total_dirs = 0
    valid_dirs = 0

    for name in os.listdir(base_dir):
        subdir_path = os.path.join(base_dir, name)
        if not os.path.isdir(subdir_path):
            continue

        total_dirs += 1
        subdir = os.path.join(subdir_path, name)
        csv_path = os.path.join(subdir, "Integrated_df.csv")

        if not os.path.exists(subdir):
            missing_combinations.append(name)
            print(f"Warning: Expected subdirectory '{subdir}' not found for combination '{name}'")
            continue

        if not os.path.isfile(csv_path):
            missing_combinations.append(name)
            print(f"Warning: No Integrated_df.csv found for combination '{name}'")
            continue

        try:
            df = pd.read_csv(csv_path)
            if df.empty:
                missing_combinations.append(name)
                print(f"Warning: Empty DataFrame for combination '{name}'")
                continue

            if 'Score' not in df.columns:
                missing_combinations.append(name)
                print(f"Warning: DataFrame in '{name}' has no 'Score' column.")
                continue

            dfs[name] = df.head(3).copy()
            valid_dirs += 1

        except Exception as e:
            missing_combinations.append(name)
            print(f"Error loading '{csv_path}': {str(e)}")

    print(f"Processed {total_dirs} directories: {valid_dirs} valid, {len(missing_combinations)} missing/invalid")

    return dfs, missing_combinations


def summarize_scores(dfs):
    """
    Compute average 'Score' for each combination DataFrame.
    """
    if not dfs:
        print("No valid combinations to summarize")
        return pd.DataFrame(columns=['Combination', 'Average_Score', 'Rank'])

    avg_scores = {name: df['Score'].mean() for name, df in dfs.items()}
    summary_df = pd.DataFrame(
        list(avg_scores.items()),
        columns=['Combination', 'Average_Score']
    )
    summary_df['Rank'] = summary_df['Average_Score'].rank(method='min', ascending=True).astype(int)
    summary_df = summary_df.sort_values('Rank')
    return summary_df


def filter_mutual_exclusive(summary_df):
    """
    Keep only combos where elements are not the same.
    Elements parsed from name like 'Task_#_<A>_x_<B>'.
    """
    # Store the original summary for returning as full_summary
    full_summary = summary_df.copy()

    pattern = re.compile(r'^Task_\d+_.+?_x_.+$')
    valid = [c for c in summary_df['Combination'] if pattern.match(c) and not re.search(r'_(\w+)_x_\1', c)]

    filtered_df = summary_df[summary_df['Combination'].isin(valid)].reset_index(drop=True)

    if filtered_df.empty and not full_summary.empty:
        print("Warning: No combinations match the expected pattern after filtering")

    return filtered_df, full_summary


def select_best_matching(summary_df):
    """
    Select the best set of combos (one per unique A and B) maximizing total score.
    Returns list of best combinations and total score.
    """
    if summary_df.empty:
        return [], 0  # Handle empty summary case

    pattern = re.compile(r'^Task_\d+_(.+?)_x_(.+)$')
    entries, As, Bs = [], set(), set()

    for _, row in summary_df.iterrows():
        combo, score = row['Combination'], row['Average_Score']
        m = pattern.match(combo)
        if not m:
            continue
        a, b = m.group(1), m.group(2)
        As.add(a)
        Bs.add(b)
        entries.append((a, b, combo, score))

    if not entries:
        return [], 0  # Handle case with no valid entries

    As, Bs = sorted(As), sorted(Bs)
    Ai = {a: i for i, a in enumerate(As)}
    Bi = {b: j for j, b in enumerate(Bs)}

    # Initialize cost matrix
    W = [[0] * len(Bs) for _ in range(len(As))]
    combo_map = {}
    combo_map = {}

    for a, b, combo, score in entries:
        i, j = Ai[a], Bi[b]
        W[i][j] = score
        combo_map[(i, j)] = combo

    # Handle case where some (i,j) positions have no corresponding combo
    cost = [[-W[i][j] for j in range(len(Bs))] for i in range(len(As))]  # Use negative for maximization

    row_ind, col_ind = linear_sum_assignment(cost)
    best_combos = [combo_map.get((i, j)) for i, j in zip(row_ind, col_ind) if (i, j) in combo_map]
    total_score = sum(W[i][j] for i, j in zip(row_ind, col_ind) if (i, j) in combo_map)

    return best_combos, total_score


def get_unique_sheet_name(writer, sheet_name, max_length=31):
    """
    Ensures unique sheet names by adding a suffix if needed.
    """
    if not sheet_name:
        sheet_name = "Sheet"

    # Remove invalid characters for Excel sheet names
    invalid_chars = r'[\\/*?:[\]]'
    sheet_name = re.sub(invalid_chars, '_', sheet_name)

    original = sheet_name[:max_length]
    unique_name = original
    suffix = 1

    while unique_name in writer.sheets:
        suffix_str = f"_{suffix}"
        unique_name = original[:(max_length - len(suffix_str))] + suffix_str
        suffix += 1

    return unique_name


def integrate_outputs(input_dir: str, output_dir: str):
    """
    Runs the full summary+selection pipeline on `input_dir`,
    writing CSV/Excel into `output_dir`.
    Returns (summary_path, excel_path).
    """
    # Create output directory if it doesn't exist
    os.makedirs(output_dir, exist_ok=True)

    # Load data
    dfs, missing_combinations = load_head_dfs(input_dir)
    if not dfs:
        print("No valid data found. Exiting.")
        return None, None

    # Summarize
    summary_df = summarize_scores(dfs)
    filtered_df, full_summary_df = filter_mutual_exclusive(summary_df)

    # Determine output paths
    out_dir = output_dir if os.path.isdir(output_dir) else os.path.dirname(output_dir) or os.getcwd()
    csv_path = os.path.join(out_dir, 'summary_scores.csv')
    excel_path = os.path.join(out_dir, 'Selected_combinations.xlsx')

    # Save summary
    full_summary_df.to_csv(csv_path, index=False)
    print(f"Summary written to {csv_path}")

    # Save missing combinations to file
    if missing_combinations:
        missing_path = os.path.join(out_dir, 'missing_combinations.txt')
        with open(missing_path, 'w') as f:
            f.write('\n'.join(missing_combinations))
        print(f"List of {len(missing_combinations)} missing combinations written to {missing_path}")

    # Select best combos
    best_combos, total_score = select_best_matching(filtered_df)
    if best_combos:
        print(f"Best {len(best_combos)} combinations:")
        for combo in best_combos:
            print(f"  {combo}")
        print(f"Total {'min' if all(filtered_df['Average_Score'] >= 0) else 'max'} score: {total_score:.4f}")
    else:
        print("No valid combinations found for selection.")

    # Write Excel with: Summary, Missing Combinations, Filtered Summary, Priority sheets, and best combos
    try:
        with pd.ExcelWriter(excel_path, engine='xlsxwriter') as writer:
            # Write full summary
            full_summary_df.to_excel(writer, sheet_name='Summary', index=False)

            # Add missing combinations sheet
            if missing_combinations:
                pd.DataFrame({'Missing_Combinations': missing_combinations}).to_excel(
                    writer, sheet_name='Missing_Combinations', index=False)

            # Mark the filtered combinations in a separate sheet
            if not filtered_df.empty:
                filtered_df.to_excel(writer, sheet_name='Filtered_Summary', index=False)

            # Add Priority sheets (top 10)
            for name, df in dfs.items():
                if name.startswith('Priority'):
                    df = df.copy()
                    df['Rank'] = df['Score'].rank(method='min', ascending=True).astype(int)
                    sheet = get_unique_sheet_name(writer, name)
                    df.to_excel(writer, sheet_name=sheet, index=False)

            # Add best combo sheets, ensuring no duplicates
            for combo in best_combos:
                if combo in dfs:
                    df = dfs[combo].copy()
                    df['Rank'] = df['Score'].rank(method='min', ascending=True).astype(int)
                    sheet = get_unique_sheet_name(writer, combo)
                    df.to_excel(writer, sheet_name=sheet, index=False)
                else:
                    print(f"Warning: Selected combo '{combo}' not found in loaded DataFrames")

        print(f"Selected combinations and Priority data written to {excel_path}")
    except Exception as e:
        print(f"Error writing Excel file: {str(e)}")
        excel_path = None

    return csv_path, excel_path


def main():
    parser = argparse.ArgumentParser(
        description="Select best mutually-exclusive task combos covering all elements"
    )
    parser.add_argument('-d', '--directory', required=True,
                        help="Parent directory with subdirs containing Integrated_df.csv")
    parser.add_argument('-o', '--output', default='output_folder',
                        help="Directory or file prefix for summary CSV and Excel output")
    args = parser.parse_args()

    try:
        integrate_outputs(args.directory, args.output)
    except Exception as e:
        print(f"Error during integration: {str(e)}")
        import traceback
        traceback.print_exc()


if __name__ == '__main__':
    main()