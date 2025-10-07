import pandas as pd
from SecStructure import analyze_primers_interaction
import os

def load_primer_files(file_dict):
    """
    Load primer data from CSV files given in a dictionary.

    Args:
        file_dict (dict): Dictionary with keys "FP" and "RP_PC" where
            file_dict["FP"] is the path to the forward primers CSV file
            and file_dict["RP_PC"] is the path to the reverse primer/probe pairs CSV file.

    Returns:
        tuple: (forward_df, reverse_df) DataFrames.
            forward_df is expected to have columns like "Sequence", "Tm", and "ΔG_pm (25 ℃)".
            reverse_df is expected to have columns like "Reverse_Primer", "Tm_RP", "ΔG_RP",
            "Probe_Carrier", "Tm_PC", "ΔG_PC", and "ΔG_HET".
    """
    if not isinstance(file_dict, dict):
        raise TypeError("Expected a dictionary for file_dict, but got: " + str(type(file_dict)))
    forward_df = pd.read_csv(file_dict["FP"], encoding='gbk')
    reverse_df = pd.read_csv(file_dict["RP_PC"], encoding='gbk')
    return forward_df, reverse_df


def integrate_primer_data(forward_df, reverse_df, delta_G_min=-9, delta_G_max=9, temperature=25):
    """
    Integrate forward primers with reverse primer/probe pairs.

    For each reverse primer candidate, this function selects forward primers whose Tm (Tm_FP)
    is within ±1°C of the reverse primer’s Tm (Tm_RP). For each pairing, it calculates the
    heterodimer ΔG values between forward and reverse as well as forward and probe using the
    analyze_primers_interaction() function. Only combinations with all cross-interaction ΔG values
    within the target range [delta_G_min, delta_G_max] are preserved.

    The integrated result contains:
      - Forward_Primer, Reverse_Primer, Probe_Carrier
      - Their Tm values (suffixes _FP, _RP, _PC)
      - The absolute difference in Tm between forward and reverse (Tm_diff)
      - The self-dimer ΔG values for the three primers (ΔG_FP, ΔG_RP, ΔG_PC)
      - The heterodimer ΔG values among each pair (ΔG_FP_RP, ΔG_FP_PC, ΔG_RP_PC)

    Args:
        forward_df (pd.DataFrame): DataFrame with forward primer data.
        reverse_df (pd.DataFrame): DataFrame with reverse primer/probe pair data.
        delta_G_min (float): Minimum acceptable heterodimer ΔG value.
        delta_G_max (float): Maximum acceptable heterodimer ΔG value.
        temperature (float): Temperature for thermodynamic calculations (default 25°C).

    Returns:
        pd.DataFrame: The integrated primer pairs DataFrame.
    """
    integrated_results = []
    # print(reverse_df.columns)
    # Iterate over each reverse primer / probe pair candidate.
    for _, rev_row in reverse_df.iterrows():
        rp = rev_row["Reverse_Primer"]
        tm_rp = rev_row["Tm_RP"]
        dg_rp = rev_row[f'ΔG_RP_HD_{temperature}C']  # self-dimer ΔG for reverse primer
        pc = rev_row["Probe_Carrier"]
        tm_pc = rev_row["Tm_PC"]
        dg_pc = rev_row[f'ΔG_PC_HD_{temperature}C']  # self-dimer ΔG for probe carrier
        dg_rp_pc = rev_row[f"ΔG_HET_{temperature}C"]  # heterodimer ΔG for reverse-probe pair

        # Select forward primers with Tm_FP within ±1°C of the reverse primer’s Tm.
        candidates = forward_df[
            (forward_df["Tm_FP"] >= tm_rp - 1) & (forward_df["Tm_FP"] <= tm_rp + 1)
            ]

        for _, fwd_row in candidates.iterrows():
            fp = fwd_row["Forward_Primer"]
            tm_fp = fwd_row["Tm_FP"]
            dg_fp = fwd_row[f'ΔG_FP_HD_{temperature}C']  # self-dimer ΔG for forward primer

            # Calculate heterodimer ΔG for forward–reverse.
            df_fp_rp = analyze_primers_interaction(fp, rp)
            dg_fp_rp = df_fp_rp[f'Delta_G_{temperature}C_HET (kcal/mol)'].iloc[0]

            # Calculate heterodimer ΔG for forward–probe.
            df_fp_pc = analyze_primers_interaction(fp, pc)
            dg_fp_pc = df_fp_pc[f'Delta_G_{temperature}C_HET (kcal/mol)'].iloc[0]

            # Check that all heterodimer ΔG values are within the target range.
            if (delta_G_min <= dg_fp_rp <= delta_G_max and
                    delta_G_min <= dg_fp_pc <= delta_G_max and
                    delta_G_min <= dg_rp_pc <= delta_G_max):
                integrated_results.append({
                    "Forward_Primer": fp,
                    "Reverse_Primer": rp,
                    "Probe_Carrier": pc,
                    "Tm_FP": tm_fp,
                    "Tm_RP": tm_rp,
                    "Tm_PC": tm_pc,
                    "Tm_diff": abs(tm_rp - tm_fp),
                    f'ΔG_FP_HD_{temperature}C': dg_fp,
                    f'ΔG_RP_HD_{temperature}C': dg_rp,
                    f'ΔG_PC_HD_{temperature}C': dg_pc,
                    f"ΔG_FP_RP_{temperature}C": dg_fp_rp,
                    f"ΔG_FP_PC_{temperature}C": dg_fp_pc,
                    f"ΔG_RP_PC_{temperature}C": dg_rp_pc
                })

    integrated_df = pd.DataFrame(integrated_results)
    return integrated_df

def save_integrated_data(integrated_df, output_dir):
    """
    Save the integrated primer data to a CSV file, sorted by proximity to target Tm (59°C).

    Args:
        integrated_df (pd.DataFrame): The integrated primer pairs DataFrame.
        output_dir (str): Directory path to save the output CSV file.
    """
    # Calculate mean Tm
    integrated_df["Mean_Tm"] = (integrated_df["Tm_FP"] + integrated_df["Tm_RP"]) / 2

    # Calculate absolute difference from target Tm (59°C) and sort
    integrated_df["Tm_Deviation"] = abs(integrated_df["Mean_Tm"] - 59)
    integrated_df = integrated_df.sort_values(by="Tm_Deviation")

    # Remove temporary deviation column (optional)
    # integrated_df = integrated_df.drop(columns=["Tm_Deviation"])

    # Save to CSV
    output_file = os.path.join(output_dir, "Integrated_df.csv")
    integrated_df.to_csv(output_file, encoding="gbk", index=False)
    print(f"Integrated primer pairs saved to {output_file}")

# def save_integrated_data(integrated_df, output_dir):
#     """
#     Save the integrated primer data to a CSV file.
#
#     Args:
#         integrated_df (pd.DataFrame): The integrated primer pairs DataFrame.
#         output_file (str): Path to the output CSV file.
#     """
#     integrated_df["Mean_Tm"] = (integrated_df["Tm_FP"]+integrated_df["Tm_RP"])/2
#
#     output_file = os.path.join(output_dir,"Integrated_df.csv")
#     integrated_df.to_csv(output_file,encoding ="gbk", index=False)
#     print(f"Integrated primer pairs saved to {output_file}")


def interaction_check(file_dict, delta_G_min=-9, delta_G_max=9, temperature=25):
    """
    Main function to load primer files from the given dictionary, integrate them,
    and save the output CSV.

    Args:
        file_dict (dict): Dictionary with keys "FP" and "RP_PC" pointing to the respective CSV files.
        output_file (str): Path for the integrated output CSV.
        delta_G_min (float): Minimum acceptable heterodimer ΔG value.
        delta_G_max (float): Maximum acceptable heterodimer ΔG value.
        temperature (float): Temperature for thermodynamic calculations (default 25°C).
    """
    forward_df, reverse_df = load_primer_files({
    "FP": os.path.join(file_dict,"designed_FP.csv"),  # forward primers CSV file
    "RP_PC": os.path.join(file_dict,"designed_RP_PC.csv")  # reverse primer/probe pairs CSV file
    })
    integrated_df = integrate_primer_data(forward_df, reverse_df, delta_G_min, delta_G_max, temperature)
    save_integrated_data(integrated_df, file_dict)



# If running as a script for testing:
if __name__ == "__main__":
    class Args:
        delta_G_min = -8  # Example threshold for heterodimer/hairpin filtering (kcal/mol)
        delta_G_max = 8
        temperature = 25
        output_dir = r"D:/Codes"

    args = Args()
    # Example file dictionary. Update these paths as needed.
    file_dict = {
        "FP": "FP.csv",  # forward primers CSV file
        "RP_PC": "RP_PC.csv"  # reverse primer/probe pairs CSV file
    }
    interaction_check(args.output_dir,args.delta_G_min, args.delta_G_max, args.temperature)