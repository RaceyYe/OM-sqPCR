import primer3
import pandas as pd

def analyze_primers_selves(sequences, experimental_temperature=37):
    # Ensure experimental_temperature is an integer to avoid decimal in column names
    experimental_temperature = int(experimental_temperature)

    # Define default columns to include even if no sequences
    default_columns = [
        'Sequence',
        'Structure_Found_HD', 'Melting_Temp_Tm_HD', f'dG_HD_{experimental_temperature}C',
        'dH_HD', 'dS_HD', 'Structure_Found_HP', 'Melting_Temp_Tm_HP',
        f'dG_HP_{experimental_temperature}C', 'dH_HP', 'dS_HP'
    ]

    if not sequences:
        # Return an empty DataFrame with the expected columns
        return pd.DataFrame(columns=default_columns)

    """
    Analyze DNA sequences for secondary structures (homodimers & hairpins).

    Parameters:
        sequences (list): A list of DNA sequences.

    Returns:
        pd.DataFrame: A DataFrame containing thermodynamic properties.
    """
    # Create an empty list to store results
    results = []

    # Loop through each sequence and calculate secondary structure properties
    for seq in sequences:
        seq = str(seq)  # <--- Add this line
        # Calculate homodimer properties
        hd_result = primer3.calc_homodimer(seq, output_structure=True)

        # Calculate hairpin properties
        hp_result = primer3.calc_hairpin(seq, output_structure=True)

        # Store results in a dictionary with suffixes (_HD for homodimer, _HP for hairpin)
        result_dict = {
            'Sequence': seq,

            # Homodimer properties
            'Structure_Found_HD': hd_result.structure_found,
            'Melting_Temp_Tm_HD': hd_result.tm,
            f'dG_HD_{experimental_temperature}C': hd_result.dg/1000,  # kcal/mol, ΔG at 37°C
            'dH_HD': hd_result.dh,
            'dS_HD': hd_result.ds,

            # Hairpin properties
            'Structure_Found_HP': hp_result.structure_found,
            'Melting_Temp_Tm_HP': hp_result.tm,
            f'dG_HP_{experimental_temperature}C': hp_result.dg/1000,  # kcal/mol, ΔG at 37°C
            'dH_HP': hp_result.dh,
            'dS_HP': hp_result.ds
        }

        # Append the dictionary to results list
        results.append(result_dict)

    # Convert results list into a pandas DataFrame
    df = pd.DataFrame(results)

    # Convert ΔG from J/mol to kcal/mol at specified temperature
    # df[f'dG_HD_{experimental_temperature}C_check'] = (df['dH_HD'] - ((273.15 + experimental_temperature) * df["dS_HD"])) / 1000
    # df[f'dG_HP_{experimental_temperature}C_check'] = (df['dH_HP'] - ((273.15 + experimental_temperature) * df["dS_HP"])) / 1000

    return df  # Return the DataFrame


def analyze_primers_interaction(seq_1, seq_2, experimental_temperature=37):
    """
    Analyze heterodimer formation between two DNA sequences.

    Parameters:
        seq_1 (str): First DNA sequence.
        seq_2 (str): Second DNA sequence.

    Returns:
        pd.DataFrame: A DataFrame containing heterodimer thermodynamic properties.
    """
    # Calculate heterodimer properties
    result = primer3.calc_heterodimer(seq_1, seq_2, output_structure=True)

    # Store results in a dictionary
    result_dict = {
        'Sequence_1': seq_1,
        'Sequence_2': seq_2,
        'Structure_Found_HET': result.structure_found,
        'Melting_Temp_Tm_HET': result.tm,
        f'dG_HET_{experimental_temperature}C': result.dg/1000,  # kcal/mol, ΔG at 37°C
        'dH_HET': result.dh,
        'dS_HET': result.ds
    }

    # Convert ΔG from J/mol to kcal/mol at specified temperature
    # result_dict[f'dG_HET_{experimental_temperature}C'] = (result_dict['dH_HET'] - (
            # (273.15 + experimental_temperature) * result_dict["dS_HET"])) / 1000

    # Convert dictionary into a DataFrame
    df = pd.DataFrame([result_dict])

    return df  # Return the DataFrame


# Example usage:
if __name__ == "__main__":
    # Define a list of sequences
    sequences = ["CGCCGTTGAGGTAGTAGATTGT","GGGAGCCCACACTCTGGTC","ACGAATATCTGGACCCGACCGCTC","TGCCGCTGAGGTAGTAGATTGTAT"]
    # Run the analysis
    df = analyze_primers_selves(sequences)
    # Display the results
    print(df.columns)
    print(df[['dG_HD_37C','dG_HP_37C']])
    # df.to_csv("/home/racey/Desktop/Projects/sRNA/Outputs/test/SecStrucutre_1.csv")

    # Define two sequences for interaction testing
    # seq1 = "GGTGTCGTGGAGTCAGTGCA"
    # seq2 = "TGCACT"
    # seq1="TTTTTAAAAATTTTTTAAAAATTTTAAAAATTTTTTTTTT"
    # seq2="TTTTTTTT"
    seq1="CCCCAAA"
    seq2="GGGG"
    # seq1="TTTTTTTTTTCCCCTTTTTTTTTT"
    # seq2="TTTTTTTTGGGGTTTTTTTTTTTT"
    # Analyze heterodimer formation
    print("[PART 2]")
    df_2 = analyze_primers_interaction(seq1, seq2)
    # Display results
    print(df_2.columns)
    print(df_2['dG_HET_37C'])
    # df_2.to_csv("/home/racey/Desktop/Projects/sRNA/Outputs/test/SecStrucutre_2.csv")