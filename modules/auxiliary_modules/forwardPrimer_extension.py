import itertools
import multiprocessing
import csv
import os
from functools import partial

from modules.auxiliary_modules.parameter_calculator import Primer5ThermodynamicCalculator
from modules.auxiliary_modules.inputCheck import AlignmentTransformer
from modules.auxiliary_modules.SecStructure import analyze_primers_selves


def tm_worker(calculator, target_range, variant_tuple):
    """
    Worker function for parallel Tm calculation.
    variant_tuple is (primer_sequence, FP_WOE) where FP_WOE is the length of the part without 5' extension in designed forward primer.
    """
    primer, FP_WOE = variant_tuple
    tm_value = calculator.calculate_tm(primer)
    return (primer, tm_value, FP_WOE) if target_range[0] <= tm_value <= target_range[1] else None


def compute_tm_parallel(seq_variants, target_range):
    """
    Computes Tm for all possible sequences in parallel.
    Filters sequences with Tm in the target range.
    Each variant is a tuple: (primer_sequence, gap_1)
    """
    calculator = Primer5ThermodynamicCalculator()
    worker_func = partial(tm_worker, calculator, target_range)

    with multiprocessing.Pool() as pool:
        results = pool.map(worker_func, seq_variants)

    # Filter out None results and return valid primers as (primer, tm_value, gap_1)
    return [result for result in results if result is not None]


def generate_primer_variants(srna_dna, min_5prime_ext=0, max_5prime_ext=6, max_3prime_ext=6):
    """
    Generate all possible primer variants based on sRNA sequence.
    Returns a list of tuples (primer_sequence, gap_1) where gap_1 is the length
    of the 5' extension.
    """
    bases = ['A', 'T', 'C', 'G']
    core_seq = srna_dna[:12]  # First 12 bases from 5' end
    variants = []

    # For each 3' extension length (including 0)
    for ext3_len in range(min(max_3prime_ext + 1, len(srna_dna) - 12)):
        fixed_ext = srna_dna[12:12 + ext3_len] if ext3_len > 0 else ""
        seq_with_fixed = core_seq + fixed_ext

        # Generate all possible 5' extensions and record the extension length as gap_1
        for ext5_len in range(min_5prime_ext, max_5prime_ext + 1):
            for flex_ext in itertools.product(bases, repeat=ext5_len):
                primer = ''.join(flex_ext) + seq_with_fixed
                variants.append((primer, ext5_len))
    return variants


def optimize_tm_with_iteration(srna_seq, Tm_min=57, Tm_max=61,
                               delta_G_min=-9, delta_G_max=9, temperature=25): #, stemloop_seq
    """
    Design forward primers with proper Tm range and favorable ΔG values.
    Propagates the 5' extension length as gap_1.
    """
    # transformer = AlignmentTransformer()
    # transformer.preliminary_check(srna_seq, stemloop_seq)

    # Convert RNA to DNA
    srna_dna = srna_seq.replace('U', 'T')
    target_tm_range = (Tm_min, Tm_max)

    # Generate and filter primer variants; each variant is (primer, gap_1)
    primer_variants = generate_primer_variants(srna_dna)
    valid_primers = compute_tm_parallel(primer_variants, target_tm_range)

    # Extract primer sequences for secondary structure analysis
    valid_sequences = [primer for primer, _, _ in valid_primers]

    # Analyze secondary structures
    dg_df = analyze_primers_selves(valid_sequences, experimental_temperature=temperature)
    dg_df.rename(columns={f'dG_HD_{temperature}C': f'dG_FP_HD_{temperature}C',
                          f'dG_HP_{temperature}C': f'dG_FP_HP_{temperature}C'}, inplace=True)

    # Build a dictionary mapping primer -> (tm_value, gap_1)
    tm_dict = {primer: (tm_value, gap_1) for primer, tm_value, gap_1 in valid_primers}

    # Filter based on ΔG values
    filtered_primers = []
    for _, row in dg_df.iterrows():
        seq = row['Sequence']
        if seq not in tm_dict:
            continue
        tm_value, gap_1 = tm_dict[seq]

        # Use standardized column names
        hd_col = f'dG_FP_HD_{temperature}C'
        hp_col = f'dG_FP_HP_{temperature}C'
        if hd_col not in row or hp_col not in row:
            raise ValueError(f"Expected columns {hd_col} and {hp_col} not found in DataFrame")

        worst_delta_g = min(row[hd_col], row[hp_col])
        if delta_G_min <= worst_delta_g <= delta_G_max:
            filtered_primers.append((
                seq,
                tm_value,
                row[hd_col],
                row[hp_col],
                gap_1
            ))
    return filtered_primers

# def optimize_tm_cdna_with_iteration(
#     cdna_seq,
#     Tm_min=57,
#     Tm_max=61,
#     delta_G_min=-9,
#     delta_G_max=9,
#     temperature=25
# ):
#     """
#     Design forward primers based on a cDNA sequence (DNA), analogous to optimize_tm_with_iteration but
#     without stem-loop alignment checks. The 5' end of cdna_seq should correspond to the sRNA-derived region
#     (T substituted for U).
#
#     Returns a list of tuples:
#         (primer_sequence, tm_value, dG_CD_HD_{temp}C, dG_CD_HP_{temp}C, gap_1)
#     """
#     # cdna_seq is already DNA (T instead of U)
#     target_tm_range = (Tm_min, Tm_max)
#
#     # Generate and screen primer variants
#     primer_variants = generate_primer_variants(cdna_seq)
#     valid_primers = compute_tm_parallel(primer_variants, target_tm_range)
#
#     # Extract sequences for secondary-structure analysis
#     valid_sequences = [primer for primer, _, _ in valid_primers]
#
#     # Analyze secondary structures
#     dg_df = analyze_primers_selves(valid_sequences, experimental_temperature=temperature)
#     # Rename to distinguish cDNA-based ΔG columns
#     dg_df.rename(
#         columns={
#             f'dG_HD_{temperature}C': f'dG_CD_HD_{temperature}C',
#             f'dG_HP_{temperature}C': f'dG_CD_HP_{temperature}C'
#         },
#         inplace=True
#     )
#
#     # Map primer -> (tm_value, gap_1)
#     tm_dict = {primer: (tm_value, gap_1) for primer, tm_value, gap_1 in valid_primers}
#
#     # Filter by ΔG
#     filtered_primers = []
#     hd_col = f'dG_CD_HD_{temperature}C'
#     hp_col = f'dG_CD_HP_{temperature}C'
#     for _, row in dg_df.iterrows():
#         seq = row['Sequence']
#         if seq not in tm_dict:
#             continue
#         tm_value, gap_1 = tm_dict[seq]
#         worst_delta_g = min(row[hd_col], row[hp_col])
#         if delta_G_min <= worst_delta_g <= delta_G_max:
#             filtered_primers.append((
#                 seq,
#                 tm_value,
#                 row[hd_col],
#                 row[hp_col],
#                 gap_1
#             ))
#     return filtered_primers


def save_primers_to_csv(primers, output_dir, temperature=25):
    """
    Save primer sequences with their characteristics to a CSV file.
    Includes the new column 'gap_1' for the 5' extension length.
    """
    transformer = AlignmentTransformer()
    calculator = Primer5ThermodynamicCalculator()
    output_csv_path = os.path.join(output_dir, "designed_FP.csv")

    with open(output_csv_path, 'w', newline='') as csvfile:
        fieldnames = ['Forward_Primer', 'Tm_FP', f'dG_HD_{temperature}C_FP','Length',
                      f'dG_HP_{temperature}C_FP', 'FP_WOE'] #f'dH ({temperature} C)', f'dS ({temperature} C)',
                      # f'dG_pm ({temperature} C)', 'A_Count', 'T_Count',
                      # 'G_Count', 'C_Count', 'GC%',
        writer = csv.DictWriter(csvfile, fieldnames=fieldnames)
        writer.writeheader()

        for primer_seq, tm_value, delta_g_hd, delta_g_hp, gap_1 in primers:
            # Calculate sequence characteristics
            characteristics = transformer.get_characteristics([primer_seq])
            seq_chars = characteristics[primer_seq]

            # Calculate thermodynamic parameters
            delta_h, delta_s = calculator.get_thermodynamic_params(primer_seq)
            delta_g = calculator.delta_G_calculator(delta_h, delta_s, temperature)

            writer.writerow({
                'Forward_Primer': primer_seq,
                'Tm_FP': round(tm_value, 2),
                # f'dH ({temperature} C)': round(delta_h, 2),
                # f'dS ({temperature} C)': round(delta_s, 2),
                # f'dG_pm ({temperature} C)': round(delta_g, 2),
                'Length': seq_chars['Length'],
                # 'A_Count': seq_chars['Counts'].get('A', 0),
                # 'T_Count': seq_chars['Counts'].get('T', 0),
                # 'G_Count': seq_chars['Counts'].get('G', 0),
                # 'C_Count': seq_chars['Counts'].get('C', 0),
                # 'GC%': round(seq_chars['GC%'], 2),
                f'dG_HD_{temperature}C_FP': round(delta_g_hd, 2),
                f'dG_HP_{temperature}C_FP': round(delta_g_hp, 2),
                'FP_WOE': seq_chars['Length'] - gap_1 # gap_1 is the 5 extension part in the forward strand
            })

    print(f"Primer data saved to {output_csv_path}")
    return output_csv_path


def main():
    """Test function for the primer selection process"""
    srna_seq = "UGAGGUAGUAGAUUGUAUAGUU"  # Example RNA sequence (5' to 3')
    # stemloop_seq = "GTATGTGGGAGCCCACACTCTGGTCGACAGATACGAATATCTGGACCCGACCGCTCCCACATAC"  # Example stem-loop sequence

    temperature = 37
    valid_primers = optimize_tm_with_iteration(
        srna_seq,
        # stemloop_seq,
        Tm_min=57,
        Tm_max=61,
        delta_G_min=-4.5,
        delta_G_max=10,
        temperature=temperature
    )

    print(f"Found {len(valid_primers)} valid primers within the target Tm range.")

    output_file = save_primers_to_csv(
        valid_primers,
        output_dir="/home/racey/Desktop/primerDesign_Outputs",
        temperature=temperature
    )

    print(f"Results saved to {output_file}")


if __name__ == "__main__":
    main()
