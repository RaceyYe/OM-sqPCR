from modules.auxiliary_modules.helper_functions import helpers
import pandas as pd
import os
from modules.auxiliary_modules.cDNA_generation import cDNA_generator
from modules.criteria import primer_design
from modules.auxiliary_modules.SecStructure import analyze_primers_selves, analyze_primers_interaction
from modules.auxiliary_modules.parameter_calculator import Primer5ThermodynamicCalculator
from typing import List, Tuple
from intervaltree import IntervalTree, Interval
# from modules.auxiliary_modules.PE_API import detect_secondary_structures
# from modules.auxiliary_modules.PE_checker import PrimerRangeChecker

helper = helpers()
calculator = Primer5ThermodynamicCalculator()


class RP_designer(primer_design):
    def __init__(self, selected_seq, GAP=25, PR_gap_max=6, min_RP_len=18, max_RP_len=26, min_probe_len=18, max_probe_len=26,
                 **kwargs):  # revise GAP from 22 into 25
        super().__init__(**kwargs)
        # self.min_gap = GAP + self.get_gap(self.output_dir, "min")
        # self.max_gap = GAP + self.get_gap(self.output_dir, "max")
        # self.FR_gap = GAP  # redundant design # the distance from forward primer's 3' end to the reverse primer's 3' end should larger than this value for designing the probe carrier
        self.selected_seq = selected_seq
        self.PR_gap = PR_gap_max
        self.min_RP_len = min_RP_len
        self.max_RP_len = max_RP_len
        self.min_probe_len = min_probe_len
        self.max_probe_len = max_probe_len
        # Build exclusion region checker for negative strand
        self.lower_tree = self._build_lower_strand_checker()

    # def _build_lower_strand_checker(self) -> PrimerRangeChecker:
    #     """Create exclusion checker with coordinates transformed to negative strand"""
    #     upper = self.cDNA.get("upper", "")
    #     lower = self.cDNA.get("lower", "")
    #     total_len = len(upper)
    #
    #     # Detect exclusion regions in positive strand coordinates
    #     intra, inter = detect_secondary_structures(upper, lower)
    #
    #     # Transform regions to negative strand coordinates
    #     converted_regions = []
    #     for region in intra + inter:
    #         # Transform [start, end] to [total_len - end - 1, total_len - start - 1]
    #         new_start = total_len - region[1] - 1
    #         new_end = total_len - region[0] - 1
    #         converted_regions.append([new_start, new_end])
    #
    #     return PrimerRangeChecker(converted_regions, [])

    def get_gap(self, output_dir, argument: str):
        FP_path = os.path.join(output_dir, "designed_FP.csv")
        try:
            min_gap, max_gap, position_range = helper.generate_gap_range(FP_path)
        except FileNotFoundError:
            raise FileNotFoundError(f"FP design file not found at {FP_path}")

        if argument == "min":
            return min_gap
        elif argument == "max":
            return max_gap
        else:
            return None

    def sliding(self):
        '''Generate all valid RP-probe pairs meeting constraints'''
        # print(self.cDNA)
        reverse_strand = self.cDNA.get("lower", "")
        if not reverse_strand:
            print("failed to fetch the reverse strand")
            return pd.DataFrame()

        fp_min_gap = self.get_gap(self.output_dir, "min")
        fp_max_gap = self.get_gap(self.output_dir, "max")
        print(f"fp_min_gap: {fp_min_gap}; fp_max_gap: {fp_max_gap}")
        results = []

        # Get total strand length for easier reference
        strand_length = len(reverse_strand)
        print(f"length of reverse strand: {strand_length}")
        print(f"the reverse strand sequence: 5'- {reverse_strand} -3'.")
        print(f"Forbidden region: [{fp_min_gap},{fp_max_gap}]")

        # Count total candidates considered and accepted
        total_considered = 0
        total_accepted = 0

        # Track counts for each forbidden region
        forbidden_region_counts = {}

        # Iterate through all possible forbidden regions (FP integration areas)
        for forbidden in range(fp_max_gap , fp_min_gap -1  , -1):
            print(f"Current forbidden region length is {forbidden}")
            forbidden_region_considered = 0
            forbidden_region_accepted = 0

            # Instead of calculating available length once, we'll check constraints directly
            # for each potential RP position

            # Maximum possible start position for RP
            max_rp_start = strand_length - forbidden - self.min_probe_len - self.min_RP_len
            print(f"available length: {max_rp_start}")

            # Iterate through possible RP starting positions
            for rp_start in range(0, max_rp_start + 1):  # +1 to include the max position
                if rp_start % 5 == 0:  # Reduce output verbosity
                    print(f"    Current rp_start is at {rp_start}")

                # Generate RP sequences of different lengths
                for rp_len in range(self.min_RP_len, self.max_RP_len + 1):
                    rp_end = rp_start + rp_len

                    # Check if RP would overflow into space needed for probe
                    if rp_end > strand_length - forbidden - self.min_probe_len:
                        continue

                    rp_seq = reverse_strand[rp_start:rp_end]
                    total_considered += 1
                    forbidden_region_considered += 1

                    # Define reverse primer range (inclusive)
                    rp_range = (rp_start, rp_start + rp_len - 1)



                    valid = self.lower_tree.is_valid(rp_range)
                    print(f"[DEBUG] rp_range: {rp_range}, is_valid: {valid}")
                    if not valid:
                        print({"[INFO] This range failed"})
                        continue
                    print({"[INFO] This range succeeded"})

                    # # Check reverse primer against exclusion regions
                    # if not self.lower_tree.is_valid(rp_range):
                    #     print("[INFO] This range failed")
                    #     continue
                    # print("[INFO] This range succeeded")


                    # Generate probe sequences with gap constraint
                    for gap in range(0, self.PR_gap + 1):
                        probe_start = rp_end + gap

                        # Check if there's still room for probe
                        if probe_start > strand_length - forbidden - self.min_probe_len:
                            continue

                        # Calculate maximum allowed probe length based on remaining space
                        max_allowed_probe_len = min(self.max_probe_len, strand_length - probe_start - forbidden)

                        # Skip if we can't fit a minimum sized probe
                        if max_allowed_probe_len < self.min_probe_len:
                            continue

                        # Iterate through valid probe lengths
                        for probe_len in range(self.min_probe_len, max_allowed_probe_len + 1):
                            # print(f"current probe_len: {probe_len}")
                            probe_end = probe_start + probe_len - 1

                            # Define probe range (inclusive)
                            probe_range = (probe_start, probe_end)

                            # Check probe against exclusion regions
                            if not self.lower_tree.is_valid(probe_range):
                                continue

                            # Final boundary check
                            if probe_end <= strand_length - forbidden:
                                probe_seq = reverse_strand[probe_start:probe_end]
                                results.append({
                                    'RP_Sequence': rp_seq,
                                    'Probe_Sequence': probe_seq,
                                    'RP_Length': rp_len,
                                    'Probe_Length': probe_len,
                                    'RP_Start': rp_start,
                                    'Probe_Start': probe_start,
                                    'Gap': gap,
                                    'Forbidden_Region': forbidden
                                })
                                total_accepted += 1
                                forbidden_region_accepted += 1

            # Store and display counts for this forbidden region
            forbidden_region_counts[forbidden] = {
                'considered': forbidden_region_considered,
                'accepted': forbidden_region_accepted
            }
            print(
                f"Forbidden region {forbidden}: Considered {forbidden_region_considered}, Accepted {forbidden_region_accepted}")

        print(f"Total candidates considered: {total_considered}")
        print(f"Total candidates accepted: {total_accepted}")

        # Print summary of all forbidden regions
        print("\nSummary of candidates by forbidden region:")
        for forbidden, counts in forbidden_region_counts.items():
            print(f"Forbidden region {forbidden}: Considered {counts['considered']}, Accepted {counts['accepted']}")

        return pd.DataFrame(results)

    def design_probes(self):
        '''Main design workflow'''
        df = self.sliding()
        if df.empty:
            print("there are nothing within the dataframe")
            return df

        # # Save results to file
        # output_path = os.path.join(self.output_dir, "RP_probe_pairs.csv")
        # df.to_csv(output_path, index=False)
        # print(f"Results saved to {output_path}")

        # Return filtered DataFrame
        return df  # [['RP_Sequence', 'Probe_Sequence', 'Gap', 'Forbidden_Region']]

    def secondary_structure_check(self, df):
        # Step 0: Pre-check and sequence string conversion
        if df.empty:
            print("[INFO] Input DataFrame is empty.")
            return pd.DataFrame()

        # Adjust ΔG values for RP_Start == 0
        mask = df['RP_Start'] == 0
        df.loc[mask, f'dG_HP_{self.temp}C_RP'] -= 1
        df.loc[mask, f'dG_HET_{self.temp}C_HET'] -= 1

        # Ensure sequences are strings
        df['RP_Sequence'] = df['RP_Sequence'].apply(lambda x: str(x) if hasattr(x, 'translate') else x)
        df['Probe_Sequence'] = df['Probe_Sequence'].apply(lambda x: str(x) if hasattr(x, 'translate') else x)

        # Step 1: Remove duplicates while preserving the one with the largest Forbidden_Region
        before_dedup = len(df)
        df = df.sort_values('Forbidden_Region', ascending=False)
        df = (
            df
            .drop_duplicates(subset=['RP_Sequence', 'Probe_Sequence', 'RP_Start', 'Probe_Start', 'Gap'], keep='first')
            .copy()
        )
        after_dedup = len(df)
        print(f"[INFO] Removed {before_dedup - after_dedup} duplicates by keeping the largest 'Forbidden_Region'.")

        # Step 2: Extract unique sequences
        rp_sequences = df['RP_Sequence'].unique().tolist()
        probe_sequences = df['Probe_Sequence'].unique().tolist()
        if not rp_sequences or not probe_sequences:
            print("[INFO] No sequences found for analysis.")
            return pd.DataFrame()

        # Step 3: Analyze sequences
        try:
            # Self-structure analysis
            rp_results = analyze_primers_selves([str(seq) for seq in rp_sequences])
            rp_results.columns = [f"{col}_RP" for col in rp_results.columns]
            rp_results['RP_Sequence'] = rp_sequences

            probe_results = analyze_primers_selves([str(seq) for seq in probe_sequences])
            probe_results.columns = [f"{col}_Probe" for col in probe_results.columns]
            probe_results['Probe_Sequence'] = probe_sequences

            # Cross-dimer analysis
            cross_dimer = []
            for _, row in df.iterrows():
                res = analyze_primers_interaction(str(row['RP_Sequence']), str(row['Probe_Sequence']))
                cross_dimer.append(res.iloc[0].to_dict())
            cross_dimer_df = pd.DataFrame(cross_dimer)
            cross_dimer_df.columns = [f"{col}_HET" for col in cross_dimer_df.columns]

        except Exception as e:
            print(f"[ERROR] Analysis failed: {e}")
            return pd.DataFrame()

        # Merge all results
        merged_df = (
            df
            .merge(rp_results, on='RP_Sequence', how='left')
            .merge(probe_results, on='Probe_Sequence', how='left')
        )
        merged_df = pd.concat([merged_df.reset_index(drop=True), cross_dimer_df.reset_index(drop=True)], axis=1)

        # Calculate melting temperatures (Tm)
        merged_df['Tm_RP'] = merged_df['RP_Sequence'].apply(lambda s: calculator.calculate_tm(seq=str(s).upper()))
        merged_df['Tm_Probe'] = merged_df['Probe_Sequence'].apply(lambda s: calculator.calculate_tm(seq=str(s).upper()))
        merged_df['Tm_RP'] = pd.to_numeric(merged_df['Tm_RP'], errors='coerce').fillna(0)
        merged_df['Tm_Probe'] = pd.to_numeric(merged_df['Tm_Probe'], errors='coerce').fillna(0)

        # Prepare ΔG columns
        dg_cols = [
            f'dG_HD_{self.temp}C_RP', f'dG_HP_{self.temp}C_RP',
            f'dG_HD_{self.temp}C_Probe', f'dG_HP_{self.temp}C_Probe'
        ]
        het_col = f'dG_HET_{self.temp}C_HET'
        merged_df[dg_cols] = merged_df[dg_cols].fillna(-999)
        merged_df[het_col] = merged_df[het_col].fillna(-999)

        # Step 4: Filtering with logging
        # ΔG filter
        mask_dg = (
                          (merged_df[dg_cols] >= self.dg_min) & (merged_df[dg_cols] <= self.dg_max)
                  ).all(axis=1) & (
                          (merged_df[het_col] >= self.dg_min_PP) & (merged_df[het_col] <= self.dg_max_PP)
                  )
        df_dg = merged_df[mask_dg]
        print(f"[INFO] Removed {len(merged_df) - len(df_dg)} rows failing ΔG filter.")

        # Tm_RP filter
        mask_tm_rp = (df_dg['Tm_RP'] >= self.tm_min) & (df_dg['Tm_RP'] <= self.tm_max)
        df_rp = df_dg[mask_tm_rp]
        print(f"[INFO] Removed {len(df_dg) - len(df_rp)} rows failing Tm_RP filter.")

        # Determine which temperature range to use based on selected_seq comparison
        print(f"[REPORT]: received selected_seq: {self.selected_seq}.")
        if self.selected_seq is not None and self.preserved_stemloop == self.selected_seq:
            print("[INFO] Using second temperature range (tm_min_probe_2/tm_max_probe_2) for probe filtering")
            mask_tm_probe = (df_rp['Tm_Probe'] >= self.tm_min_probe_2) & (df_rp['Tm_Probe'] <= self.tm_max_probe_2)
        else:
            print("[INFO] Using first temperature range (tm_min_probe_1/tm_max_probe_1) for probe filtering")
            mask_tm_probe = (df_rp['Tm_Probe'] >= self.tm_min_probe_1) & (df_rp['Tm_Probe'] <= self.tm_max_probe_1)

        df_probe = df_rp[mask_tm_probe]
        print(f"[INFO] Removed {len(df_rp) - len(df_probe)} rows failing Tm_Probe filter.")

        # Step 5: Final filtering and export
        before_g = len(df_probe)
        result_df = df_probe[~df_probe['Probe_Sequence'].str.startswith('G')]
        print(f"[INFO] Removed {before_g - len(result_df)} sequences starting with 'G'.")

        output_file = os.path.join(self.output_dir, "designed_RP_PC.csv")
        result_df.to_csv(output_file, index=False)
        print(f"[INFO] Saved {len(result_df)} final rows to {output_file}.")

        return result_df


def main():
    # Example initialization with selected_seq parameter
    selected_seq = "GTATGTGGGAGCCCACACTCTGGTCGACAGATACGAATATCTGGACCCGACCGCTCCCACATAC"

    designer = RP_designer(
        u="UGAGGUAGUAGAUUGUAUAGUU",
        l="GTATGTGGGAGCCCACACTCTGGTCGACAGATACGAATATCTGGACCCGACCGCTCCCACATAC",
        o="/home/racey/Desktop/primerDesign/Outputs/let7f",
        PR_gap_max=6
    )

    # Generate all possible pairs
    result_df = designer.design_probes()
    # Save raw results for debugging
    result_df.to_csv(os.path.join(designer.output_dir, "raw_RP_probe_pairs.csv"), index=False)
    print(f"Raw candidate count: {len(result_df)}")

    # Pass selected_seq to the secondary_structure_check method
    final_df = designer.secondary_structure_check(result_df, selected_seq=selected_seq)

    print(final_df['Forbidden_Region'].value_counts())


if __name__ == "__main__":
    main()