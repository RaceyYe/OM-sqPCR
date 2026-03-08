"""
PrimerProbeDesigner.py
Module for designing primers and probes from cDNA sequences
"""

import os
import sys
import csv
import argparse
from typing import List, Tuple, Dict, Union
# # from modules.auxiliary_modules.helper_functions import helpers
# To this (relative import):
# from .auxiliary_modules.cDNA_generation import cDNA_generator
# from .criteria import primer_design
# # (sRNA) racey@ROG:~/Desktop$ python -m primerDesign.modules.FP_design -u TGAGGTAGTAGATTGTATAGTT -l GTATGTGGGAGCCCACACTCTGGTCGACAGATACGAATATCTGGACCCGACCGCTCCCACATAC -o "/home/racey/Desktop/"


# from modules.auxiliary_modules.cDNA_generation import cDNA_generator
from modules.criteria import primer_design

sys.path.append("/home/racey/Desktop/primerDesign/auxiliary_modules")
sys.path.append(os.path.abspath(os.path.join(os.path.dirname(__file__), "..")))

from modules.auxiliary_modules.forwardPrimer_extension import (
    optimize_tm_with_iteration,
    # optimize_tm_cdna_with_iteration,
    save_primers_to_csv,
)
from modules.auxiliary_modules.SecStructure import analyze_primers_selves
# from modules.auxiliary_modules.PE_API import detect_secondary_structures
# from modules.auxiliary_modules.PE_checker import PrimerRangeChecker


class FP_designer(primer_design):
    """Main class for primer and probe design"""

    def __init__(self, min_5ext=None, max_5ext=6, **kwargs):
        super().__init__(**kwargs)
        self.FP_base = self.cDNA.get("upper", "Unknown")[:12]
        self.max_5ext = max_5ext

        # Store temperature for later use
        self.temp = kwargs.get('temp', 37)

        # ###  Validity check
        # intra, inter = detect_secondary_structures(self.cDNA.get(("upper", "")),self.cDNA.get(("lower","")))
        # checker = PrimerRangeChecker(intra, inter)
        # ###

        # Handle min_5ext logic
        if min_5ext is not None:
            self.min_5ext = min_5ext
        elif hasattr(self, 'min_5ext') and self.min_5ext == self.max_5ext:
            pass
        else:
            self.min_5ext = 0

    def design_primers(self, output_dir, temp) -> List[Tuple]:
        """Main primer design workflow using optimize_tm_with_iteration"""

        # Extract the preserved sRNA sequence (convert to RNA format for the function)
        srna_seq = self.cDNA.get("upper", "Unknown")[:len(self.preserved_sRNA)]
        # Convert DNA back to RNA format for the optimize_tm_with_iteration function
        srna_rna = srna_seq.replace('T', 'U')

        # Get the stem-loop sequence (lower strand)
        stemloop_seq = self.cDNA.get("lower", "Unknown")

        print(f"Designing primers for sRNA sequence: {srna_rna}")
        print(f"Using stem-loop sequence: {stemloop_seq}")
        print(f"Base sequence for forward primers: {self.FP_base}")
        print(f"Temperature: {temp}°C")
        print(f"Tm range: {self.tm_min} - {self.tm_max}")
        print(f"ΔG range: {self.dg_min} - {self.dg_max}")

        # Use the optimize_tm_with_iteration function with proper parameters
        try:
            valid_primers = optimize_tm_with_iteration(
                srna_seq=srna_rna,
                # stemloop_seq=stemloop_seq,
                Tm_min=self.tm_min,
                Tm_max=self.tm_max,
                delta_G_min=self.dg_min,
                delta_G_max=self.dg_max,
                temperature=temp
            )

            print(f"Found {len(valid_primers)} valid primers within the specified criteria")

            # Save primers to CSV using the function from forwardPrimer_extension
            if valid_primers:
                save_primers_to_csv(valid_primers, output_dir, temp)
                print(f"Primers saved to {output_dir}")
            else:
                print("No valid primers met the specified criteria.")

            return valid_primers

        except Exception as e:
            print(f"Error in primer design: {str(e)}")
            raise


    def _filter_by_secondary(self, primers: List[Tuple]) -> List[Tuple]:
        """Apply secondary structure filtering including range validation"""
        sequences = [p[0] for p in primers]
        print(f"Quantity of filtered sequences based on Tm selection: {len(sequences)}")

        # Check secondary structures of primers themselves
        dg_df = analyze_primers_selves(sequences, experimental_temperature=self.temp)
        dg_df.columns = [f"{col}_FP" for col in dg_df.columns]

        # Check primer ranges against exclusion regions
        range_validity = self.check_primer_ranges(primers)

        filtered = []
        for i, (primer, tm, ext_len) in enumerate(primers):
            if not range_validity[i]:
                continue  # Skip primers in exclusion regions

            matches = dg_df[dg_df['Sequence_FP'] == primer]
            if not matches.empty:
                row = matches.iloc[0]
                if (self.dg_min <= row[f'dG_HD_{self.temp}C_FP'] <= self.dg_max and
                        self.dg_min <= row[f'dG_HP_{self.temp}C_FP'] <= self.dg_max):
                    filtered.append((
                        primer,
                        tm,
                        row[f'dG_HD_{self.temp}C_FP'],
                        row[f'dG_HP_{self.temp}C_FP'],
                        ext_len
                    ))
        return filtered

def main():
    parser = argparse.ArgumentParser(description='Design PCR primers and probes')
    parser.add_argument('-u', '--upper_strands',type=str, required=True,
                        help='upper one of the cDNA strands (both 5\'->3\')') # dest='u',
    parser.add_argument('-l', '--lower_strands',type=str, required=True,
                        help='lower one of the cDNA strands (both 5\'->3\')') # dest='l',
    parser.add_argument('-o', '--output', required=True, help='Output directory')
    parser.add_argument('-b', '--tm_base', type=float, default=59, help='Base TM for variance')
    parser.add_argument('-v', '--tm_var',type=float, default=2, help='TM variance')
    parser.add_argument('--dg_min', type=float, default=-4.5, help='Minimum ΔG')
    parser.add_argument('--dg_max', type=float, default=4.5, help='Maximum ΔG')
    parser.add_argument('--temp', type=float, default=37, help='Reaction temperature')
    parser.add_argument('--min_5ext', type=int, default=0, help='Minimum 5\' extension length')
    parser.add_argument('--max_5ext', type=int, default=6, help='Maximum 5\' extension length')

    args = parser.parse_args()

    # Validate and create output directory
    if not os.path.exists(args.output):
        os.makedirs(args.output, exist_ok=True)

    # Instantiate FP_designer with extracted u and l, and other parameters from args
    designer = FP_designer(
        u=args.upper_strands,
        l=args.lower_strands,
        o=args.output,
        b=args.tm_base,
        v=args.tm_var,
        dg_min=args.dg_min,
        dg_max=args.dg_max,
        temp=args.temp,
        min_5ext=args.min_5ext,
        max_5ext=args.max_5ext
    )

    print(f"[INFO] upper strand: {args.upper_strands}")
    print(f"[INFO] lower strand: {args.lower_strands}")
    # Proceed with primer design
    try:
        valid_primers = designer.design_primers(output_dir=args.output, temp=args.temp)
        print(f"Successfully designed {len(valid_primers)} primers")
    except ValueError as e:
        print(f"Design failed: {str(e)}")
        sys.exit(1)
    except Exception as e:
        print(f"Unexpected error: {str(e)}")
        sys.exit(1)


if __name__ == '__main__':
    main()