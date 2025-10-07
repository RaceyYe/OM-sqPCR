import pandas as pd
import argparse
import os
import re
from modules.auxiliary_modules.SecStructure import analyze_primers_selves


class RatingCheck:
    __Types = {"FP", "RP"}

    def __init__(self, input_data, seq_types=None, temp=37, threshold=90):
        self.input_data = input_data
        self.seq_types = seq_types
        self.temp = temp
        self.threshold = threshold

    def check_sequence_type(self, colnames):
        detected = []
        for t in self.__Types:
            pattern = re.compile(rf'.*_{t}$')
            for colname in colnames:
                if pattern.search(colname):
                    detected.append(t)
                    break
        print(f"Detected types are as follows: {detected}")
        return detected

    def _filter_by_rating(self, df, seq_type):
        """Filter the dataframe by rating score based on dG values."""
        dg_hd_col = f'dG_HD_{self.temp}C_{seq_type}'
        dg_hp_col = f'dG_HP_{self.temp}C_{seq_type}'

        if dg_hd_col not in df.columns or dg_hp_col not in df.columns:
            raise ValueError(f"Missing columns {dg_hd_col} or {dg_hp_col} for rating filter")

        df[f'Rating_{seq_type}'] = 100 + 1.8 * df[dg_hd_col] + 1.4 * df[dg_hp_col]
        filtered_df = df[df[f'Rating_{seq_type}'] > self.threshold].copy()

        print(f"{len(filtered_df)} sequences passed the rating filter for {seq_type}")
        return filtered_df

    def run_from_dataframe(self, df):
        cols = df.columns.tolist()
        inferred_types = self.check_sequence_type(cols) if not self.seq_types else self.seq_types

        for seq_type in inferred_types:
            matched_cols = [col for col in cols if col.endswith(f"_{seq_type}")]
            if not matched_cols:
                print(f"No column found ending with '_{seq_type}'")
                continue
            colname = matched_cols[0]

            sequences = df[colname].dropna().astype(str).tolist()
            print(f"Analyzing {len(sequences)} sequences for {seq_type} from column '{colname}'...")

            dg_df = analyze_primers_selves(sequences, experimental_temperature=self.temp)

            for col in dg_df.columns:
                if col != "Sequence":
                    suffix_col = f"{col}_{seq_type}"
                    df[suffix_col] = dg_df.set_index("Sequence").reindex(df[colname].astype(str)).reset_index(drop=True)[col]

            df = self._filter_by_rating(df, seq_type)

        return df

    def run_from_list(self, sequences, seq_type):
        print(f"Running sequence list for {seq_type}")
        df = analyze_primers_selves(sequences, experimental_temperature=self.temp)

        # Compute rating
        df[f'Rating_{seq_type}'] = 100 + 1.8 * df[f'dG_HD_{self.temp}C_{seq_type}'] + 1.4 * df[f'dG_HP_{self.temp}C_{seq_type}']
        df = df[df[f'Rating_{seq_type}'] > self.threshold].copy()

        print(df.head())
        return df

class CombinationCheck:
    '''
    This check evaluates combinations based on:
    1. Longer forward primers are better (FP_WOE + 6)
    2. Smaller reverse-probe gaps are better (Gap)
    3. Closer RP starts to the 5' end are better (RP_Start)
    '''

    def __init__(self, output_dir, weights=(1.0, -2.0, -2.0)):
        self.output_dir = output_dir
        self.weights = weights
        self.df = pd.read_csv(os.path.join(output_dir, "combined_df.csv"))

    def _systematic_score(self):
        w1, w2, w3 = self.weights

        # Sanity checks for required columns
        for col in ['FP_WOE', 'Gap', 'RP_Start']:
            if col not in self.df.columns:
                raise ValueError(f"Missing required column: {col}")

        # Calculate components
        self.df['Forward_Extension'] = self.df['FP_WOE'] + 6
        self.df['RP_Probe_Gap'] = self.df['Gap']
        # self.df['RP_Start'] = self.df['RP_Start']

        # Normalize values (to 0–1 range) to ensure fair weighted sum
        norm = lambda x: (x - x.min()) / (x.max() - x.min() + 1e-6)

        self.df['norm_ext'] = norm(self.df['Forward_Extension'])      # higher better
        self.df['norm_gap'] = 1 - norm(self.df['RP_Probe_Gap'])       # lower better
        self.df['norm_rp_start'] = 1 - norm(self.df['RP_Start'])      # lower better

        # Compute final score
        self.df['systematic_score'] = (
            w1 * self.df['norm_ext'] +
            w2 * self.df['norm_gap'] +
            w3 * self.df['norm_rp_start']
        )

        # Optional: sort by score descending
        self.df = self.df.sort_values(by='systematic_score', ascending=False)

        return self.df

    def save_scored_dataframe(self, filename="scored_combined_df.csv"):
        scored_df = self._systematic_score()
        output_path = os.path.join(self.output_dir, filename)
        scored_df.to_csv(output_path, index=False)
        print(f"Saved scored dataframe to {output_path}")
        return scored_df




def main():
    parser = argparse.ArgumentParser(description="Check thermodynamic properties of primers")
    parser.add_argument('--input', type=str, required=True, help="CSV file path or comma-separated sequence list")
    parser.add_argument('--type', nargs='+', help="Type(s) to check: FP, RP, Probe")
    parser.add_argument('--temp', type=int, default=37, help="Temperature for delta G calculation")
    parser.add_argument('--output', type=str, help="Output CSV path (if input is a file)")

    args = parser.parse_args()

    if os.path.isfile(args.input):
        print(f"Reading from file: {args.input}")
        df = pd.read_csv(args.input)
        checker = RatingCheck(df, args.type, temp=args.temp)
        updated_df = checker.run_from_dataframe(df)

        if args.output:
            updated_df.to_csv(args.output, index=False)
            print(f"Output saved to {args.output}")
        else:
            print(updated_df.head())

    else:
        try:
            sequences = [s.strip() for s in args.input.split(',')]
            if not args.type or len(args.type) != 1:
                raise ValueError("For sequence list input, please specify exactly one --type")
            checker = RatingCheck(None, args.type, temp=args.temp)
            df = checker.run_from_list(sequences, args.type[0])
            print(f"the file was checked by the ")
            print(df)
        except Exception as e:
            print(f"Error parsing input: {e}")


if __name__ == "__main__":
    main()
