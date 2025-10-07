#!/usr/bin/env python3
"""
Sequence Combination Generator

This module generates combinations of sRNA and stemloop sequences
and stores the results in a pandas DataFrame.
"""

import pandas as pd
import itertools
import argparse
from typing import Set, List, Tuple


def generate_combinations(srna_sequences: Set[str], stemloop_sequences: Set[str]) -> pd.DataFrame:
    """
    Generate all possible combinations of sRNA and stemloop sequences.

    Args:
        srna_sequences: A set of sRNA sequences
        stemloop_sequences: A set of stemloop sequences

    Returns:
        A pandas DataFrame with columns 'Task', 'sRNA_sequence', and 'Stemloop_sequence'
    """
    # Generate all combinations
    combinations = list(itertools.product(srna_sequences, stemloop_sequences))

    # Create a list of dictionaries for the DataFrame
    data = [
        {"Task": i + 1, "sRNA_sequence": srna, "Stemloop_sequence": stemloop}
        for i, (srna, stemloop) in enumerate(combinations)
    ]

    # Create DataFrame
    result_df = pd.DataFrame(data)

    return result_df


def save_dataframe(df: pd.DataFrame, output_file: str = "combinations.csv") -> None:
    """
    Save the DataFrame to a CSV file.

    Args:
        df: The DataFrame to save
        output_file: The output file path (default: combinations.csv)
    """
    df.to_csv(output_file, index=False)
    print(f"Results saved to {output_file}")


def main():
    """
    Main function for testing the module.
    """
    parser = argparse.ArgumentParser(description="Generate combinations of sRNA and stemloop sequences")
    parser.add_argument("--srna", nargs="+", help="List of sRNA sequences")
    parser.add_argument("--stemloop", nargs="+", help="List of stemloop sequences")
    parser.add_argument("--output", default="combinations.csv", help="Output CSV file")

    args = parser.parse_args()

    # If no arguments provided, use test data
    if args.srna is None or args.stemloop is None:
        print("No input sequences provided. Using test data...")

        # Example test data
        srna_sequences = {
            "ACGGUUACGGUUACGGA",
            "UUACGGACCGUUACGGU",
            "ACGGUUACGGACCGUUA"
        }

        stemloop_sequences = {
            "ACGGUUACGGUUACGGACCG",
            "UUACGGACCGUUACGGUACG",
            "ACGGUUACGGACCGUUACGG"
        }
    else:
        srna_sequences = set(args.srna)
        stemloop_sequences = set(args.stemloop)

    print(f"Processing {len(srna_sequences)} sRNA sequences and {len(stemloop_sequences)} stemloop sequences")

    # Generate combinations
    result_df = generate_combinations(srna_sequences, stemloop_sequences)

    # Display results
    print(f"Generated {len(result_df)} combinations")
    print(result_df.head(3))

    # Save results
    save_dataframe(result_df, args.output)

    return result_df


if __name__ == "__main__":
    main()