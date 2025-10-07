import pandas as pd
import primer3
from typing import List, Tuple, Dict, Optional

class cDNASecondaryStructureDetector:
    """
    A module to detect potential secondary structure forming regions in cDNA sequences
    for primer design optimization.
    """

    def __init__(self, forward_strand, reverse_strand, min_length: int = 4, max_length: int = 15):
        self.FS = forward_strand
        # now complement_SM returns the complemented string, and reverse_SM reverses it
        inferred_RS = self.reverse_SM(self.complement_SM(forward_strand))
        if reverse_strand != inferred_RS:
            raise ValueError(
                "reverse_strand is not the exact reverse complement of forward_strand;"
                " inputs are not cDNA in correct format"
            )
        self.min_length = min_length
        self.max_length = max_length

    @staticmethod
    def complement_SM(seq: str) -> str:
        comp_map = str.maketrans('ATGC', 'TACG')
        return seq.translate(comp_map)

    @staticmethod
    def reverse_SM(seq: str) -> str:
        return seq[::-1]

    def merge_intervals(intervals: List[Tuple[int, int]]) -> List[Tuple[int, int]]:
        if not intervals:
            return []
        sorted_intervals = sorted(set(intervals), key=lambda x: x[0])
        merged = []
        start_prev, end_prev = sorted_intervals[0]
        for start_curr, end_curr in sorted_intervals[1:]:
            if start_curr <= end_prev:
                end_prev = max(end_prev, end_curr)
            else:
                merged.append((start_prev, end_prev))
                start_prev, end_prev = start_curr, end_curr
        merged.append((start_prev, end_prev))
        return merged

    def find_complementary_regions(
            self,
            seq1: str,
            seq2: str,
            exclude_direct: bool = True
    ) -> List[Tuple[int, int, int, str]]:
        """
        Find reverse-complementary regions between two 5'→3' strands.

        Returns:
            List of (pos1_in_seq1, mapped_pos2_in_seq2, match_length, query_seq)
        """
        matches = []
        seq2_35 = self.reverse_SM(seq2)  # Reverse of 5'→3' = 3'→5'
        len1, len2 = len(seq1), len(seq2_35)

        for i in range(len1):
            for L in range(self.min_length, self.max_length + 1):
                if i + L > len1:
                    continue
                query = seq1[i:i + L]
                complementary_query = self.complement_SM(query)

                for j in range(len2 - L + 1):
                    candidate = seq2_35[j:j + L]
                    if candidate == complementary_query:
                        mapped_pos2 = len(seq2) - (j + L) # from 5' terminal #  j  # directly use the position from the 3' terminal
                        if i == -( mapped_pos2 + L - len(seq2)):
                            continue
                        matches.append((i, mapped_pos2, L, query))

        return matches


    def find_self_complementary_regions(
            self,
            target_seq: str,
            target_seq_pos: int,
            target_seq_length: int,
            seq: str
    ) -> List[Tuple[int, int, int, str]]:
        """
        Check whether a specific subsequence of `seq` has any exact
        reverse-complement elsewhere in `seq`.  Returns at most one match.

        Args:
            target_seq:       The subsequence you're testing (5'->3').
            target_seq_pos:   Its start index in `seq` (0-based).
            target_seq_length:Length of that subsequence.
            seq:              Full strand (5'->3').

        Returns:
            A singleton list [(target_pos, match_pos, length, target_seq)]
            or [] if no self‑complement is found.
        """
        # sanity check
        assert seq[target_seq_pos:target_seq_pos + target_seq_length] == target_seq, (
            "target_seq must match seq at the given position and length"
        )

        rc = self.reverse_SM(target_seq)
        seq_len = len(seq)
        end_pos = target_seq_pos + target_seq_length

        # slide window over entire seq, skip the original region
        for j in range(seq_len - target_seq_length + 1):
            if target_seq_pos <= j < end_pos:
                continue
            if self.complement_SM(seq[j:j + target_seq_length]) == rc:
                return [(target_seq_pos, j, target_seq_length, target_seq)]

        # no match found
        return []

    def extend_sequence(self, seq: str, start: int, end: int, extend_by: int = 1) -> Tuple[str, int, int]:
        """
        Extend a sequence region by specified bases on both 5' and 3' sides.

        Args:
            seq: Original sequence
            start: Start position
            end: End position (exclusive)
            extend_by: Number of bases to extend on each side

        Returns:
            Tuple of (extended_sequence, new_start, new_end)
        """
        new_start = max(0, start - extend_by)  # Extend on 5' side
        new_end = min(len(seq), end + extend_by)  # Extend on 3' side
        extended_seq = seq[new_start:new_end]

        return extended_seq, new_start, new_end

    def get_counterpart_sequence(self, seq: str, pos: int, length: int) -> Tuple[str, int, int, int]:
        """
        Extend by exactly 1 base on each side (unless hitting seq ends).
        """
        start_pos = max(0, pos - 1)
        end_pos = min(len(seq), pos + length + 1)
        extended_seq = seq[start_pos:end_pos]

        touched_5 = (start_pos == 0)
        touched_3 = (end_pos == len(seq))
        end_adjustments = int(touched_5) + int(touched_3)

        return extended_seq, start_pos, end_pos, end_adjustments

    def calculate_heterodimer_score(
            self,
            seq1: str,
            seq2: str
    ) -> float:
        """
        Calculate raw heterodimer ΔG (kcal/mol) via primer3.
        """
        try:
            result = primer3.calc_heterodimer(seq1, seq2)
            dg = getattr(result, 'dg', getattr(result, 'delta_g', 0.0)) / 1000.0
        except Exception:
            dg = 0.0
        return dg

    def process_sequences(self, upper_strand: str, lower_strand: str) -> pd.DataFrame:
        records = []

        # Inter-strand complementarity
        inter_matches = self.find_complementary_regions(upper_strand, lower_strand, exclude_direct=True)
        merged_inter = []
        if inter_matches:
            inter_sorted = sorted(inter_matches, key=lambda x: (x[0], x[1]))
            current = list(inter_sorted[0])
            for i in range(1, len(inter_sorted)):
                m = inter_sorted[i]
                if current[0] + current[2] == m[0] and current[1] + current[2] == m[1]:
                    current[2] += m[2]
                    current[3] = upper_strand[current[0]:current[0] + current[2]]
                else:
                    merged_inter.append(tuple(current))
                    current = list(m)
            merged_inter.append(tuple(current))

        for pos1, pos2, length, match_seq in merged_inter:
            target_window = upper_strand[pos1:pos1 + length]
            counterpart_seq, start2, end2, end_adjustments = self.get_counterpart_sequence(
                lower_strand, pos2, length
            )

            dg = self.calculate_heterodimer_score(target_window, counterpart_seq)
            revised = dg - end_adjustments

            records.append({
                'type': 'inter_strand',
                'pos1_start': pos1,
                'pos1_end': pos1 + length,
                'pos2_start': pos2,
                'pos2_end': pos2 + length,
                'counterpart_start': start2,
                'counterpart_end': end2,
                'match_length': length,
                'match_seq': match_seq,
                'counterpart_seq': counterpart_seq,
                'dg_kcal': dg,
                'end_adjustments': end_adjustments,
                'revised_dg': revised
            })

        # Intra-strand complementarity
        all_self = []
        for i in range(len(upper_strand)):
            for L in range(self.min_length, self.max_length + 1):
                if i + L > len(upper_strand):
                    continue
                tgt = upper_strand[i:i + L]
                hits = self.find_self_complementary_regions(tgt, i, L, upper_strand)
                all_self.extend(hits)

        # Merge intra-strand matches
        intra_dict = {}
        for pos1, pos2, length, match_seq in all_self:
            key = (min(pos1, pos2), max(pos1 + length, pos2 + length))
            if key not in intra_dict or length > intra_dict[key][2]:
                intra_dict[key] = (pos1, pos2, length, match_seq)

        for pos1, pos2, length, match_seq in intra_dict.values():
            target_window = upper_strand[pos1:pos1 + length]
            counterpart_seq, start2, end2, end_adjustments = self.get_counterpart_sequence(
                upper_strand, pos2, length
            )

            dg = self.calculate_heterodimer_score(target_window, counterpart_seq)
            revised = dg - end_adjustments

            records.append({
                'type': 'intra_strand',
                'pos1_start': pos1,
                'pos1_end': pos1 + length,
                'pos2_start': pos2,
                'pos2_end': pos2 + length,
                'counterpart_start': start2,
                'counterpart_end': end2,
                'match_length': length,
                'match_seq': match_seq,
                'counterpart_seq': counterpart_seq,
                'dg_kcal': dg,
                'end_adjustments': end_adjustments,
                'revised_dg': revised
            })

        df = pd.DataFrame.from_records(records)
        return df[df['revised_dg'] <= -5.0].reset_index(drop=True) if not df.empty else df

    def get_exclusion_regions(self, df: pd.DataFrame, score_threshold: float = -5.0) -> Dict[
        str, List[Tuple[int, int]]]:
        """
        Get regions to exclude from primer design based on secondary structure potential.

        Args:
            df: DataFrame from process_sequences
            score_threshold: Minimum score threshold for exclusion in kcal/mol (default: -5.0)

        Returns:
            Dictionary with strand names as keys and list of (start, end) tuples as values
        """
        exclusion_regions = {'upper': [], 'lower': []}

        # Filter based on score threshold
        significant_matches = df[df['revised_dg'] <= score_threshold]

        for _, row in significant_matches.iterrows():
            # Add regions from both positions of each match
            if row['type'] == 'inter_strand':
                # Upper strand region
                exclusion_regions['upper'].append((row['pos1_start'], row['pos1_end']))
                # Lower strand region
                exclusion_regions['lower'].append((row['pos2_start'], row['pos2_end']))
            elif row['type'] == 'intra_strand':
                # Both regions are on upper strand
                exclusion_regions['upper'].append((row['pos1_start'], row['pos1_end']))
                exclusion_regions['upper'].append((row['pos2_start'], row['pos2_end']))

        # Remove duplicates and sort
        for strand in exclusion_regions:
            exclusion_regions[strand] = sorted(list(set(exclusion_regions[strand])))

        return exclusion_regions


def main():
    # Example sequences (~30 nucleotides each)
    upper_strand = "AACTATACAATCTACTACCTCAGTATGTGGGAGCGGTCGGGTCCAGATATTCGTATCTGTCGACCAGAGTGTGGGCTCCCACATAC"
    lower_strand = "GTATGTGGGAGCCCACACTCTGGTCGACAGATACGAATATCTGGACCCGACCGCTCCCACATACTGAGGTAGTAGATTGTATAGTT"
    # upper_strand = "AAAAAAAAAAAAACCCC"
    # lower_strand = "GGGGTTTTTTTTTTTTT"
    # upper_strand = "AAAAAAAAAAAAAACCCCC"
    # lower_strand = "GGGGGTTTTTTTTTTTTTT"

    detector = cDNASecondaryStructureDetector(upper_strand, lower_strand, min_length=4, max_length=10)
    complementary_sequence = detector.reverse_SM(detector.complement_SM(upper_strand))
    print(f"Upper strand: 5' - {upper_strand} - 3'")
    print(f"Complementary: 5' - {complementary_sequence} - 3'")

    # Process sequences
    results_df = detector.process_sequences(upper_strand, lower_strand)

    # Save results if not empty
    if not results_df.empty:
        results_df.to_csv("/home/racey/Desktop/primerDesign_Outputs/test_outputs/preliminary_exclusion.csv",
                          index=False)

        # Display results
        print("\nDetected secondary structure regions:")
        print(results_df.to_string(index=False))

        # Get exclusion regions for primer design
        exclusion_regions = detector.get_exclusion_regions(results_df)

        print("\nRegions to exclude from primer design:")
        for strand, regions in exclusion_regions.items():
            print(f"{strand} strand: {regions}")
    else:
        print("\nNo significant secondary structure regions detected.")


if __name__ == "__main__":
    main()