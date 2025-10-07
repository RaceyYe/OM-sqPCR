from typing import List, Tuple

from modules.auxiliary_modules.preliminary_exclusion import cDNASecondaryStructureDetector  # Adjust the import path accordingly


def detect_secondary_structures(
    forward_strand: str,
    reverse_strand: str,
    min_len: int = 4,
    max_len: int = 15,
    score_threshold: float = -5.0
) -> Tuple[List[List[int]], List[List[int]]]:
    """
    Detect potential secondary structure regions in cDNA sequences.

    Args:
        forward_strand (str): 5'→3' sequence of the forward strand.
        reverse_strand (str): 5'→3' sequence of the reverse strand (should be reverse complement of forward).
        min_len (int): Minimum complementary region length to consider.
        max_len (int): Maximum complementary region length to consider.
        score_threshold (float): Threshold for revised ΔG to be considered significant.

    Returns:
        Tuple[
            List[List[int]],  # intra_strand matches: [pos1_start, pos1_end, pos2_start, pos2_end]
            List[List[int]]   # inter_strand matches: [pos1_start, pos1_end, pos2_start, pos2_end]
        ]
    """
    detector = cDNASecondaryStructureDetector(forward_strand, reverse_strand, min_length=min_len, max_length=max_len)
    result_df = detector.process_sequences(forward_strand, reverse_strand)

    if result_df.empty:
        return [], []

    # Filter by score threshold
    result_df = result_df[result_df['revised_dg'] <= score_threshold]

    intra_matches: List[List[int]] = []
    inter_matches: List[List[int]] = []

    for _, row in result_df.iterrows():
        match = [row['pos1_start'], row['pos1_end'], row['pos2_start'], row['pos2_end']]
        if row['type'] == 'intra_strand':
            intra_matches.append(match)
        elif row['type'] == 'inter_strand':
            inter_matches.append(match)

    return intra_matches, inter_matches
