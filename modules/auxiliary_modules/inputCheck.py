import re

class AlignmentTransformer:
    """Class for DNA/RNA sequence manipulation, alignment, and analysis."""

    def __init__(self, alignment_type="RD"):
        """
        Initialize the transformer with the desired alignment type.
        Supported types: DD (DNA-DNA), DR (DNA-RNA), RD (RNA-DNA), RR (RNA-RNA)
        """
        if alignment_type not in {"DD", "DR", "RD", "RR"}:
            raise ValueError("Invalid alignment type. Use DD, DR, RD, or RR.")
        self.alignment_type = alignment_type

    @staticmethod
    def get_complement(sequence, seq_type="D2D"):
        """Get the complementary sequence while maintaining 5' to 3' annotation."""
        complement_rules = {
            'D2D': {'A': 'T', 'T': 'A', 'C': 'G', 'G': 'C'},
            'R2R': {'A': 'U', 'U': 'A', 'C': 'G', 'G': 'C'},
            'D2R': {'A': 'U', 'T': 'A', 'C': 'G', 'G': 'C'},
            'R2D': {'A': 'T', 'U': 'A', 'C': 'G', 'G': 'C'}
        }

        if seq_type not in complement_rules:
            raise ValueError(f"Invalid seq_type: {seq_type}. Choose from D2D, R2R, D2R, R2D.")

        comp_dict = complement_rules[seq_type]
        sequence = sequence.upper()

        try:
            complement = ''.join(comp_dict[base] for base in sequence)  # Keep 5' to 3' order
        except KeyError as e:
            raise ValueError(f"Invalid base {e} in sequence for {seq_type} mapping.") from e

        return complement  # Return in correct 5' to 3' orientation

    @staticmethod
    def get_characteristics(sequences):
        """Output and return basic characteristics of the input sequences."""
        characteristics = {}
        for seq in sequences:
            length = len(seq)
            counts = {base: seq.upper().count(base) for base in "ATCGU"}
            gc_content = ((counts.get('G', 0) + counts.get('C', 0)) / length) * 100 if length > 0 else 0
            characteristics[seq] = {
                'Length': length,
                'Counts': counts,
                'GC%': gc_content
            }
            # print(f"Sequence: 5'-{seq}-3'")
            # print(f"  Length: {length}, Counts: {counts}, GC%: {gc_content:.2f}%")
        return characteristics

    def determine_type(self, sequence):
        """Determine if the sequence is RNA or DNA."""
        return 'R' if 'U' in sequence else 'D'

    def preliminary_check(self, seq1, seq2):
        """Check if sequences match the expected DNA/RNA types."""
        seq1, seq2 = seq1.upper(), seq2.upper()


def main():
    transformer = AlignmentTransformer("RD")

    # Example sequences
    rna_seq = "UGUCAGUUUGUCAAAUACCCGA"
    dna_seq = "GTTGGCTCTGGTGCAGGGTCCGAGGTATTCGCACCAGAGCCAACTCGGGT"

    # Complementary Sequence
    print("\nComplementary Sequence:")
    transformer.get_complement(dna_seq, "D2D")

    # Sequence Characteristics
    print("\nSequence Characteristics:")
    transformer.get_characteristics([rna_seq, dna_seq])

    # Validate Alignment
    print("\nValidate Alignment:")
    transformer.validate_alignment(rna_seq, dna_seq, mode="33", align_len=6)


if __name__ == "__main__":
    main()
