from typing import List, Tuple
from intervaltree import IntervalTree, Interval
from modules.auxiliary_modules.PE_API import detect_secondary_structures


class PrimerRangeChecker:
    """
    Checks whether candidate primer ranges are fully contained in any exclusion regions,
    using IntervalTree for accurate and efficient interval search.
    """

    def __init__(self, intra_matches: List[List[int]], inter_matches: List[List[int]]):
        """
        Initialize with intra and inter exclusion regions.

        Args:
            intra_matches: List of [start, end, ..., ...] for intra-matches.
            inter_matches: Same as above for inter-matches.
        """
        all_regions = intra_matches + inter_matches
        self.tree = IntervalTree(
            Interval(r[0], r[1]+1 , True) for r in all_regions
        )

    def is_valid(self, candidate: Tuple[int, int]) -> bool:
        s, e = candidate
        print(f"[DEBUG] Checking candidate range {s}-{e}")
        for interval in self.tree.overlap(s, e + 1):  # Consider adding +1 if off-by-one
            print(f"[DEBUG] Overlaps with interval: {interval.begin}-{interval.end}")
            if interval.begin <= s and interval.end >= e + 1:
                print("[INFO] This range failed (fully contained)")
                return False
        print("[INFO] This range succeeded")
        return True

    # def is_valid(self, candidate: Tuple[int, int]) -> bool:
    #     """
    #     Return True if candidate range is NOT fully contained within any exclusion region.
    #
    #     Args:
    #         candidate: Tuple (start, end) of the primer candidate.
    #
    #     Returns:
    #         bool: True if not fully contained in any exclusion region.
    #     """
    #     s, e = candidate
    #     for interval in self.tree.overlap(s, e):
    #         if interval.begin <= s and interval.end >= e:
    #             return False  # Fully contained → invalid
    #     return True

    def are_valid(self, candidates: List[Tuple[int, int]]) -> List[bool]:
        """
        Check a list of candidate primer ranges and return a list of booleans
        indicating whether each candidate is valid (not fully contained in any exclusion region).

        Args:
            candidates: List of candidate primer ranges as (start, end) tuples.

        Returns:
            List[bool]: List of validity statuses corresponding to the input candidates.
        """
        validity = []
        for candidate in candidates:
            s, e = candidate
            contained = False
            # Check if any overlapping interval fully contains candidate
            for interval in self.tree.overlap(s, e):
                if interval.begin <= s and interval.end >= e:
                    contained = True
                    break  # No need to check other overlaps
            validity.append(not contained)
        return validity


# Example usage with multiple candidates
if __name__ == '__main__':
    upper = "AACTATACAATCTACTACCTCAGTATGTGGGAGCGGTCGGGTCCAGATATTCGTATCTGTCGACCAGAGTGTGGGCTCCCACATAC"
    lower = "GTATGTGGGAGCCCACACTCTGGTCGACAGATACGAATATCTGGACCCGACCGCTCCCACATACTGAGGTAGTAGATTGTATAGTT"
    intra, inter = detect_secondary_structures(upper, lower)

    checker = PrimerRangeChecker(intra, inter)

    # Test with multiple candidates
    candidates = [(1, 30), (10, 40), (50, 70)]
    validity = checker.are_valid(candidates)

    print(f"Validity results: {validity}")
    # Output example: [False, True, True] (if (1,30) is invalid but others are valid)

# from typing import List, Tuple
# from intervaltree import IntervalTree, Interval
# from modules.auxiliary_modules.PE_API import detect_secondary_structures
# 
# class PrimerRangeChecker:
#     """
#     Checks whether candidate primer ranges are fully contained in any exclusion regions,
#     using IntervalTree for accurate and efficient interval search.
#     """
# 
#     def __init__(self, intra_matches: List[List[int]], inter_matches: List[List[int]]):
#         """
#         Initialize with intra and inter exclusion regions.
# 
#         Args:
#             intra_matches: List of [start, end, ..., ...] for intra-matches.
#             inter_matches: Same as above for inter-matches.
#         """
#         all_regions = intra_matches + inter_matches
#         self.tree = IntervalTree(
#             Interval(r[0], r[1], True) for r in all_regions
#         )
# 
#     def is_valid(self, candidate: Tuple[int, int]) -> bool:
#         """
#         Return True if candidate range is NOT fully contained within any exclusion region.
# 
#         Args:
#             candidate: Tuple (start, end) of the primer candidate.
# 
#         Returns:
#             bool: True if not fully contained in any exclusion region.
#         """
#         s, e = candidate
#         for interval in self.tree.overlap(s, e):
#             if interval.begin <= s and interval.end >= e:
#                 return False  # Fully contained → invalid
#         return True
# 
# 
# # Example usage
# if __name__ == '__main__':
#     upper = "AACTATACAATCTACTACCTCAGTATGTGGGAGCGGTCGGGTCCAGATATTCGTATCTGTCGACCAGAGTGTGGGCTCCCACATAC"
#     lower = "GTATGTGGGAGCCCACACTCTGGTCGACAGATACGAATATCTGGACCCGACCGCTCCCACATACTGAGGTAGTAGATTGTATAGTT"
#     intra, inter = detect_secondary_structures(upper, lower)
# 
#     checker = PrimerRangeChecker(intra, inter)
#     candidate = (1, 30)
# 
#     if checker.is_valid(candidate):
#         print("Candidate is valid")
#     else:
#         print("Candidate is fully contained in an exclusion region")
