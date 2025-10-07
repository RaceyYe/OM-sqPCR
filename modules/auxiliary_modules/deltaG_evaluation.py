import math


class WatsonCrickNN:
    def __init__(self):
        # Table 1 data from the paper
        self.nn_params = {
            ('AA', 'TT'): {'ΔH°': -7.6, 'ΔS°': -21.3, 'ΔG°37': -1.00},
            ('AT', 'TA'): {'ΔH°': -7.2, 'ΔS°': -20.4, 'ΔG°37': -0.88},
            ('TA', 'AT'): {'ΔH°': -7.2, 'ΔS°': -21.3, 'ΔG°37': -0.58},
            ('CA', 'GT'): {'ΔH°': -8.5, 'ΔS°': -22.7, 'ΔG°37': -1.45},
            ('GT', 'CA'): {'ΔH°': -8.4, 'ΔS°': -22.4, 'ΔG°37': -1.44},
            ('CT', 'AG'): {'ΔH°': -7.8, 'ΔS°': -21.0, 'ΔG°37': -1.28},
            ('GA', 'CT'): {'ΔH°': -8.2, 'ΔS°': -22.2, 'ΔG°37': -1.30},
            ('CG', 'GC'): {'ΔH°': -10.6, 'ΔS°': -27.2, 'ΔG°37': -2.17},
            ('GC', 'CG'): {'ΔH°': -9.8, 'ΔS°': -24.4, 'ΔG°37': -2.24},
            ('GG', 'CC'): {'ΔH°': -8.0, 'ΔS°': -19.9, 'ΔG°37': -1.84},
            'Initiation': {'ΔH°': +0.2, 'ΔS°': -5.7, 'ΔG°37': +1.96},
            'Terminal AT penalty': {'ΔH°': +2.2, 'ΔS°': +6.9, 'ΔG°37': +0.05},
            'Symmetry correction': {'ΔH°': 0.0, 'ΔS°': -1.4, 'ΔG°37': +0.43}
        }

    def get_params(self, base1, base2):
        key = f"{base1}{base2}/{base2}{base1}"
        if key in self.nn_params:
            return self.nn_params[key]
        else:
            raise KeyError(f"Base pair {base1}{base2} not found in table 1")


class InternalMismatch:
    def __init__(self):
        # Table 2 data from the paper would be initialized here
        # For simplicity, I'm creating a placeholder structure
        self.mismatch_params = {
            # This would need to be implemented based on the specific rules from Table 2 in the paper
            # The paper's Table 2 has 44 unique parameters for internal single mismatches
            # For the sake of this example, I'm creating a simplified version
            # the raw format is as follows, "CX/GY", in which the former two represent the paired bases, and in the dictionary, XY are sequential
            'GC': {'AA': 0.17, 'AC': 0.81, 'AG': -0.25, 'CA': 0.47,'CC': 0.79, 'CT': 0.62, 'GA': -0.52, 'GG': -1.11, "GT":0.08, 'TC': 0.98, 'TG': -0.59, "TT":0.45},
            'CG': {'AA': 0.43, 'AC': 0.75, 'AG': 0.03, 'CA': 0.79,'CC': 0.70, 'CT': 0.62, 'GA': -0.11, 'GG': -0.11, "GT":-0.47, 'TC': 0.40, 'TG': -0.32, "TT":-0.12},
            'AT': {'AA': 0.61, 'AC': 0.88, 'AG': 0.14, 'CA': 0.77,'CC': 1.33, 'CT': 0.64, 'GA': 0.02, 'GG': -0.13, "GT":0.71, 'TC': 0.73, 'TG': 0.07, "TT":0.69},
            'TA': {'AA': 0.69, 'AC': 0.92, 'AG': 0.42, 'CA': 1.33,'CC': 1.05, 'CT': 0.97, 'GA': 0.74, 'GG': 0.44, "GT":0.43, 'TC': 0.75, 'TG': 0.34, "TT":0.68}
        }

    def calculate(self, upper_strand, lower_strand):
        # Alignment algorithm would be implemented here
        # For demonstration, I'm assuming the strands are already aligned
        # and finding the first mismatch position

        # Reverse complement the lower strand for comparison
        complement = {'A': 'T', 'T': 'A', 'C': 'G', 'G': 'C'}
        aligned_lower = ''.join([complement[base] for base in reversed(lower_strand)])

        # Find mismatch positions
        mismatches = []
        for i in range(min(len(upper_strand), len(aligned_lower))):
            if upper_strand[i] != aligned_lower[i]:
                mismatches.append((i, upper_strand[i], aligned_lower[i]))

        # Calculate values for each mismatch
        results = []
        for pos, upper_base, lower_base in mismatches:
            # Determine neighboring bases
            left_upper = upper_strand[pos - 1] if pos > 0 else None
            right_upper = upper_strand[pos + 1] if pos < len(upper_strand) - 1 else None
            left_lower = aligned_lower[pos - 1] if pos > 0 else None
            right_lower = aligned_lower[pos + 1] if pos < len(aligned_lower) - 1 else None

            # Determine the context for lookup in Table 2
            # This would need to be implemented based on the specific rules from the paper
            context_key = f"{upper_base}{lower_base}"

            # Retrieve parameters from Table 2 (placeholder implementation)
            if context_key in self.mismatch_params:
                value = self.mismatch_params[context_key]
                results.append({
                    'position': pos,
                    'upper_base': upper_base,
                    'lower_base': lower_base,
                    'delta_G': value
                })
            else:
                raise KeyError(f"Mismatch {upper_base}{lower_base} not found in table 2")

        return results


class DanglingEnds:
    def __init__(self):
        # Table 3 data from the paper
        self.penalty_values = {
            '5_prime': {
                'AA/T': {'ΔH°': 0.2, 'ΔG°37': -0.51},
                'AC/G': {'ΔH°': -6.3, 'ΔG°37': -0.96},
                'AG/C': {'ΔH°': -3.7, 'ΔG°37': -0.58},
                'AT/A': {'ΔH°': -2.9, 'ΔG°37': -0.50},
                'CA/T': {'ΔH°': 0.6, 'ΔG°37': -0.42},
                'CC/G': {'ΔH°': -4.4, 'ΔG°37': -0.52},
                'CG/C': {'ΔH°': -4.0, 'ΔG°37': -0.43},
                'CT/A': {'ΔH°': -4.1, 'ΔG°37': -0.02},
                'GA/T': {'ΔH°': -1.1, 'ΔG°37': -0.62},
                'GC/G': {'ΔH°': -5.1, 'ΔG°37': -0.72},
                'GG/C': {'ΔH°': -3.9, 'ΔG°37': -0.56},
                'GT/A': {'ΔH°': -4.2, 'ΔG°37': 0.48},
                'TA/T': {'ΔH°': -6.9, 'ΔG°37': -0.71},
                'TC/G': {'ΔH°': -4.0, 'ΔG°37': -0.58},
                'TG/C': {'ΔH°': -4.9, 'ΔG°37': -0.61},
                'TT/A': {'ΔH°': -0.2, 'ΔG°37': -0.10}
            },
            '3_prime': {
                'AA/T': {'ΔH°': -0.5, 'ΔG°37': -0.12},
                'CA/G': {'ΔH°': -5.9, 'ΔG°37': -0.82},
                'GA/C': {'ΔH°': -2.1, 'ΔG°37': -0.92},
                'TA/A': {'ΔH°': -0.7, 'ΔG°37': -0.48},
                'AC/C': {'ΔH°': 4.7, 'ΔG°37': 0.28},
                'CC/T': {'ΔH°': -2.6, 'ΔG°37': -0.31},
                'GC/C': {'ΔH°': -0.2, 'ΔG°37': -0.23},
                'TC/T': {'ΔH°': 4.4, 'ΔG°37': -0.19},
                'AG/C': {'ΔH°': -4.1, 'ΔG°37': -0.01},
                'CG/T': {'ΔH°': -3.2, 'ΔG°37': -0.01},
                'GG/C': {'ΔH°': -3.9, 'ΔG°37': -0.44},
                'TG/T': {'ΔH°': -1.6, 'ΔG°37': -0.50},
                'AT/C': {'ΔH°': -3.8, 'ΔG°37': 0.13},
                'CT/T': {'ΔH°': -5.2, 'ΔG°37': -0.52},
                'GT/C': {'ΔH°': -4.4, 'ΔG°37': -0.35},
                'TT/T': {'ΔH°': 2.9, 'ΔG°37': -0.29},

            }
        }

    def get_penalty(self, end_type, base_pair):
        if end_type not in ['5_prime', '3_prime']:
            raise ValueError("End type must be '5_prime' or '3_prime'")

        if base_pair in self.penalty_values[end_type]:
            return self.penalty_values[end_type][base_pair]['ΔG°37']
        else:
            raise KeyError(f"Base pair {base_pair} not found for {end_type} in table 3")


class LoopDatabase:
    def __init__(self):
        # Table 4 data from the paper
        self.loop_params = {
            'hairpin': {
                3: {'ΔG°37': 3.5},
                4: {'ΔG°37': 3.5},
                5: {'ΔG°37': 3.3},
                6: {'ΔG°37': 4.0},
                7: {'ΔG°37': 4.2},
                8: {'ΔG°37': 4.3},
                9: {'ΔG°37': 4.5},
                10: {'ΔG°37': 4.6},
                12: {'ΔG°37': 5.0},
                14: {'ΔG°37': 5.1},
                16: {'ΔG°37': 5.3},
                18: {'ΔG°37': 5.5},
                20: {'ΔG°37': 5.7},
                25: {'ΔG°37': 6.1},
                30: {'ΔG°37': 6.3}
            },
            'bulge': {
                1: {'ΔG°37': 4.0},
                2: {'ΔG°37': 2.9},
                3: {'ΔG°37': 3.1},
                4: {'ΔG°37': 3.2},
                5: {'ΔG°37': 3.3},
                6: {'ΔG°37': 3.5},
                7: {'ΔG°37': 3.7},
                8: {'ΔG°37': 3.9},
                9: {'ΔG°37': 4.1},
                10: {'ΔG°37': 4.3},
                12: {'ΔG°37': 4.5},
                14: {'ΔG°37': 4.8},
                16: {'ΔG°37': 5.0},
                18: {'ΔG°37': 5.2},
                20: {'ΔG°37': 5.3},
                25: {'ΔG°37': 5.6},
                30: {'ΔG°37': 5.9}
            },
            'internal': {
                3: {'ΔG°37': 3.2},
                4: {'ΔG°37': 3.6},
                5: {'ΔG°37': 4.0},
                6: {'ΔG°37': 4.4},
                7: {'ΔG°37': 4.6},
                8: {'ΔG°37': 4.8},
                9: {'ΔG°37': 4.9},
                10: {'ΔG°37': 4.9},
                12: {'ΔG°37': 5.2},
                14: {'ΔG°37': 5.4},
                16: {'ΔG°37': 5.6},
                18: {'ΔG°37': 5.8},
                20: {'ΔG°37': 5.9},
                25: {'ΔG°37': 6.3},
                30: {'ΔG°37': 6.6}
            },
            'multibranched': {
                # This would need to be implemented based on the specific rules from the paper
                # The paper mentions that multiloops remain the least verified parameters
                # and that parameters for larger multiloops are calculated with Equation 7
            }
        }

    def get_value(self, loop_type, length):
        if loop_type not in self.loop_params:
            raise KeyError(f"Loop type {loop_type} not found in table 4")

        if length in self.loop_params[loop_type]:
            return self.loop_params[loop_type][length]['ΔG°37']
        else:
            # Implement extrapolation based on Equation 7 from the paper
            # This is a placeholder implementation
            x = max(self.loop_params[loop_type].keys())
            value_x = self.loop_params[loop_type][x]['ΔG°37']
            return value_x + 2.44 * 1.9872 * 310.15 * (math.log(length) - math.log(x))