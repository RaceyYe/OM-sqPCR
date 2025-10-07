"""
Primer6-Specific Thermodynamic Calculator
Implements exact ΔG and Tm calculations from Primer6's documentation
"""

import math

class primer_calcualtor:
    pass

class Primer5ThermodynamicCalculator:
    """Class for calculating thermodynamic properties of DNA sequences using Primer6's parameters"""
    __kelvin_base = 273.15
    __kelvin_standard = 273.15 + 25

    # Primer6's Nearest-Neighbor Parameters (adjusted from Allawi & SantaLucia, 1997)
    # Values pre-negated to match thermodynamic conventions
    NN_PARAMS = {
        'AA': {'ΔH': -9100, 'ΔS': -24},
        'TT': {'ΔH': -9100, 'ΔS': -24},
        'AT': {'ΔH': -8600, 'ΔS': -23.9},
        'TA': {'ΔH': -6000, 'ΔS': -16.9},
        'CA': {'ΔH': -5800, 'ΔS': -12.9},
        'TG': {'ΔH': -5800, 'ΔS': -12.9},
        'GT': {'ΔH': -6500, 'ΔS': -17.3},
        'AC': {'ΔH': -6500, 'ΔS': -17.3},
        'CT': {'ΔH': -7800, 'ΔS': -20.8},
        'AG': {'ΔH': -7800, 'ΔS': -20.8},
        'GA': {'ΔH': -5600, 'ΔS': -13.5},
        'TC': {'ΔH': -5600, 'ΔS': -13.5},
        'CG': {'ΔH': -11900, 'ΔS': -27.8},
        'GC': {'ΔH': -11100, 'ΔS': -26.7},
        'GG': {'ΔH': -11000, 'ΔS': -26.6},
        'CC': {'ΔH': -11000, 'ΔS': -26.6},
    }

    # Primer6-specific adjustments (pre-negated)
    INITIATION = {'ΔH': 0, 'ΔS': -15.1}
    SALT_ADJUSTMENT = 16.6  # log[Sodium] multiplier

    def __init__(self, default_sodium=50.0, default_primer_conc=500.0, default_temp=25.0):
        """
        Initialize the calculator with default parameters

        Args:
            default_sodium: Default sodium concentration in mM (default: 50.0)
            default_primer_conc: Default primer concentration in nM (default: 500.0)
            default_temp: Default temperature in °C (default: 25.0)
        """
        self.default_sodium = default_sodium
        self.default_primer_conc = default_primer_conc
        self.default_temp = default_temp

    def calculate_thermodynamic_values(self, seq):
        """
        Calculate the thermodynamic values (ΔH and ΔS) for a DNA sequence

        Args:
            seq: DNA sequence

        Returns:
            Tuple of (delta_h, delta_s) in cal/mol and cal/mol·K
        """
        delta_h = self.INITIATION['ΔH']
        delta_s = self.INITIATION['ΔS']

        # Nearest-neighbor contributions
        for i in range(len(seq) - 1):
            dimer = seq[i:i + 2].upper()
            delta_h += self.NN_PARAMS[dimer]['ΔH']
            delta_s += self.NN_PARAMS[dimer]['ΔS']

        return delta_h, delta_s

    def get_thermodynamic_params(self,seq, temp=None, sodium=None):
        if temp is None:
            temp = self.default_temp
        if sodium is None:
            sodium = self.default_sodium

        delta_h, delta_s = self.calculate_thermodynamic_values(seq)
        # print(f"\n delta H: {delta_h}; delta S: {delta_s} \n")
        # for calculating delta G, the initialization values should be subtracted
        delta_h += 0
        delta_s += 15.1
        return delta_h, delta_s

    def calculate_delta_g(self, seq, temp=None, sodium=None):
        """
        Calculate the Gibbs free energy (ΔG) for a DNA sequence using Primer6's formula

        Args:
            seq: DNA sequence
            temp: Temperature in °C (if None, uses default)
            sodium: Sodium concentration in mM (if None, uses default)

        Returns:
            ΔG value in kcal/mol
        """
        if temp is None:
            temp = self.default_temp
        if sodium is None:
            sodium = self.default_sodium

        delta_h, delta_s = self.calculate_thermodynamic_values(seq)
        # print(f"\n delta H: {delta_h}; delta S: {delta_s} \n")
        # for calculating delta G, the initialization values should be subtracted
        delta_h += 0
        delta_s += 15.1

        # Convert to kcal/mol
        return (delta_h  - (self.__kelvin_base + temp) * delta_s) / 1000

    def calculate_tm(self, seq, primer_conc=None, sodium=None):
        """
        Calculate the melting temperature (Tm) for a DNA sequence using Primer6's formula

        Args:
            seq: DNA sequence
            primer_conc: Primer concentration in nM (if None, uses default)
            sodium: Sodium concentration in mM (if None, uses default)

        Returns:
            Tm value in °C
        """
        if primer_conc is None:
            primer_conc = self.default_primer_conc
        if sodium is None:
            sodium = self.default_sodium

        delta_h, delta_s = self.calculate_thermodynamic_values(seq)

        try:
            return (delta_h) / (delta_s + 1.987 * math.log(250e-12 / 4)) + self.SALT_ADJUSTMENT * math.log10(
                0.205 / (1 + 0.7 * 0.205)) - 273.15
        except ZeroDivisionError:
            return float('nan')

    @staticmethod
    def reverse_complement(seq):
        """
        Get the reverse complement of a DNA sequence

        Args:
            seq: DNA sequence

        Returns:
            Reverse complement of the sequence
        """
        comp = {'A': 'T', 'T': 'A', 'C': 'G', 'G': 'C'}
        return ''.join(comp.get(base.upper(), 'N') for base in reversed(seq))

    def is_self_complementary(self, seq):
        """
        Check if a sequence is self-complementary according to Primer6's definition
        (3' 5-mer match)

        Args:
            seq: DNA sequence

        Returns:
            Boolean indicating if the sequence is self-complementary
        """
        rc = self.reverse_complement(seq)
        return any(seq[i:i + 5] == rc[-i - 5:-i] for i in range(len(seq) - 4))

    def C2K(self, celsius: float) -> float:
        """Convert Celsius to Kelvin"""
        return celsius + self.__kelvin_base

    def delta_G_calculator(self, delta_H, delta_S, temperature) -> float:
        """
        Calculate delta G from delta H and delta S

        Args:
            delta_H: Enthalpy in cal/mol
            delta_S: Entropy in cal/mol·K

        Returns:
            Delta G in kcal/mol
        """
        # the unit should all be as kcal/mol
        # according to the thermal-dynamic formula
        return (delta_H - (self.__kelvin_base + temperature) * delta_S) / 1000


def main():
    """Test the Primer6ThermodynamicCalculator with example sequences"""
    calculator = Primer5ThermodynamicCalculator()
    test_sequence = "ATCGATCG"  # Example from Primer6 docs

    delta_g = calculator.calculate_delta_g(test_sequence)
    tm = calculator.calculate_tm(test_sequence)

    print(f"Sequence: {test_sequence}")
    print(f"ΔG (Primer6): {delta_g:.2f} kcal/mol")
    print(f"Tm (Primer6): {tm:.1f} °C")

    # Additional example
    test_sequence2 = "GCGCGCGC"
    delta_g2 = calculator.calculate_delta_g(test_sequence2)
    tm2 = calculator.calculate_tm(test_sequence2)

    print(f"\nSequence: {test_sequence2}")
    print(f"ΔG (Primer6): {delta_g2:.2f} kcal/mol")
    print(f"Tm (Primer6): {tm2:.1f} °C")

    # Additional example
    test_sequence2 = "ATATATATAT" # "CATATATG"
    delta_g2 = calculator.calculate_delta_g(test_sequence2)
    tm2 = calculator.calculate_tm(test_sequence2)
    print(f"\nSequence: {test_sequence2}")
    print(f"ΔG (Primer6): {delta_g2:.2f} kcal/mol")
    print(f"Tm (Primer6): {tm2:.1f} °C")

if __name__ == "__main__":
    main()