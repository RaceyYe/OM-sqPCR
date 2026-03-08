import csv

from Bio.Seq import Seq


class helpers:


    def strand_preprocess(self, seq):
        if seq:
            if 'U' in seq:
                print("there are U, therefore it's RNA")
                return seq.upper().replace('U', 'T')
            else:
                return seq.upper()
        return seq

    def generate_cDNA(self, upper_strand_5_prime, lower_strand_5_prime, overlap=False):
        """
        根据上链和下链的 5' 端序列生成完整的 cDNA 序列。
        :param upper_strand_5_prime: 上链的 5' 端序列
        :param lower_strand_5_prime: 下链的 5' 端序列
        :param overlap: 是否有重叠部分，默认为 False
        :return: 完整的 cDNA 序列，以字典形式返回，包含上链和下链
        """
        upper_strand_5_prime = self.strand_preprocess(upper_strand_5_prime)
        lower_strand_5_prime = self.strand_preprocess(lower_strand_5_prime)
        base_pairing = {'A': 'T', 'T': 'A', 'C': 'G', 'G': 'C'}

        if overlap:
            # 找到互补部分
            for i in range(2, min(len(upper_strand_5_prime), len(lower_strand_5_prime)) + 1):
                upper_suffix = upper_strand_5_prime[-i:]
                lower_prefix = lower_strand_5_prime[:i]
                reverse_complement = ''.join([base_pairing.get(base, base) for base in upper_suffix[::-1]])
                if lower_prefix == reverse_complement:
                    # 拼接生成完整序列
                    upper = upper_strand_5_prime
                    lower = lower_strand_5_prime[::-1]  # 转换为 3' -> 5' 方向用于拼接
                    combined_upper = upper
                    combined_lower = lower[i - 1::-1] + upper_strand_5_prime[-i:][::-1] + lower[i:]
                    return {
                        "upper": combined_upper,
                        "lower": combined_lower
                    }
            print("警告：未找到足够的互补序列。")
        else:
            reverse_complement = ''.join([base_pairing.get(base, base) for base in upper_strand_5_prime[::-1]])
            if lower_strand_5_prime == reverse_complement:
                return {
                    "upper": upper_strand_5_prime,
                    "lower": lower_strand_5_prime
                }
            else:
                print("警告：下链序列与上链的反向互补序列不匹配。")

        return {
            "upper": upper_strand_5_prime,
            "lower": lower_strand_5_prime
        }

    def check_RC(self, Seq_1, Seq_2):
        """Ensure strands are in correct orientation and format"""
        strand1 = self.strand_preprocess(Seq_1)
        strand2 = self.strand_preprocess(Seq_2)

        if Seq(strand1).reverse_complement() == Seq(strand2):
            return strand1, strand2
        elif Seq(strand1[::-1]).reverse_complement() == Seq(strand2):
            print("the input two strands are complementary, meaning either one being reversed; the direction has been corrected")
            return strand1, strand2[::-1]  # make sure the sequence is from 5' to 3' for both strands
        else:
            raise ValueError("Input strands are not complementary in any orientation")

    @staticmethod
    def generate_gap_range(csv_file_path):
        """
        Generates gap range from the "FP_WOE" column in a CSV file.

        Args:
            csv_file_path (str): Path to the CSV file

        Returns:
            tuple: (min_gap, max_gap, list_of_gaps)
        """
        gap_values = []

        with open(csv_file_path, 'r') as file:
            reader = csv.DictReader(file)

            if 'FP_WOE' not in reader.fieldnames:
                raise ValueError("CSV file does not contain a 'FP_WOE' column")

            for row in reader:
                try:
                    gap = int(row['FP_WOE'])
                    gap_values.append(gap)
                except ValueError:
                    raise ValueError(f"Non-integer value in FP_WOE: {row['FP_WOE']}")

        if not gap_values:
            return 0, 0, []

        min_gap = min(gap_values)
        max_gap = max(gap_values)
        gap_range = list(range(min_gap, max_gap + 1))  # Include max_gap

        return min_gap, max_gap, gap_range  # Return 3 values

