import argparse
from Bio.Seq import Seq
import sys
sys.path.append("/home/racey/Desktop/primerDesign")
from modules.auxiliary_modules.helper_functions import helpers
helper = helpers()

class cDNA_generator:
    '''
    All input and output (except the printed sequence) should be in 5' -> 3'
    '''

    def __init__(self, sRNA, stem_loop):
        self.sRNA = helper.strand_preprocess(sRNA)
        self.stem_loop = helper.strand_preprocess(stem_loop)
        self.cDNA = self.generate_cDNA()

    def generate_cDNA(self):
        # Convert sRNA to DNA (replace U with T if necessary)
        sRNA_dna = self.sRNA.replace('U', 'T')
        upper_strand = Seq(sRNA_dna) + Seq(self.stem_loop).reverse_complement()

        # Generate lower strand using Biopython's reverse_complement method
        lower_strand = upper_strand.reverse_complement()

        return (upper_strand, lower_strand)

    def display_cDNA(self):
        upper, lower = self.cDNA
        print(f"forward strand: 5' {upper} 3'")
        print(f"reverse strand: 3' {lower[::-1]} 5'")

    def fetch_sequence(self, strand: str = "both"):
        upper, lower = self.cDNA  # Need to unpack the cDNA tuple first
        if strand == "both":
            return {"upper": upper, "lower": lower}
        elif strand == "upper":
            return {"upper": upper}
        elif strand == "lower":
            return {"lower": lower}
        else:
            return None


def main():
    parser = argparse.ArgumentParser(description='cDNA generator based on sRNA and stemloop sequence')
    parser.add_argument("-s", "--sRNA", required=True, help='sRNA sequence (containing "U" except none is)')
    parser.add_argument('-sl', '--stemloop', required=True, help='stemloop primer sequence')

    args = parser.parse_args()

    designer = cDNA_generator(args.sRNA, args.stemloop)
    designer.display_cDNA()
    fs = designer.fetch_sequence()
    print(fs)
    upper = fs.get("upper","Unknown")
    lower = fs.get("lower","Unknown")

    print(Seq(upper).reverse_complement()==Seq(lower))

if __name__ == "__main__":
    main()