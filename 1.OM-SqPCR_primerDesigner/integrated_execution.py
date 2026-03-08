from modules.FP_design import FP_designer
from modules.auxiliary_modules.forwardPrimer_extension import *
from modules.RP_design_v2_copy import RP_designer
from modules.validityCheck import interaction_check
from modules.pairing import generate_combinations
import os
import argparse
import pandas


class execution:

    def __init__(self, args, task="", criteria=None, GF_Seq=None):
        # Preserve original args immutably
        self.args = argparse.Namespace(**vars(args))  # Create a copy
        # self.args = args # nested directories would generate rather than parallel
        self.prefix = task
        # Correct output directory setup
        self.output_dir = os.path.join(args.output, task)
        os.makedirs(self.output_dir, exist_ok=True)
        # Ensure args use the new output directory for this task
        self.args.output = self.output_dir
        self.criteria = criteria
        self.GF_Seq = GF_Seq


    def __FPDesign(self):
        # Instantiate FP_designer with extracted u and l, and other parameters from args
        designer = FP_designer(
            u=self.args.upper_strands,
            l=self.args.lower_strands,
            o=self.args.output,
            b=self.args.tm_base,
            v=self.args.tm_var,
            # b1=self.args.tm_base_probe_1,
            # v1=self.args.tm_var_probe_1,
            dg_min=self.args.dg_min,
            dg_max=self.args.dg_max,
            temp=self.args.temp,
            min_5ext=self.args.min_5ext,
            max_5ext=self.args.max_5ext
        )

        # Proceed with primer design
        try:
            valid_primers = designer.design_primers(self.args.output, self.args.temp)
            print(f"Forward primers designed successfully for {self.prefix}")
        except ValueError as e:
            print(f"Forward primer design failed for {self.prefix}: {str(e)}")
            return 1

    # def   __FPDesign(self):
    #     valid_primers = optimize_tm_with_iteration(
    #         self.args.upper_strands,
    #         # stemloop_seq,
    #         Tm_min=57,
    #         Tm_max=61,
    #         delta_G_min=-4.5,
    #         delta_G_max=10,
    #         temperature=temperature
    #     )
    #     # pass

    def __RPDesign(self):
        designer = RP_designer(
            selected_seq=self.GF_Seq,
            u=self.args.upper_strands,
            l=self.args.lower_strands,
            o=self.args.output,
            b=self.args.tm_base,
            v=self.args.tm_var,
            b1=self.args.tm_base_probe_1,
            v1=self.args.tm_var_probe_1,
            b2=self.args.tm_base_probe_2,
            v2=self.args.tm_var_probe_2,
            dg_min=self.args.dg_min,
            dg_max=self.args.dg_max,
            temp=self.args.temp,
            GAP=self.args.GAP,
            PR_gap_max=self.args.PR_gap_max,
            min_RP_len=self.args.min_RP_len,
            max_RP_len=self.args.max_RP_len,
            min_probe_len=self.args.min_probe_len,
            max_probe_len=self.args.max_probe_len
        )

        # Generate all possible pairs
        try:
            result_df = designer.design_probes()
            final_df = designer.secondary_structure_check(result_df)
            print(f"Reverse primers and probes designed successfully for {self.prefix}")
            print(final_df.columns)
        except Exception as e:
            print(f"Reverse primer/probe design failed for {self.prefix}: {str(e)}")
            return 1

    # def __interaction(self):
    #     try:
    #         # Check if prerequisite files exist
    #         fp_file = os.path.join(self.args.output, "designed_FP.csv")
    #         rp_file = os.path.join(self.args.output, "designed_RP_PC.csv")
    #         if self.args.output.startWith("Priority"):
    #             task_type="type1"
    #         else:
    #             task_type="type2"


    def __interaction(self):
        try:
            # Paths to prerequisite files
            fp_file = os.path.join(self.args.output, "designed_FP.csv")
            rp_file = os.path.join(self.args.output, "designed_RP_PC.csv")

            # Determine task type by folder name
            folder = os.path.basename(self.args.output.rstrip(os.sep))
            if folder.startswith("Priority"):
                task_type = "type1"
            else:
                task_type = "type2"

            # Ensure output is indeed a directory
            if not os.path.isdir(self.args.output):
                raise NotADirectoryError(f"Expected directory, got: {self.args.output}")

            # Check files exist
            if not os.path.exists(fp_file):
                raise FileNotFoundError(f"Forward primer file not found: {fp_file}")
            if not os.path.exists(rp_file):
                raise FileNotFoundError(f"Reverse primer file not found: {rp_file}")

            # Run the interaction check
            interaction_check(
                task_type=task_type,
                criteria=self.criteria,
                delta_G_min=self.args.dg_min,
                delta_G_max=self.args.dg_max,
                temperature=self.args.temp,
                delta_G_FR_min=self.args.dg_min_FR,
                delta_G_FR_max=self.args.dg_max_FR,
                delta_G_PO_min=self.args.dg_min_PP,
                delta_G_PO_max=self.args.dg_max_PP,
                directory=self.args.output
            )
            print(f"Interaction check completed for {self.prefix}")
            return 0

        except FileNotFoundError as e:
            print(f"[ERROR] File error during interaction check for {self.prefix}: {e}")
            return 1

        except Exception as e:
            print(f"[ERROR] Interaction check failed for {self.prefix}: {e}")
            import traceback;
            traceback.print_exc()
            return 1

    def run_design(self):
        """Execute all design steps in sequence"""
        print(f"\nStarting design process for {self.prefix}")
        self.__FPDesign()
        self.__RPDesign()
        self.__interaction()
        print(f"Completed design process for {self.prefix}\n")


def main():
    parser = argparse.ArgumentParser(description='Design PCR primers and probes')
    parser.add_argument('-u', '--upper_strands', type=str, required=True,
                        help='Upper strand of the cDNA (5\'->3\')')
    parser.add_argument('-l', '--lower_strands', type=str, required=True,
                        help='Lower strand of the cDNA (5\'->3\')')
    parser.add_argument('-o', '--output', required=True, help='Output directory')
    parser.add_argument('-b', '--tm_base', type=float, default=59, help='Base TM for variance')
    parser.add_argument('-v', '--tm_var', type=float, default=2, help='TM variance')
    parser.add_argument('-b1', '--tm_base_probe_1', type=float, default=58.5, help='probe Base TM for variance')
    parser.add_argument('-v1', '--tm_var_probe_1', type=float, default=2, help='probe TM variance')
    parser.add_argument('-b2', '--tm_base_probe_2', type=float, default=71, help='probe Base TM for variance (option 2)')
    parser.add_argument('-v2', '--tm_var_probe_2', type=float, default=2, help='probe TM variance (option 2)')
    parser.add_argument('--dg_min', type=float, default=-4.5, help='Minimum ΔG')
    parser.add_argument('--dg_max', type=float, default=4.5, help='Maximum ΔG')
    parser.add_argument('--temp', type=float, default=37, help='Reaction temperature')
    parser.add_argument('--min_5ext', type=int, default=6, help='Minimum 5\' extension length')
    parser.add_argument('--max_5ext', type=int, default=6, help='Maximum 5\' extension length')
    parser.add_argument('--GAP', type=int, default=25,
                        help='Minimum gap between forward and reverse primers, ensuring the space for probe')
    parser.add_argument('--PR_gap_max', type=int, default=6,
                        help='Maximum gap between reverse primer and probe carrier')
    parser.add_argument('--min_RP_len', type=int, default=18, help='Minimum reverse primer length')
    parser.add_argument('--max_RP_len', type=int, default=26, help='Maximum reverse primer length')
    parser.add_argument('--min_probe_len', type=int, default=18, help='Minimum probe carrier length')
    parser.add_argument('--max_probe_len', type=int, default=26, help='Maximum probe carrier length')

    args = parser.parse_args()


    # Set default output if not provided
    if not args.output:
        args.output ="./Combination_execution"


def run_directly(output_dir, generator=False):
    info_path = os.path.join(output_dir, "tasks_info")
    os.makedirs(info_path,exist_ok=True)
    """Run the script with hardcoded values without requiring command line arguments"""
    if not generator:
        """Run the script with hardcoded values without requiring command line arguments"""
        tasks = {
            "Task1": {
                "upper": "UGUCAGUUUGUCAAAUACCCCA",
                "lower": "GTTGGCTCTGGTGCAGGGTCCGAGGTATTCGCACCAGAGCCAAC"
            },
            # "Task2": {
            #     "upper": "UCCUCGUUAGUAUAGUGGUGAGUAUCCCCGCC",
            #     "lower": "GTCGTATCCAGTGCGTGTCGTGGAGTCGGCAATTGCACTGGATACGAC"
            # },
            # "Task3": {
            #     "upper": "GCGACCUCAGAUCAGACGUGGCGACCCGCUGAAUUU",
            #     "lower": "CTCAACTGGTGTCGTGGAGTCAGTGCAATTCCAGTTGAG"
            # },
            # "Task4": {
            #     "upper": "UGAGGUAGUAGAUUGUAUAGUU",
            #     "lower": "GTATGTGGGAGCCCACACTCTGGTCGACAGATACGAATATCTGGACCCGACCGCTCCCACATAC"
            # }
        }
    else:
        Seqs = {
            'sRNAs': {
                "UGUCAGUUUGUCAAAUACCCCA",
                # "UCCUCGUUAGUAUAGUGGUGAGUAUCCCCGCC",
                # "GCGACCUCAGAUCAGACGUGGCGACCCGCUGAAUUU",
                # "UGAGGUAGUAGAUUGUAUAGUU"
            },
            'stemloops': {
                "GTTGGCTCTGGTGCAGGGTCCGAGGTATTCGCACCAGAGCCAAC",
                "GTCGTATCCAGTGCGTGTCGTGGAGTCGGCAATTGCACTGGATACGAC",
                "CTCAACTGGTGTCGTGGAGTCAGTGCAATTCCAGTTGAG",
                "GTATGTGGGAGCCCACACTCTGGTCGACAGATACGAATATCTGGACCCGACCGCTCCCACATAC"
            }
        }

        result_df = generate_combinations(Seqs.get('sRNAs', set()), Seqs.get("stemloops", set()))
        result_df.to_csv(os.path.join(info_path,"combinations.csv"), index=False)

        # Display results
        print(f"Generated {len(result_df)} combinations")

        # Convert result_df to tasks dictionary
        tasks = {
            f"Task{row['Task']}": {
                "upper": row['sRNA_sequence'],
                "lower": row['Stemloop_sequence']
            }
            for _, row in result_df.iterrows()
        }

    # # Print the resulting dictionary
    # print(tasks)

    class Args:
        pass

    args = Args()
    args.output = "/home/racey/Desktop/primerDesign_Outputs/Combination_execution"
    args.tm_base = 59
    args.tm_var = 2
    args.tm_base_probe_1 = 58.5
    args.tm_var_probe_1 = 2
    args.tm_base_probe_2 = 71
    args.tm_var_probe_2 = 2
    args.dg_min = -4.5
    args.dg_max = 4.5
    args.temp = 37
    args.min_5ext = 6
    args.max_5ext = 6
    args.GAP = 25
    args.PR_gap_max = 6
    args.min_RP_len = 18
    args.max_RP_len = 22
    args.min_probe_len = 18
    args.max_probe_len = 26
    args.dg_min_FR = -5.5
    args.dg_max_FR = 5.5
    args.dg_min_PP = -7.5
    args.dg_max_PP = 7.5

    for task_name, task_data in tasks.items():
        args.upper_strands = task_data["upper"]
        args.lower_strands = task_data["lower"]
        executor = execution(args, task_name)
        executor.run_design()


if __name__ == '__main__':
    import sys
    if len(sys.argv) > 1:  # If command-line args provided
        main()
    else:  # No arguments -> use hardcoded values
        run_directly(output_dir="/home/racey/Desktop/primerDesign/Outputs/Combination_execution")