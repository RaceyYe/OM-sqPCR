import pandas as pd
import sys

if __name__ == "__main__":
    calibration_file = pd.read_csv(sys.argv[1], usecols=["Model", "Threshold"])
    output_file = sys.argv[2]
    calibration_file.to_csv(output_file, index=False)

