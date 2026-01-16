#/usr/bin/env python3

# this tool is used to combine multiple CSV files (replicates) into a single CSV file
# it looks at all the files matching the pattern rep*.csv in the current directory, combining the recorded work values
# from the 3rd column of each file, using the first column as the index.
# the output file is named <prefix>_combined.csv, where <prefix> is provided as a command-line argument.


import sys
import numpy as np
import glob
import pandas as pd
import os

def main():
    # Check for argument
    if len(sys.argv) != 2:
        print("Usage: python combine_reps.py <prefix>")
        print("Example: python combine_reps.py initial_attempt")
        sys.exit(1)

    prefix = sys.argv[1]

    # Find matching files
    pattern = f"rep*.csv"
    files = sorted(glob.glob(pattern))

    if not files:
        print(f"No files found matching pattern '{pattern}'")
        sys.exit(1)

    combined_df = None

    for file in files:
        # Read file and extract first (index) and 3rd column
        df = pd.read_csv(file)

        # Check for at least 3 columns
        if df.shape[1] < 3:
            print(f"Warning: {file} has fewer than 3 columns — skipping")
            continue

        # Take index (col 0) and needed column (col 2)
        temp = df.iloc[:, [0, 2]].copy()

        # Rename the value column based on the file (e.g., rep1_1.5.csv → rep1)
        base_name = os.path.basename(file).split("_")[0]
        temp.columns = ["index", base_name]

        # Merge
        if combined_df is None:
            combined_df = temp
        else:
            combined_df = pd.merge(combined_df, temp, on="index", how="outer")

    if combined_df is None:
        print("No valid data to combine.")
        sys.exit(1)

    # Save combined CSV
    output_filename = f"{prefix}_combined.csv"
    combined_df = combined_df.replace('NAN', np.nan)
    combined_df = combined_df.drop(columns=['index'])
    combined_df = combined_df.dropna(how='all')
    combined_df.to_csv(output_filename, index=False)
    print(f"✅ Combined file saved as: {output_filename}")

if __name__ == "__main__":
    main()

