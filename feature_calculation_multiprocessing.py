#!/usr/bin/env python3

"""
example of calling the script:
python3.10 /dir/CombCalculator/feature_calculation_multiprocessing.py \
--input_csv /dir/CombCalculator/example_inputs/fc_TEST.csv \
--batch_size 100 \
--num_threads 8 \
--output_dir /dir/fc__multi_OUTPUT
"""
import os
import subprocess
import pandas as pd
from glob import glob
from multiprocessing import Pool, cpu_count
import numpy as np
import argparse

def process_batch(batch_file, script_path):
    """Process a single batch file by calling the feature calculation script."""
    batch_name = os.path.basename(batch_file)
    output1_file = os.path.join(output1_folder, batch_name.replace(".csv", "_OUTPUT.csv"))
    output2_file = os.path.join(output2_folder, batch_name.replace(".csv", "_OUTPUT_RAW.csv"))
    
    subprocess.run(["python", script_path, batch_file, output1_file, output2_file], check=True)

def split_into_batches(input_csv, batch_size, output_folder):
    """Split the input CSV into batches and save them to the output folder."""
    # Read the input CSV
    df = pd.read_csv(input_csv)

    # Split the dataframe into batches
    batches = np.array_split(df, np.ceil(len(df) / batch_size))

    # Save each batch as a new CSV file
    batch_files = []
    for i, batch in enumerate(batches):
        batch_filename = os.path.join(output_folder, f"batch_{i+1}.csv")
        batch.to_csv(batch_filename, index=False)
        batch_files.append(batch_filename)

    return batch_files

def main():
    # Argument parsing
    parser = argparse.ArgumentParser(
        description="Process a large CSV by splitting it into batches and applying feature calculations."
    )
    
    parser.add_argument(
        "--input_csv", type=str, required=True,
        help="Path to the input CSV file."
    )
    parser.add_argument(
        "--batch_size", type=int, required=True,
        help="Number of rows per batch."
    )
    parser.add_argument(
        "--num_threads", type=int, default=4,
        help="Number of threads to use for processing. Default is 4."
    )
    parser.add_argument(
        "--output_dir", type=str, required=True,
        help="Directory where the batches and output files will be saved."
    )

    args = parser.parse_args()

    # Validate inputs
    input_csv = args.input_csv
    if not os.path.isfile(input_csv):
        print(f"Error: The file {input_csv} does not exist.")
        return

    # Get the script path dynamically (feature_calculation.py should be in the same directory as this script)
    script_path = os.path.join(os.path.dirname(os.path.abspath(__file__)), "feature_calculation.py")

    # Define the output directories for batches and results
    input_folder = args.output_dir  # Set the output directory for input files and batch files
    global output1_folder, output2_folder
    output1_folder = os.path.join(input_folder, "output1")
    output2_folder = os.path.join(input_folder, "output2")

    # Ensure output folders exist
    os.makedirs(output1_folder, exist_ok=True)
    os.makedirs(output2_folder, exist_ok=True)

    # Split the input CSV into batches
    batch_files = split_into_batches(input_csv, args.batch_size, input_folder)
    
    # Set number of threads (from arguments)
    num_threads = min(args.num_threads, cpu_count())  # Limit the number of threads to the available CPUs

    # Process batch files in parallel
    with Pool(num_threads) as pool:
        pool.starmap(process_batch, [(batch_file, script_path) for batch_file in batch_files])
    
    print("All batch files processed. Now concatenating outputs...")

    # Concatenate all OUTPUT.csv files
    output1_files = sorted(glob(os.path.join(output1_folder, "*_OUTPUT.csv")))
    df_output1 = pd.concat([pd.read_csv(f) for f in output1_files], ignore_index=True)
    df_output1.to_csv(os.path.join(input_folder, "FINAL_OUTPUT.csv"), index=False)

    # Concatenate all OUTPUT_RAW.csv files
    output2_files = sorted(glob(os.path.join(output2_folder, "*_OUTPUT_RAW.csv")))
    df_output2 = pd.concat([pd.read_csv(f) for f in output2_files], ignore_index=True)
    df_output2.to_csv(os.path.join(input_folder, "FINAL_OUTPUT_RAW.csv"), index=False)

    print("Concatenation completed. FINAL_OUTPUT.csv and FINAL_OUTPUT_RAW.csv saved.")

if __name__ == "__main__":
    main()
