#!/usr/bin/env python3

"""

examples of calling the script:
1. FULL UNIPROT
python3.10 /dir/CombCalculator/calculate_uniprot.py \
--mode fetch_full \
--output_file /dir/mapping_OUTPUT/uniprot_mapped_disorder_TEST.csv \
--cached_disorder_file /dir/CombCalculator/Datasets/uniprot_mapped_disorder.csv

2. SAMPLE UNIPROT
python3.10 /dir/CombCalculator/calculate_uniprot.py \
--mode fetch_sample \
--sample_size 20 \
--output_file /dir/mapping_OUTPUT/uniprot_mapped_disorder_TEST.csv \
--cached_disorder_file /dir/CombCalculator/Datasets/uniprot_mapped_disorder.csv 

3. LIST OF UIDs FOR MAPPING
python3.10 /dir/CombCalculator/calculate_uniprot.py \
--mode map_existing \
--input_file /dir/CombCalculator/example_inputs/uid_list_TEST.csv \
--output_file /dir/mapping_OUTPUT/uniprot_mapped_disorder_TEST.csv \
--cached_disorder_file /dir/CombCalculator/Datasets/uniprot_mapped_disorder.csv

"""
import os
import sys
sys.path.append(os.path.join(os.path.dirname(__file__), 'iupred2a'))
from iupred2a import iupred2a_lib
import pandas as pd
import numpy as np
from unipressed import UniprotkbClient
from tqdm import tqdm
import re
import argparse
from pathlib import Path

# Function to parse command-line arguments
def parse_args():
    parser = argparse.ArgumentParser(description="UniProt Mapping & Disorder Prediction")
    parser.add_argument('--mode', choices=['fetch_full', 'fetch_sample', 'map_existing'], required=True, help="Mode of operation.")
    parser.add_argument('--sample_size', type=int, default=10, help="Sample size if mode is 'fetch_sample'.")
    parser.add_argument('--input_file', type=str, help="Input CSV file for 'map_existing' mode. (Should contain a single 'uid' column of UniProt IDs)")
    parser.add_argument('--output_file', type=str, required=True, help="Output CSV file.")
    parser.add_argument('--cached_disorder_file', type=str, required=True, help="File with disorder data.")

    return parser.parse_args()

# Fetch data from UniProtKB
def fetch_uniprot_data(organism_id="9606", sample_size=None):
    data = {'uid': [], 'refseq_id': [], 'seq': []}
    track = 0
    all_records = UniprotkbClient.search(query={"organism_id": organism_id}).each_record()
    print("Fetching UniProt entries...")

    all_entries = []  # List to store all fetched entries

    for record in all_records:
        if (record.get('entryType') == 'UniProtKB reviewed (Swiss-Prot)') and ('-' not in record.get('primaryAccession')):
            uid_value = record.get('primaryAccession')
            crosrefs = record.get('uniProtKBCrossReferences', [])
            refseq_list = [i['id'] for i in crosrefs if i['database'] == 'RefSeq']
            refseq_id = re.sub(r'\.\d+', '', refseq_list[0]) if refseq_list else np.nan
            sequence_field = record.get('sequence', {})
            seq_key = sequence_field.get('value', '')

            all_entries.append({'uid': uid_value, 'refseq_id': refseq_id, 'seq': seq_key})
            track += 1
            print(f"Uniprot IDs fetched: {track} / ~20,4k")

    # Convert the list of entries to a DataFrame
    df = pd.DataFrame(all_entries)
    print(f"Fetched {len(df)} UniProt entries.")
    
    # Randomly sample from the fetched entries if sample_size is specified
    if sample_size:
        df = df.sample(n=sample_size, random_state=42).reset_index(drop=True)  # random_state for reproducibility
    
    return df

# Map existing dataset with UniProt IDs and include all relevant columns
def map_existing_uniprot_ids(input_df, cached_disorder_file):
    if os.path.exists(cached_disorder_file):
        cached_df = pd.read_csv(cached_disorder_file)
    else:
        cached_df = pd.DataFrame(columns=['uid', 'refseq_id', 'seq', 'disorder'])

    # Merge the input data with the cached data based on 'uid'
    input_df = input_df.copy()
    input_df = input_df.merge(cached_df, on='uid', how='left')

    # Identify missing UIDs (those without disorder data)
    missing_uids = input_df[input_df['disorder'].isna()]['uid']
    return input_df, missing_uids

# Calculate overall disorder score
def calculate_overall_disorder(iupred_result):
    disorder_scores = iupred_result[0]
    return sum(disorder_scores) / len(disorder_scores)

# Run disorder predictions and map results
def run_disorder_predictions(df, cached_disorder_file):
    if os.path.exists(cached_disorder_file):
        cached_df = pd.read_csv(cached_disorder_file)
    else:
        cached_df = pd.DataFrame(columns=['uid', 'disorder'])

    disorder_scores = []
    for _, row in tqdm(df.iterrows(), total=len(df), desc="Predicting Disorder", unit="sequence"):
        uid, seq = row['uid'], row['seq']
        cached_row = cached_df[cached_df['uid'] == uid]
        if not cached_row.empty:
            overall_disorder_score = cached_row.iloc[0]['disorder']
            # Check if the cached value is NaN or not a valid number
            if pd.isna(overall_disorder_score) or not isinstance(overall_disorder_score, (int, float)):
                print(f"{uid} - Cached disorder score is invalid. Recalculating...")
                iupred2_result = iupred2a_lib.iupred(seq, mode="long")
                overall_disorder_score = calculate_overall_disorder(iupred2_result)
                print(f"{uid} - Recalculated Disorder Score: {overall_disorder_score:.4f}")

                # Update the cache
                cached_df = pd.concat([
                    cached_df,
                    pd.DataFrame({'uid': [uid], 'disorder': [overall_disorder_score]})
                ], ignore_index=True)
        else:
            # If no cached value, calculate disorder score
            iupred2_result = iupred2a_lib.iupred(seq, mode="long")
            overall_disorder_score = calculate_overall_disorder(iupred2_result)
            print(f"{uid} - New Disorder Score: {overall_disorder_score:.4f}")

            # Add to cache
            cached_df = pd.concat([
                cached_df,
                pd.DataFrame({'uid': [uid], 'disorder': [overall_disorder_score]})
            ], ignore_index=True)

        disorder_scores.append(overall_disorder_score)

    # Save updated cache
    cached_df.to_csv(cached_disorder_file, index=False)
    print(f"Updated cache saved to: {cached_disorder_file}")

    df['disorder'] = disorder_scores
    return df

# Main Execution Flow
if __name__ == "__main__":
    args = parse_args()

    mode = args.mode
    sample_size = args.sample_size
    input_file = args.input_file
    output_final_csv = args.output_file
    cached_disorder_file = args.cached_disorder_file

    # Ensure output directory exists
    output_path = Path(output_final_csv)
    output_path.parent.mkdir(parents=True, exist_ok=True)

    if mode == "fetch_full":
        uniprot_df = fetch_uniprot_data()
        uniprot_df = run_disorder_predictions(uniprot_df, cached_disorder_file)
    elif mode == "fetch_sample":
        uniprot_df = fetch_uniprot_data(sample_size=sample_size)
        uniprot_df = run_disorder_predictions(uniprot_df, cached_disorder_file)
    elif mode == "map_existing":
        input_df = pd.read_csv(input_file)
        uniprot_df, missing_uids = map_existing_uniprot_ids(input_df, cached_disorder_file)

        if not missing_uids.empty:
            print(f"Found {len(missing_uids)} UIDs without cached disorder scores. Running predictions...")
            missing_data = fetch_uniprot_data()
            missing_data = missing_data[missing_data['uid'].isin(missing_uids)]
            missing_data = run_disorder_predictions(missing_data, cached_disorder_file)
            uniprot_df = pd.concat([uniprot_df, missing_data], ignore_index=True)
    else:
        raise ValueError("Invalid mode. Set mode to 'fetch_full', 'fetch_sample', or 'map_existing'.")

    # Save the final CSV
    uniprot_df.to_csv(output_final_csv, index=False)
    print("Final CSV with disorder values saved to:", output_final_csv)