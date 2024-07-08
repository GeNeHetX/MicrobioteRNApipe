import pandas as pd
import glob
import argparse
import os

def merge_otu_tables(first_folder_path, second_folder_path, output_file_path):
    # Ensure the first folder path ends with a slash
    if not first_folder_path.endswith(os.sep):
        first_folder_path += os.sep

    # Ensure the second folder path ends with a slash
    if not second_folder_path.endswith(os.sep):
        second_folder_path += os.sep

    # Read the first file
    first_files = glob.glob(first_folder_path + "*.csv")
    if not first_files:
        raise FileNotFoundError(f"No CSV files found in the directory: {first_folder_path}")
    
    first_file = first_files[0]
    df1 = pd.read_csv(first_file)

    # List of OTU files to merge from the second folder
    otu_files = glob.glob(second_folder_path + "*.csv")
    if not otu_files:
        raise FileNotFoundError(f"No CSV files found in the directory: {second_folder_path}")

    # Iterate over the list of OTU files and merge each with the first file
    for file in otu_files:
        df_temp = pd.read_csv(file)
        df1 = pd.merge(df1, df_temp, on='X', how='outer')

    # Fill missing values with 0
    df1.fillna(0, inplace=True)

    # Write the merged dataframe to a new CSV file
    df1.to_csv(output_file_path, index=False)

if __name__ == "__main__":
    parser = argparse.ArgumentParser(description='Merge OTU tables.')
    parser.add_argument('first_folder', type=str, help='Path to the folder containing the first OTU table')
    parser.add_argument('second_folder', type=str, help='Path to the folder containing the rest of the OTU tables')
    parser.add_argument('output_file', type=str, help='Path to the output CSV file')

    args = parser.parse_args()

    merge_otu_tables(args.first_folder, args.second_folder, args.output_file)
