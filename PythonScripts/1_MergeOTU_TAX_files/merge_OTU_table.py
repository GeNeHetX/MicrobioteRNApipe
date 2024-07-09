import pandas as pd
import glob
import argparse
import os

def merge_otu_tables(folder_path, output_file_path):
    # Ensure the folder path ends with a slash
    if not folder_path.endswith(os.sep):
        folder_path += os.sep

    # List of OTU files in the folder
    otu_files = glob.glob(folder_path + "*.csv")
    if not otu_files:
        raise FileNotFoundError(f"No CSV files found in the directory: {folder_path}")

    # Read the first file
    first_file = otu_files[0]
    df1 = pd.read_csv(first_file)

    # Iterate over the list of OTU files, starting from the second file
    for file in otu_files[1:]:
        df_temp = pd.read_csv(file)
        df1 = pd.merge(df1, df_temp, on='X', how='outer')

    # Fill missing values with 0
    df1.fillna(0, inplace=True)

    # Write the merged dataframe to a new CSV file
    df1.to_csv(output_file_path, index=False)

if __name__ == "__main__":
    parser = argparse.ArgumentParser(description='Merge OTU tables.')
    parser.add_argument('folder', type=str, help='Path to the folder containing the OTU tables')
    parser.add_argument('output_file', type=str, help='Path to the output CSV file')

    args = parser.parse_args()

    merge_otu_tables(args.folder, args.output_file)
