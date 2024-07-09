import pandas as pd
import glob
import argparse
import os

def merge_tax_tables(folder_path, output_file_path):
    # Ensure the folder path ends with a slash
    if not folder_path.endswith(os.sep):
        folder_path += os.sep

    # List of TAX files in the folder
    tax_files = glob.glob(folder_path + "*.csv")
    if not tax_files:
        raise FileNotFoundError(f"No CSV files found in the directory: {folder_path}")

    # Read the first file
    first_file = tax_files[0]
    df1 = pd.read_csv(first_file)

    # Iterate over the list of TAX files, starting from the second file
    for file in tax_files[1:]:
        df_temp = pd.read_csv(file)
        new_rows_df = df_temp[~df_temp['X'].isin(df1['X'])]
        df1 = pd.concat([df1, new_rows_df], ignore_index=True)

    # Write the merged dataframe to a new CSV file
    df1.to_csv(output_file_path, index=False)

if __name__ == "__main__":
    parser = argparse.ArgumentParser(description='Merge TAX tables.')
    parser.add_argument('folder', type=str, help='Path to the folder containing the TAX tables')
    parser.add_argument('output_file', type=str, help='Path to the output CSV file')

    args = parser.parse_args()

    merge_tax_tables(args.folder, args.output_file)
