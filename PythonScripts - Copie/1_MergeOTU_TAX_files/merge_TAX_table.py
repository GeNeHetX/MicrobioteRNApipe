import pandas as pd
import glob
import argparse
import os

def merge_tax_tables(first_folder_path, second_folder_path, output_file_path):
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

    # List of TAX files to merge from the second folder
    tax_files = glob.glob(second_folder_path + "*.csv")
    if not tax_files:
        raise FileNotFoundError(f"No CSV files found in the directory: {second_folder_path}")

    # Iterate over the list of TAX files and merge each with the first file
    for file in tax_files:
        df_temp = pd.read_csv(file)
        new_rows_df = df_temp[~df_temp['X'].isin(df1['X'])]
        df1 = pd.concat([df1, new_rows_df])

    # Write the merged dataframe to a new CSV file
    df1.to_csv(output_file_path, index=False)

if __name__ == "__main__":
    parser = argparse.ArgumentParser(description='Merge TAX tables.')
    parser.add_argument('first_folder', type=str, help='Path to the folder containing the first TAX table')
    parser.add_argument('second_folder', type=str, help='Path to the folder containing the rest of the TAX tables')
    parser.add_argument('output_file', type=str, help='Path to the output CSV file')

    args = parser.parse_args()

    merge_tax_tables(args.first_folder, args.second_folder, args.output_file)
