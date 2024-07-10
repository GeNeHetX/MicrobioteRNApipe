import os
import pandas as pd
import argparse


def merge_kreport_files(kreport_directory, output_file, taxonomic_level):
    # List all the files ending with "kreport" in the specified directory
    file_list = [os.path.join(kreport_directory, f) for f in os.listdir(kreport_directory) if f.endswith("kreport")]

    # Initialize an empty dictionary to store data from each file
    data_dict = {}

    # Step 2: Read each file, filter it, and store data in the dictionary
    for file in file_list:
        # Extract sample name from file name
        sample_name = os.path.basename(file).split('.')[0]
        
        # Read file (assuming no headers and tab-separated)
        data = pd.read_csv(file, sep='\t', header=None)
        
        # Filter the data on the specified taxonomic level ("G" or "S")
        filtered_data = data[(data[5] == taxonomic_level)]
        
        # Remove spaces from the 7th column
        filtered_data[7] = filtered_data[7].str.replace(' ', '')
        
        # Extract the second column (abundance) and use the 7th column (taxonomic ID) as the index
        filtered_abundance = filtered_data.set_index(7)[2]
        
        # Store data in dictionary
        data_dict[sample_name] = filtered_abundance

    # Step 3: Create DataFrame from the dictionary
    merged_data = pd.DataFrame(data_dict)
    
    # Rename the index based on the taxonomic level
    index_name = "Genus" if taxonomic_level == "G" else "Species"
    merged_data.index.name = index_name
    
    # Fill NaN values with 0
    merged_data = merged_data.fillna(0)
    
    # Write the merged data to a file
    merged_data.to_csv(output_file)
    print(f"Merging complete. Output saved to {output_file}")

if __name__ == "__main__":
    # Set up argument parsing
    parser = argparse.ArgumentParser(description="Merge kreport files by taxonomic level")
    parser.add_argument("-i", "--input", required=True, help="Input directory containing kreport files")
    parser.add_argument("-o", "--output", required=True, help="Output file path")
    parser.add_argument("-t", "--taxonomic_level", required=True, choices=['G', 'S'], help="Taxonomic level to filter by (G for Genus, S for Species)")

    # Parse arguments
    args = parser.parse_args()

    # Call the function with parsed arguments
    merge_kreport_files(args.input, args.output, args.taxonomic_level)