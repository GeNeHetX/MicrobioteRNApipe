import os
import csv
import argparse

def filter_and_compare(file1_path, file2_dir, output_file, threshold):
    # Read the first file and filter based on threshold
    removed_taxa = set()
    filtered_lines = []
    with open(file1_path, 'r') as file1:
        lines = file1.readlines()
        for line in lines[1:]:
            try:
                spearman_correlation = float(line.split()[1])
                if spearman_correlation < threshold:
                    removed_taxa.add(line.split()[0].strip())
                else:
                    filtered_lines.append(line)
            except ValueError:
                pass
    
    # Define the path for taxa_to_remove.csv in the output directory
    taxa_to_remove_file = os.path.join(os.path.dirname(os.path.abspath(output_file)), 'taxa_to_remove.csv')
    
    # Write removed taxon names to a CSV file in the output directory
    with open(taxa_to_remove_file, 'w', newline='') as csvfile:
        csv_writer = csv.writer(csvfile)
        csv_writer.writerow(['Taxon'])
        for taxon in removed_taxa:
            csv_writer.writerow([taxon])

    # Iterate over the second files
    for filename in os.listdir(file2_dir):
        filepath = os.path.join(file2_dir, filename)
        with open(filepath, 'r') as file2:
            lines = file2.readlines()

        # Remove rows where taxon name matches
        filtered_lines = [line for line in lines if line.split('\t')[5].strip() not in removed_taxa]

        # Write filtered lines back to the second file
        with open(filepath, 'w') as file2:
            file2.writelines(filtered_lines)

    # Write filtered lines to the output file
    with open(output_file, 'w') as output:
        output.writelines(filtered_lines)

if __name__ == "__main__":
    # Argument parser
    parser = argparse.ArgumentParser(description="Filter and compare files")
    parser.add_argument("file1", help="Path to the first file")
    parser.add_argument("file2_dir", help="Path to the directory containing second files")
    parser.add_argument("output_file", help="Path to the output file")
    parser.add_argument("--threshold", type=float, default=0.6, help="Threshold for filtering (default: 0.6)")
    args = parser.parse_args()

    # Call the function with arguments
    filter_and_compare(args.file1, args.file2_dir, args.output_file, args.threshold)

