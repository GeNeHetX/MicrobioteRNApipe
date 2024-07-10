import os
import csv
from collections import defaultdict
import matplotlib.pyplot as plt
import pandas as pd
import argparse

# Function to read data from a kreport file and accumulate reads count for each taxon
def read_kreport_file(file_path, taxonomic_level):
    taxon_counts = defaultdict(int)
    with open(file_path, 'r') as file:
        reader = csv.reader(file, delimiter='\t')
        for row in reader:
            # Filter by the specified taxonomic level
            if row[5] == taxonomic_level:
                taxon = row[7].strip().replace(' ', '')  # Remove leading/trailing white spaces and replace spaces with underscores
                reads_count = int(row[2])
                taxon_counts[taxon] += reads_count
    return taxon_counts

# Function to calculate mean abundance for taxons across files with the same prefix
def calculate_mean_abundance(input_folder, file_prefixes, taxonomic_level):
    total_taxon_counts = defaultdict(int)
    file_count = 0

    for prefix in file_prefixes:
        # Accumulate reads count for each taxon across files with the same prefix
        for filename in os.listdir(input_folder):
            if filename.startswith(prefix) and filename.endswith("kreport"):
                file_path = os.path.join(input_folder, filename)
                taxon_counts = read_kreport_file(file_path, taxonomic_level)
                for taxon, count in taxon_counts.items():
                    total_taxon_counts[taxon] += count
                file_count += 1

    if file_count == 0:
        print("No files found matching the specified prefixes.")
        return None

    # Calculate mean abundance
    mean_abundance = {taxon: count / file_count for taxon, count in total_taxon_counts.items()}
    return mean_abundance

def generate_pie_chart(abundance_percentages, title, output_folder, top_n):
    # Sort taxa based on abundance and consider only the top N
    sorted_abundance = dict(sorted(abundance_percentages.items(), key=lambda item: item[1], reverse=True)[:top_n])
    
    # Calculate total abundance of taxa not in the top N
    others_abundance = sum(abundance_percentages.values()) - sum(sorted_abundance.values())
    
    # Add "Others" category to the sorted abundance dictionary 
    sorted_abundance["Others"] = others_abundance

    labels = list(sorted_abundance.keys())
    sizes = list(sorted_abundance.values())

    plt.figure(figsize=(10, 10))
    patches, texts, autotexts = plt.pie(sizes, labels=labels, autopct='%1.1f%%', startangle=140, shadow=True)

    # Customize text properties
    for text in texts:
        text.set_fontsize(10)
        text.set_fontweight('bold')

    for autotext in autotexts:
        autotext.set_fontsize(10)
        autotext.set_fontweight('bold')

    plt.axis('equal')  # Equal aspect ratio ensures that pie is drawn as a circle.
    plt.title(f"{title}")
    plt.tight_layout()  # Ensures all labels are visible

    # Print labels and percentages
    print(f"Labels for {title}:")
    for label, autotext in zip(labels, autotexts):
        percentage = autotext.get_text()
        print(f"{label:<20} {percentage:>10}")

    plt.savefig(os.path.join(output_folder, f"{title.replace(' ', '')}_pie_chart.png"), dpi=300)
    plt.close()  # Close the plot to release memory

    # Plot sorted abundance
    plt.figure(figsize=(10, 6))
    plt.barh(list(sorted_abundance.keys()), list(sorted_abundance.values()))
    plt.xlabel('Abundance')
    plt.ylabel('Taxon')
    plt.title(f'Sorted Abundance for {title}')
    plt.tight_layout()
    plt.savefig(os.path.join(output_folder, f"{title.replace(' ', '')}_sorted_abundance.png"), dpi=300)
    plt.close()

# Main function
def main():
    parser = argparse.ArgumentParser(description="Merge kreport files by taxonomic level and generate pie charts")
    parser.add_argument("-i", "--input", required=True, help="Input directory containing kreport files")
    parser.add_argument("-o", "--output", required=True, help="Output folder for pie charts and CSV file")
    parser.add_argument("-t", "--taxonomic_level", required=True, choices=['G', 'S'], help="Taxonomic level to filter by (G for Genus, S for Species)")
    parser.add_argument("-n", "--top_n", type=int, default=15, help="Number of top taxa to include in the pie chart")

    args = parser.parse_args()

    input_folder = args.input
    output_folder = args.output
    taxonomic_level = args.taxonomic_level
    top_n = args.top_n

    os.makedirs(output_folder, exist_ok=True)

    # Determine the prefixes from the files in the input folder
    prefixes = set()
    for filename in os.listdir(input_folder):
        if filename.endswith("kreport"):
            prefix = filename.split('_')[0]
            prefixes.add(prefix)

    combined_abundance = defaultdict(int)
    total_files = 0

    for prefix in prefixes:
        abundance = calculate_mean_abundance(input_folder, [prefix], taxonomic_level)

        if abundance:
            total_abundance = sum(abundance.values())
            abundance_percentage = {taxon: (count / total_abundance) * 100 for taxon, count in abundance.items()}
            
            # Generate pie chart for each prefix
            generate_pie_chart(abundance_percentage, f"Mean Abundance for {prefix}", output_folder, top_n)
            
            # Save abundance percentages to CSV
            percentage_df = pd.DataFrame(list(abundance_percentage.items()), columns=['Taxon', 'Abundance'])
            percentage_df.to_csv(os.path.join(output_folder, f"{prefix}_abundance_percentage.csv"), sep='\t', index=False)

            # Combine abundances for overall summary
            for taxon, count in abundance.items():
                combined_abundance[taxon] += count
            total_files += 1

    if total_files > 0:
        combined_abundance_percentage = {taxon: (count / total_files) for taxon, count in combined_abundance.items()}
        combined_total_abundance = sum(combined_abundance_percentage.values())
        combined_abundance_percentage = {taxon: (count / combined_total_abundance) * 100 for taxon, count in combined_abundance_percentage.items()}
        
        # Generate pie chart for the combined data
        generate_pie_chart(combined_abundance_percentage, "Mean Abundance for Combined Data", output_folder, top_n)
        
        # Save combined abundance percentages to CSV
        combined_percentage_df = pd.DataFrame(list(combined_abundance_percentage.items()), columns=['Taxon', 'Abundance'])
        combined_percentage_df.to_csv(os.path.join(output_folder, "combined_abundance_percentage.csv"), sep='\t', index=False)

if __name__ == "__main__":
    main()
