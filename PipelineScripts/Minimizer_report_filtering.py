import os
import sys

def process_file(input_file, output_file, threshold):
    # Initialize variables to keep track of whether to keep lines or not
    keep_lines = False
    filtered_lines = []

    # Open the input file and read all lines
    with open(input_file, 'r') as infile:
        lines = infile.readlines()

    # Loop through each line in the input file
    for line in lines:
        # Split the line into columns
        columns = line.strip().split('\t')

        # Check if it's a D entry and the 7th column is 'Bacteria'
        if columns[5] == "D" and "Bacteria" in columns[7]:
            keep_lines = True

        # Check if it's a D entry and the 7th column is not 'Bacteria'
        if keep_lines and columns[5] == 'D' and "Bacteria" not in columns[7]:
            keep_lines = False

        # If keep_lines is True, further process the line
        if keep_lines:
            # Check if the 2nd column value meets the threshold criteria
            try:
                if float(columns[2]) < threshold:
                    continue
            except ValueError:
                # If the value in the second column cannot be converted to a float, skip this line
                continue
            
            # Add the line to the filtered lines list
            filtered_lines.append(line)

    # Write the filtered lines to the output file
    with open(output_file, 'w') as outfile:
        for line in filtered_lines:
            outfile.write(line)

    print("Processed", input_file)

if __name__ == "__main__":
    # Check if the correct number of arguments are provided
    if len(sys.argv) != 7:
        print("Usage: python script.py input_file output_file threshold")
        sys.exit(1)

    # Extract input and output file paths from command-line arguments
    input_file = sys.argv[2]
    output_file = sys.argv[4]
    threshold = float(sys.argv[6])

    # Call the function to process the file
    process_file(input_file, output_file, threshold)

# import os
# import sys

# def process_file(input_file, output_file):
#     # Initialize variables to keep track of whether to keep lines or not
#     keep_lines = False
#     filtered_lines = []

#     # Open the input file and read all lines
#     with open(input_file, 'r') as infile:
#         lines = infile.readlines()

#     # Loop through each line in the input file
#     for line in lines:
#         # Split the line into columns
#         columns = line.strip().split('\t')
#         # Check if it's a D entry and the 7th column is 'Bacteria'
#         if columns[5] == "D" and "Bacteria" in columns[7]:
#             keep_lines = True

#         # Check if it's a D entry and the 7th column is not 'Bacteria'
#         if keep_lines and columns[5] == 'D' and "Bacteria" not in columns[7]:
#             keep_lines = False

#         # If keep_lines is True, add the line to the filtered lines list
#         if keep_lines:
#             filtered_lines.append(line)

#     # Write the filtered lines to the output file
#     with open(output_file, 'w') as outfile:
#         for line in filtered_lines:
#             outfile.write(line)

#     print("Processed", input_file)

# if __name__ == "__main__":
#     # Check if the correct number of arguments are provided
#     if len(sys.argv) != 5:
#         print("Usage: python script.py input_file output_file")
#         sys.exit(1)

#     # Extract input and output file paths from command-line arguments
#     input_file = sys.argv[2]
#     output_file = sys.argv[4]

#     # Call the function to process the file
#     process_file(input_file, output_file)