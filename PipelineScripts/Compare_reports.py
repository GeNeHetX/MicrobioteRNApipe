import argparse

def read_file(file_path):
    data = []
    with open(file_path, 'r') as file:
        for line in file:
            line = line.strip()
            if line:
                parts = line.split('\t')
                data.append(parts)
    return data

def compare_files(file1_data, file2_data):
    missing_rows = []
    for row in file1_data:
        if row not in file2_data:
            missing_rows.append(row)
    return missing_rows

if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="Compare two files and find missing rows in the second file compared to the first one")
    parser.add_argument("file1", help="Path to the first file")
    parser.add_argument("file2", help="Path to the second file")
    args = parser.parse_args()

    file1_data = read_file(args.file1)
    file2_data = read_file(args.file2)

    missing_rows = compare_files(file1_data, file2_data)

    print("Missing rows in", args.file2, "compared to", args.file1, ":")
    for row in missing_rows:
        print("\t".join(row))
