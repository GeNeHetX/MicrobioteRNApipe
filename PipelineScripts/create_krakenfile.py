import sys
import shutil

## The goal of this function is to generate the kreport original output file, which will be used to generate abundance table for each sample
def remove_columns(input_file, output_file):
    with open(input_file, 'r') as f_input:
        with open(output_file, 'w') as f_output:
            for line in f_input:
                ## Remove columns 4 and 5, keeping the first three columns unchanged
                columns = line.strip().split('\t')
                new_line = '\t'.join(columns[:3] + columns[5:]) + '\n'  # Keeping the first three columns
                f_output.write(new_line)

if __name__ == "__main__":
    if len(sys.argv) != 3:
        print("Usage: python script.py <input_file> <output_file>")
        sys.exit(1)

    input_file = sys.argv[1]
    output_file = sys.argv[2]

    try:
        remove_columns(input_file, output_file)
        print("Les colonnes 4 et 5 ont été supprimées avec succès dans le fichier de sortie.")
    except FileNotFoundError:
        print("Le fichier d'entrée spécifié n'existe pas.")
    except Exception as e:
        print("Une erreur s'est produite:", e)
