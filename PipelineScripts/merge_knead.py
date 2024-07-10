import os
import pandas as pd
import argparse

parser = argparse.ArgumentParser()

parser.add_argument("-i", "--input_paths", nargs='+', type=str, help="paths to the tsv files")
parser.add_argument("-o", "--output_dir", type=str, help="output directory for the merged file")
args = parser.parse_args()

# Check if output directory exists, if not, create it
if not os.path.exists(args.output_dir):
    os.makedirs(args.output_dir)

# Lire chaque fichier et stocker les DataFrames dans une liste
dfs = [pd.read_csv(fichier, sep='\t') for fichier in args.input_paths]

# Fusionner les DataFrames sur les colonnes
resultat = pd.concat(dfs, ignore_index=True)

# Écrire le résultat dans un fichier "merge_knead.tsv"
resultat.to_csv(os.path.join(args.output_dir, 'merge_knead.tsv'), sep='\t', index=False)

# Afficher le résultat final
print(resultat)
