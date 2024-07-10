#  Generating Genus/Species abundance table files from the set of kreport files

From the set of Kreport files that you obtain after executing the pipeline, it is possible to generate an abundance table for Genus/Species taxa.

* Download the Create_abundance_file.py script available in the 3_Samples_Abundance_Table folder

* Then, open a terminal and run the following command

`python generate_pie_charts.py -i PATH_TO_YOUR_KREPORT_FILES_FOLDER -o PATH_TO_YOUR_OUTPUT_FILE -t "G"`


With:

-i (PATH_TO_YOUR_KREPORT_FILES_FOLDER): The path containing all your Kreport files

-o (PATH_TO_YOUR_OUTPUT_FILE): The path where you desire to obtain your csv file

-t ("G"): The taxa level that will be taken into account to generate the abundance table. Adjust it to "S" in case you want to retrieve the species abundance table file.

Make sure that your script is located in the directory where you run this command.

When you generate this file, you can use it to perform Bacterial set enrichment analysis. Check the following [Tutorial](https://github.com/GeNeHetX/MicrobioteRNApipe/blob/main/BacterialEnrichmentAnalysis/Tutorial_bacterial_enrichment_analysis.md)