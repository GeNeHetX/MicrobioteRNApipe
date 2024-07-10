## Generating bacterial composition pie chart

To generate the bacterial composition of the samples that were analyzed using the MicrobioteRNAPipe :

- Download the generate_pie_charts.py script that is available within this folder.

- Then, open a terminal and run the following command:

`python generate_pie_charts.py -i PATH_TO_YOUR_KREPORT_FILES_FOLDER -o PATH_TO_YOUR_OUTPUT_DIRECTORY -t "G" -n 20`

With:

**-i (PATH_TO_YOUR_KREPORT_FILES_FOLDER)**: The path containing all your Kreport files

**-o (PATH_TO_YOUR_OUTPUT_DIRECTORY)**: The path where you desire to obtain your pie charts and the abundance percentage CSV file.

**-t ("G")**: The taxa level that will be taken into account to generate the pie-chart. Adjust it to "S" in case you want to retrieve the species bacterial abundance pie chart.

**-n ("20")**: The number of top taxa that will be shown without the pie-chart. By default the value is 15, and it can be adjusted by changing the this parameter.

Make sure that your script is located in the directory where you run this command.

Also, for better visualization, use the CSV file to generate a pie-chart in case the python generated pie-chart are not clear (which can be the case)


