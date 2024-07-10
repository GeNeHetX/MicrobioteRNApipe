## Merging OTU and TAX table outputs for several runs of the pipeline.

If during the execution of the pipeline, you executed samples separately and want to assemble the different outputs of the OTU and TAX tables into:

- Download both merge_OTU_table.py and merge_TAX_table.py scripts that you can find in this folder.

- Then, you have to download all the OTU and TAX table files of each of your runs. Make sure that you put all the OTU tables on a folder, all the TAX tables on a another folder so that the script works correctly. Also, make sure to adjust the name of the different OTU tables, as well as TAX tables to avoid overwriting the files.

- Run these commands on your terminal:

`python merge_OTU_table.py PATH_TO_YOUR_OTU_TABLE_FOLDER PATH_TO_YOUR_OUTPUT_FILE`

`python merge_TAX_table.py PATH_TO_YOUR_TAX_TABLE_FOLDER PATH_TO_YOUR_OUTPUT_FILE`

With :

**PATH_TO_YOUR_OTU_TABLE_FOLDER** : The path containing the OTU tables of all the samples that were analyzed separately.

**PATH_TO_YOUR_OTU_TABLE_FOLDER** : The path containing the TAX tables of all the samples that were analyzed separately.

**PATH_TO_YOUR_OUTPUT_FILE** : The path where you want to generate the merged output.

Make sure that the scripts are present in the directory where you execute these commands.

