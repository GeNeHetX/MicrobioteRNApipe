# 1. Git cloning and requirements for the pipeline execution

First of all, you have to clone this git repository using the following command:

`git clone https://github.com/GeNeHetX/MicrobioteRNApipe.git `

Then, in case that you don't have the following databases :

- Human genome database 
- Human transcriptome database 
- Silva16s database
- Kraken2 database

Download the first three databases by running the  the following commands in your terminal : 

```
pip install kneaddata

kneaddata_database --download human_genome bowtie2 HUMAN_GENOME_DATABASE_OUTPUT_FOLDER

kneaddata_database --download human_transcriptome bowtie2 HUMAN_TRANSCRIPTOME_DATABASE_OUTPUT_FOLDER

kneaddata_database --download ribosomal_RNA bowtie2 SILVA16S_DATABASE_OUTPUT_FOLDER
```

Then download the kraken database using the following [link](https://genome-idx.s3.amazonaws.com/kraken/k2_standard_20230314.tar.gz) (64Gb)

 
# 2. Transfering the required databases and files to the IFB cluster

In the IFB cluster, you have to create a directory where you will copy all the databases and the different files required to execute the pipeline.

To do so, use the following command :

```
cd /shared/projects/PROJECT_NAME/

mkdir DIRECTORY_NAME
```

Inside of this directory, create a folder for each of the databases using the following commands:

```
cd DIRECTORY_NAME

mkdir human_genome

mkdir human_transcriptome

mkdir silva16s

mkdir KrakenDB
```

Then, copy all the databases that you downloaded to the directory directory using the following commands:

```
rsync -v -r HUMAN_GENOME_DATABASE_OUTPUT_FOLDER/* YOUR_IFB_USERNAME@core.cluster.france-bioinformatique.fr:/shared/projects/PROJECT_NAME/DIRECTORY_NAME/human_genome
```
```
rsync -v -r HUMAN_TRANSCRIPTOME_DATABASE_OUTPUT_FOLDER/* YOUR_IFB_USERNAME@core.cluster.france-bioinformatique.fr:/shared/projects/PROJECT_NAME/DIRECTORY_NAME/human_transcriptome
```
```
rsync -v -r SILVA16S_DATABASE_OUTPUT_FOLDER/* YOUR_IFB_USERNAME@core.cluster.france-bioinformatique.fr:/shared/projects/PROJECT_NAME/DIRECTORY_NAME/silva16s
```
```
rsync -v -r KRAKEN_DATABASE_OUTPUT_FOLDER/* YOUR_IFB_USERNAME@core.cluster.france-bioinformatique.fr:/shared/projects/PROJECT_NAME/DIRECTORY_NAME/KrakenDB
```
`PROJECT_NAME` is the name of your project in the IFB cluster (for example microbiotePDACRNA)


Once you have copied all the databases, you have to create a directory where you want to copy all your FastQ files to the IFB cluster.

VERY IMPORTANT : Please copy your FASTQ files on a completely different project because they take a lot of space.

To do so, use the following command:

```
rsync -v -r PATH_TO_YOUR_FASTQ_FOLDER/* YOUR_IFB_USERNAME@core.cluster.france-bioinformatique.fr:/shared/projects/SAMPLES_PROJECT/DIRECTORY_NAME/
```

Afterwards, copy the necessary files to execute the pipeline `launch_nf.job`,  `main_microbiote_pipeline.nf` and `test.config`.

First, create a directory on the IFB cluster:

```
cd /shared/projects/PROJECT_NAME/DIRECTORY_NAME/

mkdir main
```

Then, go back to the terminal and run the following commands:

```
rsync -v -r MicrobioteRNAPipe/NextflowPipeline/launch_nf.job YOUR_IFB_USERNAME@core.cluster.france-bioinformatique.fr:/shared/projects/PROJECT_NAME/main/
```
```
rsync -v -r MicrobioteRNAPipe/NextflowPipeline/main_mircrobiote_pipeline.nf YOUR_IFB_USERNAME@core.cluster.france-bioinformatique.fr:/shared/projects/PROJECT_NAME/main/
```
```
rsync -v -r MicrobioteRNAPipe/NextflowPipeline/test.config YOUR_IFB_USERNAME@core.cluster.france-bioinformatique.fr:/shared/projects/PROJECT_NAME/main/
```
The last step is to create a script folder inside your project that will contain all the necessary steps to execute the pipeline:

```
cd ..
mkdir script
```
Once the folder is created, run this command on your terminal:

```
rsync -v -r MicrobioteRNAPipe/ YOUR_IFB_USERNAME@core.cluster.france-bioinformatique.fr:/shared/projects/PROJECT_NAME/script/
```

# 3. Parameters modification

Now that you have all the requirements to execute the pipeline, the test.config file should be adjusted to successfully execute the pipeline.

This file contains several parameters :

- **fastq_dir**: the directory of your input fastq files. If you want to select only a set of samples, provide the R1 and R2 files of these samples in this directory.
- **human_genome**: reference database of the human genome
- **human_transcriptome**: reference database of the human transcriptome
- **silva16s**: reference database silva16s, 16s RNA for the procaryote species. This parameter is not included by default, but if you want to remove ribosomal 16s RNAs from the analysis, include it by adding --reference-db ${ref_silva16s} in the Kneaddata step of the pipeline.
- **kraken_db**: The kraken2 database used for taxonomic classification
- **confidence_score**: The confidence score parameter value selected to perform taxonomic classification
- **threshold**: The number of reads assigned threshold that is used to filter genus/species taxa. 
- **project_dir**: the project directory where you have all the necessary files and databases to execute the pipeline
- **output_dir**: the directory where you want to obtain your results

##  3.1  Input files parameter

The FASTQ input file parameter can be adjusted on the `test.config` configuration file of the pipeline.

Here is an example of how this parameter is defined:

- **fastq_dir** = `"/shared/projects/pancreas_microbiote_rna/set_ech_melanome"`

NB: If you have samples that end with _001 (such as retina_A1_R1_001.fastq), you have also to replace this parameter:

```
suffix = "*{1,2}"
```
with:

```
suffix = "*{1,2}_001"
```
##  3.2  Database parameters

The different databases needed to run the pipeline are present on the `test.config` configuration file of the pipeline.

- **human_genome** = `"/shared/projects/microbiote_pdacrna/anne/newMetagen/ref/human_genome"`

- **human_transcriptome** = `"/shared/projects/microbiote_pdacrna/anne/newMetagen/ref/human_transcriptome"`

- **silva16s** = `"/shared/projects/microbiote_pdacrna/anne/newMetagen/ref/silva16S"`

- **kraken_db** = `"/shared/projects/vmdc_wgte/krakenDB"`

Adjust them as needed based on the localization of your databases in the IFB cluster

##  3.3   Filtering parameters

The filtering parameters that are required for certain tools of the pipeline can be adjusted on the `test.config` configuration file of the pipeline.

- **confidence_score**: `"0.1"`

- **threshold**: `"0"`


##  3.4   Project and Output directories

- **project_dir** = `"/shared/projects/can_mic/newMetagen"`

- **output_dir** = `"/shared/projects/microbiote_pdacrna/ali/Melanoma_run"`


These parameters can be adjusted as needed to match your actual paths.


# 4. Pipeline execution on the IFB cluster


Once u have adjusted all the parameters to match what you desire to analyze, executing the pipeline will require running three commands in the following order:


First, in your IFB terminal, move to the desired folder where you want to obtain the results of the execution of the pipeline using this command:

```
cd $output_dir
```

Then, create a screen to visualize the execution state of the pipeline using this command:

```
screen -S SCREEN_NAME
```

Finally, run the following command to execute the pipeline:

```
sbatch PATH_TO_JOB_EXECUTION PATH_TO_PIPELINE_MAIN PATH_TO_CONFIGURATION_FILE
```
With:

PATH_TO_JOB_EXECUTION: The path containing the launch_nf.job execution file 

`(e.g : /shared/projects/can_mic/newMetagen/script/launch_nf.job)`

PATH_TO_PIPELINE_MAIN: The path containing the main_microbiote_pipeline.nf pipeline execution file. 

`(e.g : /shared/projects/can_mic/newMetagen/main/main_microbiote_pipeline.nf)`

PATH_TO_CONFIGURATION_FILE: The path containing the .config pipeline configuration file. 

`(e.g : /shared/projects/can_mic/newMetagen/config/test.config)`



