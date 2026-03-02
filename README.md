## Authors

- [Ali YOUNCHA](https://github.com/MrAli1582)

# Identification of bacteria in FFPE samples using RNAseq

## Description

This pipeline was designed to detect and identify non human RNA reads in RNA seq data, in order to have a separation between human and non human reads.

It uses **Kneaddata** software for quality control and decontamination. It allows us to keep paired reads that do not match (or align) the databases used as filters.  

Then it uses **Kraken2** to assign each read that passed those filtering criteria to a known taxon.

It also uses **Nextflow** and **singularity** for easy deployment on any machine or HPC and to ensure reproducibility of results.

**Warning**: This pipeline was developed on an HPC clusterr ([IFB](https://www.france-bioinformatique.fr/cluster-ifb-core/)), you might need a few adjustments to run it locally.

## Dependencies

To run this pipeline you will need to have installed on your machine:
* [Nextflow](https://www.nextflow.io/) (version 23.10.1)
* [Singularity](https://docs.sylabs.io/guides/3.5/user-guide/index.html) (version 3.5.3)
* [Kneaddata](https://github.com/biobakery/kneaddata)
* [Kraken2](https://github.com/DerrickWood/kraken2)

Please check the specific dependencies for each tool. For beginners, more information and command lines are available in this file: [README](https://github.com/GeNeHetX/MicrobioteRNApipe/blob/main/NextflowPipeline/Dependencies.md).

## Input files:

To run this pipeline you will need to provide:
* Fastq
* [Kneaddata database](https://github.com/biobakery/kneaddata?tab=readme-ov-file#download-the-database) with human genome, human transcriptome and silva16S directories
* [Kraken2 database](https://github.com/DerrickWood/kraken2/blob/master/docs/MANUAL.html)

Example command lines to download Kraken2 and Kneaddata databases are availables [here](https://github.com/GeNeHetX/MicrobioteRNApipe/blob/main/NextflowPipeline/Dependencies.md).

**Note**: If you are running this microbiome analysis in parallel with a human analysis, you can use the unmapped reads from your human alignment with STAR. This will reduce the size of the FASTQ files and decrease processing time.

## Arguments and Parameters:
To run, this pipeline requires a few arguments specific to yours files:
* fastq_dir : directory with the fastqs to process
* suffix : suffix of your fastq files
* human_genome : kneaddata database directory with human genome files
* human_transcriptome : kneaddata database directory  with human transcriptome files
* silva16s : kneaddata database directory  with silva16s files
* kraken_db : kraken2 database directory

* script_dir : path the 'PipelineScripts' directory of Microbiote_pdacrna pipeline
* output_dir : path to save output files
* work_dir : path to save the nextflow work directory

There is also two parameters you can adjust acording to your data and preferences:
* threshold : threshold for number of reads filtering
* confidence_score : confidence kraken score


All thoses parameters needs to be store in a config file. You can find an example [here](https://github.com/GeNeHetX/MicrobioteRNApipe/blob/main/NextflowPipeline/example.config).

## Pipeline execution steps
To run pipeline with nextflow:
```
nextflow -C path/to/file.config run path/to/MicrobioteRNApipe/NextflowPipeline/main_microbiote_pipeline.nf -resume
```

To run pipeline with a job file :
```
sbatch PATH_TO_JOB_EXECUTION PATH_TO_PIPELINE_MAIN PATH_TO_CONFIGURATION_FILE
```

## Output files:

Once the pipeline is executed successfully, you will obtain several results :

- Kraken kreports : Kraken step of the pipeline provides a taxonomic classification file for all the samples that were taken into account in the run.
- OTU and TAX tables : two tables will be generated at the end of the pipeline's execution. 

## Possible analysis from MicrobioteRNAPipe output results

- If you want to retrieve the microbial composition pie-chart of all your samples, check the following [README](https://github.com/GeNeHetX/MicrobioteRNApipe/blob/main/AnalysisScripts/2_Samples_Bacterial_PieCharts/README.md)

- if you didn't run all the samples on the same pipeline and want to generate one OTU and TAX table of all these samples, check this [README](https://github.com/GeNeHetX/MicrobioteRNApipe/blob/main/AnalysisScripts/1_MergeOTU_TAX_files/README.md) 

- If you want to perform Alpha and Beta diversity analyses using the OTU and TAX table outputs of the pipeline, check the first part of the [tutorial](https://github.com/GeNeHetX/MicrobioteRNApipe/blob/main/AlphaBetaAnalysis/Tutorial_Alpha_beta.md) and it's second part [tutorial_P2](https://github.com/GeNeHetX/MicrobioteRNApipe/blob/main/AlphaBetaAnalysis/Tutorial_Alpha_beta_part2.md)

- Finally, if you want to perform Bacterial set enrichment analysis, check the following [tutorial](https://github.com/GeNeHetX/MicrobioteRNApipe/blob/main/BacterialEnrichmentAnalysis/Tutorial_bacterial_enrichment_analysis.md)

