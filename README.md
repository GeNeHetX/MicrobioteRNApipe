## Authors

- [Ali YOUNCHA](https://github.com/MrAli1582)

# Identification of bacteria in FFPE samples using RNAseq

## Description

This pipeline was designed to detect and identify non human RNA reads in RNA seq data, in order to have a separation between human and non human reads.  

It uses **Kneaddata** software for quality control and decontamination. It allows us to keep paired reads that do not match (or align) the databases used as filters.  

Then it uses **Kraken2** to assign each read that passed those filtering criteria to a known taxon.

It also uses **Nextflow** and **singularity** for easy deployment on any machine or HPC and to ensure reproducibility of results.


## 00_Dependencies

To check the necessary dependencies before running the pipeline, click on the following [README](https://github.com/GeNeHetX/MicrobioteRNApipe/blob/main/NextflowPipeline/Dependencies.md) 


## 01_Pipeline execution steps

The execution of the pipeline requires following several steps that can be found in this [README](https://github.com/GeNeHetX/MicrobioteRNApipe/blob/main/NextflowPipeline/PipelineExecution.md)

## Output files:

Once the pipeline is executed successfully, you will obtain several results :

- Kraken step of the pipeline provides a taxonomic classification file for all the samples that were taken into account in the run.
- An OTU and TAX table will be generated at the end of the pipeline's execution. 

## Possible analysis from MicrobioteRNAPipe output results

- If you want to retrieve the microbial composition pie-chart of all your samples, check the following [README](https://github.com/GeNeHetX/MicrobioteRNApipe/blob/main/PythonScripts/2_Samples_Bacterial_PieCharts/README.md)

- if you didn't run all the samples on the same pipeline and want to generate one OTU and TAX table of all these samples, check this [README](https://github.com/GeNeHetX/MicrobioteRNApipe/blob/main/PythonScripts/1_MergeOTU_TAX_files/README.md) 

- If you want to perform Alpha and Beta diversity analyses using the OTU and TAX table outputs of the pipeline, check the first part of the [tutorial](https://github.com/GeNeHetX/MicrobioteRNApipe/blob/main/AlphaBetaAnalysis/Tutorial_Alpha_beta.pdf) and it's second part [tutorial_P2](https://github.com/GeNeHetX/MicrobioteRNApipe/blob/main/AlphaBetaAnalysis/Tutorial_Alpha_beta_part2.pdf)

- Finally, if you want to perform Bacterial set enrichment analysis, check the following [tutorial](https://github.com/GeNeHetX/MicrobioteRNApipe/blob/main/BacterialEnrichmentAnalysis/Tutorial_bacterial_enrichment_analysis.pdf)

