## Identification of bacteria in FFPE samples using RNAseq

## Description

This pipeline was designed to detect and identify non human DNA in RNA seq data, in order to have a separation between human and non human reads.  

It uses **Kneaddata** software for quality control and decontamination. It allows us to keep paired reads that do not match (or align) the databases used as filters.  

Then it uses **Kraken2** to assign each read that passed those filtering criteria to a known taxon.

It also uses **Nextflow** and **singularity** for easy deployment on any machine or HPC and to ensure reproducibility of results.


## Output files:

Once the pipeline is executed succesfully, you will obtain several results :

- Kraken step of the pipeline provides a taxonomic classification file for all the samples that were taken into account in the run.
- An OTU and TAX table will be generated at the end of the pipeline's execution. These files can be used then to perform several analysis, such as Bacterial diversity using Alpha and Beta diversity approaches, Bacterial enrichment set analysis, as well as the possibility to generate the  bacterial diversity pie-chart of the samples and to generate a table of abundance for all samples. Note that each folder of this repository contains all the necessary step to perform the several analysis mentionned earlier. 
