## Identification of bacteria in FFPE samples using RNAseq

## Description

This pipeline was designed to detect and identify non human DNA in RNA seq data, in order to have a separation between human and non human reads.  

It uses **Kneaddata** software for quality control and decontamination. It allows us to keep paired reads that do not match (or align) the databases used as filters.  

Then it uses **Kraken2** to assign each read that passed those filtering criterias to a known taxon.

It also uses **Nextflow** and **singularity** for easy deployment on any machine or HPC and to ensure reproducibility of results.


## Summary
The simplified workflow is described in the graph below:

## Version
Nextflow is written using DSL2

## Description of the process
- **kneaddata**: Quality control and decontamination tool
- **merge_knead**: Merge count tables that result from the Kneaddata process 
- **kraken2**: Taxonomic classification of reads that pass the quality control and decontamination step. This tool generates as an output a KREPORT file that contains several columns : Percentage of reads covered by the clade rooted at this taxon, Number of reads covered by the clade rooted at this taxon, Number of reads assigned directly to this taxon, number of total k-mers used, number of unique k-mers, A rank code indicating (U)nclassified, (R)oot, (D)omain, (K)ingdom, (P)hylum, (C)lass, (O)rder, (F)amily, (G)enus, or (S)pecies, NCBI taxonomic ID number and Indented scientific name.
- **Generate_original_kreport**: Allows to generate a kreport file that does not contain the columns of the number of unique k-mers and total used for a taxon.
- **Filter_kreport_file**: Allows to filter the output file of Kraken to keep only Bacterial reads.
- **Filter_kreport_original_file**: Allows to filter the Kreport file without the columns number of unique k-mers and total used. This file will then be used to perform the rest of the steps.
- **run_Kraken2biom**: Creates BIOM-Format tables from Kraken output, see [kraken-biom](https://github.com/smdabdoub/kraken-biom)
- **biom_tab**: Creates two tables: OTU_table and TAX_table


## Prerequisite 

### > Install your environnement - computing cluster (IFB, etc) LINUX

Load the last version of nextflow:  
```module load nextflow/23.10.1```

add containers for '**kneaddata**' and '**run_Kraken2biom**' on "withName" attributes.

### > Install your environnement - local LINUX

To run this pipeline you will need to have installed on your machine:
* [Nextflow](https://www.nextflow.io/) (version 23.10.1)
* [Singularity](https://docs.sylabs.io/guides/3.5/user-guide/index.html) (version 3.5.3)

#### Singularity:

install system dependencies:
```bash
  sudo apt-get update && sudo apt-get install -y \
  build-essential \
  libssl-dev \
  uuid-dev \
  libgpgme11-dev \
  squashfs-tools \
  libseccomp-dev \
  wget \
  pkg-config \
  git \
  cryptsetup
```

**install go**:
```bash
  install the last version of Go [Go](https://go.dev/dl/)
  # move the dezipped file:
  sudo mv go/ /usr/local/bin/
  sudo nano /etc/profile
  #ajouter cette ligne à la fin du fichier:
  export PATH=$PATH:/usr/local/go/bin
  go version
```

**install singularity**:

download last version, for example:
> singularity-ce-4.0.3.tar.gz

```bash
tar -xzf singularity-ce-4.0.3.tar.gz 
cd singularity-ce-4.0.3/
./mconfig
cd builddir/
```

```bash
./mconfig 
cd builddir/
make
sudo make -C install
sudo mv ~/singularity-ce-4.0.1 /usr/local/bin
sudo /usr/local/bin
singularity --version
```
À retester

## Required Input

## Fastq samples (paired-end)

For the pipeline to work correctly, input files needs to respect to following format: 

`*{1,2}.{ext}`   
where **_ext_** can be one of:  
`fastq`,`fastq.gz`, `fastq.bz2`, `fq`, `fq.gz`, `fq.bz2`

Note: If you have FASTQ files that contain _001 in their name after R1 or R2, use the following format:

`*{1,2}_001.{ext}`
where **_ext_** can be one of:  
`fastq`,`fastq.gz`, `fastq.bz2`, `fq`, `fq.gz`, `fq.bz2`

This format can be adjusted in the Channel section of the MicrobioteRNApipe (main_microbiote_pipeline.nf)

#### Reference Genome

Download the reference genome of your choice (fasta and GTF). Here we download several reference genomes in format bt2 (bowtie2). 

Here, we make the choice to have:  
* Human genome database (3.5 Gb)
* Human transcriptome database (254 Mb)
* silva16s database (3 Gb)

To download these databases, install kneaddata in your machine and see the Download database section in the Kneaddata tool github 
[here](https://github.com/biobakery/kneaddata/blob/master/readme.md#installation)


## Kraken database

The Kraken2 standard database can be downloaded [here](https://genome-idx.s3.amazonaws.com/kraken/k2_standard_20230314.tar.gz) (64GB).



## Execution

To ewecute the pipeline, run the following command :

```sbatch launch_nf.job main_microbiote_pipeline.nf test.config``` 

Adjust this command based on the paths to access the different files.

## Configuration file (test.config)

Please fill in the necessary parameters in the configuration file before starting the execution, according the need and the ressource of your project.

## The params: 

- **fastq_dir**: the directory of your input fastq files. If you want to select only a set of samples, provide the R1 and R2 files of these samples in this directory.
- **human_genome**: reference database of the human genome
- **human_transcriptome**: reference database of the human transcriptome
- **silva16s**: reference database silva16s, 16s RNA for the procaryote species. This parameter is not included by default, but if you want to remove ribosomal 16s RNAs from the analysis, include it by adding --reference-db ${ref_silva16s} in the Kneaddata step of the pipeline.
- **project_dir**: the directory that containes your entirely project
- **output_dir**: the directory where you will obtain your results

## Output files:

Once the pipeline is executed succesfully, you will obtain several results :

- Kraken step of the pipeline provides a taxonomic classification file for all the samples that were taken into account in the run.
- An OTU and TAX table will be generated at the end of the pipeline's execution. These files can then be used to perform microbial diversity analysis, that is made possible thanks to an R script that will be detailed in a different section of the github.

