# Identification of bacteria in FFPE samples using RNAseq

## Description
This pipeline was design to detect and identify non human DNA in RNA seq data.  
It uses **Kneaddata** software for data filtering and cleaning. It allows us to keep paired reads that do match (or align) the databases used as filters.  
Then **MetaSPAdes** is a _de novo_ assembly method.  
The pipeline then finally uses **Kraken2** to taxonomically assign them.  
Visualisation of the results are made with krona, or output as heatmap at different taxonomic level?
It uses **Nextflow** and **singularity** for easy deployment on any machine or HPC and to ensure reproducibility of results.


## Summary
The simplified workflow is described in the graph bellow:

## Version
Nextflow is written using DSL2

## Description of the process
- **kneaddata**: data filtering & data cleaning
- **merge_knead**: merge count tables that result from the Kneaddata process (uses Bowtie2 after Trimmomatic)
- **assembly_metaSPAdes**: assembly
- **kraken2**: taxonomic assignment
- **run_Kraken2Mpa**: convert the kreport file to a mpa file (visualisation of the results)
- **run_KronaReport**: convert the kreport file to a html that is a piechart representing our results
- **run_MpaCombineContig**: combine the mpa file 
- **run_Heatmap_reads**: create heatmaps
- **run_Kraken2biom**: see [kraken-biom](https://github.com/smdabdoub/kraken-biom)
- **biom_tab**: create two tables: OTU_table and TAX_table


## Prerequisite 

### > Install your environnement - computing cluster (IFB, etc) LINUX

Load the last version of nextflow:  
```module load nextflow/23.04.1```

add containers for '**kneaddata**' and '**run_Kraken2biom**' on "withName" attributes.

### > Install your environnement - local LINUX

To run this pipeline you will need to have installed on your machine:
* [Nextflow](https://www.nextflow.io/) (version 22.10.4)
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

### Required Input

#### Fastq samples (paired-end)

For the pipeline to work correctly, input files needs to
respect to following format: 

`*{1,2}.{ext}`   
where **_ext_** can be one of:  
`fastq`,`fastq.gz`, `fastq.bz2`, `fq`, `fq.gz`, `fq.bz2`

#### Reference Genome

Download the reference genome of your choice (fasta and GTF). Here we download several reference genomes in format bt2 (bowtie2). 

Here, we make the choice to have:  
* Human genome database
* Human transcriptome database
* silva16s database

#### Kraken database

The Kraken2 standard database can be downloaded [here](https://genome-idx.s3.amazonaws.com/kraken/k2_standard_20230314.tar.gz) (64GB).


#### Krona database


## Execution

```nextflow run main.nf -c config/<your_config.config>``` 

### Configuration file (test.config)

Please fill in the necessary parameters in the configuration file before starting the execution, according the need and the ressource of your project.

#### The params: 

- **fastq_dir**: the directory of your input fastq files
- **human_genome**: reference database of the human genome
- **human_transcriptome**: reference database of the human transcriptome
- **silva16s**: reference database silva16s, 16s RNA for the procaryote species
- **krona_tax**: taxonomy for krona
- **project_dir**: the directory that containes your entirely project
- **output_dir**: the directory where your results will be found


#### The process:

complete with the necessary ressources

Example:
- cpus: 10 
- memory: 165.GB


#### TODO:

# MicrobioteRNApipe