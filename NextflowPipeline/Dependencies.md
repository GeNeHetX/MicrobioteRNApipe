# Dependencies 

## > Install your environnement - local LINUX

To run this pipeline you will need to have installed on your machine:
* [Nextflow](https://www.nextflow.io/) (version 23.10.1)
* [Singularity](https://docs.sylabs.io/guides/3.5/user-guide/index.html) (version 3.5.3)

## > Install your environnement - computing cluster (IFB, etc) LINUX

Load the last version of nextflow:  
```module load nextflow/23.10.1```

### Singularity:

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

download the latest version, for example:
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
## Git cloning and requirements for the pipeline execution

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

 