
## Required Input

## Fastq samples (paired-end)

For the pipeline to work correctly, input files needs to respect to following format: 

`*{1,2}.{ext}`

or:

`*{1,2}_001.{ext}`

where **_ext_** can be one of:  
`fastq`,`fastq.gz`, `fastq.bz2`, `fq`, `fq.gz`, `fq.bz2`


This format can be adjusted in the configuration file (test.config) by changing the following line:

`suffix = "*{1,2}"`

Into:

`suffix = "*{1,2}_001"`


#### Reference Genome

Download the reference genome of your choice (fasta and GTF). Here we download several reference genomes in format bt2 (bowtie2). 

Here, we make the choice to have:  
* Human genome database (3.5 Gb)
* Human transcriptome database (254 Mb)
* silva16s database (3 Gb)

To download these databases, run the following commands:

`pip install kneaddata`

`kneaddata_database --download human_genome bowtie2 $DIR`

`kneaddata_database --download human_transcriptome bowtie2 $DIR`

`kneaddata_database --download ribosomal_RNA bowtie2 $DIR`

NB: $DIR is the output directory where you desire to store the databases.


## Kraken database

The Kraken2 standard database can be downloaded [here](https://genome-idx.s3.amazonaws.com/kraken/k2_standard_20230314.tar.gz) (64GB).


## Configuration file Parameters

Please fill in the necessary parameters in the configuration file before starting the execution, according the need and the ressource of your project.

## The params: 

- **fastq_dir**: the directory of your input fastq files. If you want to select only a set of samples, provide the R1 and R2 files of these samples in this directory.
- **human_genome**: reference database of the human genome
- **human_transcriptome**: reference database of the human transcriptome
- **silva16s**: reference database silva16s, 16s RNA for the procaryote species. This parameter is not included by default, but if you want to remove ribosomal 16s RNAs from the analysis, include it by adding --reference-db ${ref_silva16s} in the Kneaddata step of the pipeline.
- **project_dir**: the directory that contains your project
- **output_dir**: the directory where you will obtain your results

