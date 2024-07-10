// ----------------------------------------------
//             Kneaddata
// ----------------------------------------------

process kneaddata {
    label 'maxi'

    errorStrategy 'retry'

    maxRetries 5

    //publishDir path: "${params.output_dir}/kneaddata", mode: 'copy', overwrite: true, pattern: "${sample}_kneaddata/*"
    //publishDir path: "${params.output_dir}/kneaddata", mode: 'copy', overwrite: true, pattern: "${sample}_kneaddata/fastqc/*"

    input:
        //récupération du nom du fichier pour chaque paire et du chemin du fichier fastq
        //val: accède à la valeur d'entrée par son nom dans le script de processus.
        //path : Gère la valeur d'entrée comme un chemin, en plaçant correctement le fichier dans le contexte d'exécution.
        tuple val(sample), path(reads)
        val(ref_transcriptome_human)
        val(ref_genome_human)
        // val(ref_silva16s)

    output:

        //tuple val(sample), file("${sample}_kneaddata/${sample}_Read_SILVA_128_LSUParc_SSUParc_ribosomal_RNA_bowtie2_paired_contam_1_{1,2}.fastq.gz")
        //file: Gère la valeur d'entrée en tant que fichier, en la mettant correctement en scène dans le contexte d'exécution.
        tuple val(sample), path("${sample}_kneaddata/${sample}_paired_{1,2}.fastq"), emit: paired
        //tuple val(sample), file("${sample}_kneaddata/${sample}_paired_SILVA_128_LSUParc_SSUParc_ribosomal_RNA_bowtie2_contam_{1,2}.fastq"), emit: contam
        //tuple val(sample), file("${params.output_dir}/knead_count/${sample}_read_count_table.tsv")
        path("${sample}_kneaddata/${sample}_read_count_table.tsv"), emit: tabCount
        // path("${sample}_kneaddata/fastqc/${sample}{1,2}_fastqc.html")
        // path("${sample}_kneaddata/fastqc/${sample}_paired_{1,2}_fastqc.html")
        file("${sample}_kneaddata/${sample}.log")
        
        
    script:
    """
        kneaddata \
        --verbose \
        --input1 ${reads[0]} \
        --input2 ${reads[1]} \
        --reference-db ${ref_genome_human} \
        --reference-db ${ref_transcriptome_human} \
        --remove-intermediate-output \
        --bowtie2-options="--sensitive" \
        --decontaminate-pairs "lenient" \
        --bypass-trf \
        --run-fastqc-start \
        --run-fastqc-end \
        --threads ${task.cpus} \
        --output-prefix ${sample} \
        --output ${sample}_kneaddata
         
        kneaddata_read_count_table \
        --input ${sample}_kneaddata \
        --output ${sample}_kneaddata/${sample}_read_count_table.tsv
    """
}

// ----------------------------------------------
//             Merge Knead
// ----------------------------------------------


process merge_knead {
    label 'midi'

    module 'python/3.7'

    errorStrategy 'retry'

    maxRetries 5

    // publishDir "${params.output_dir}/merge_knead", mode: 'copy', overwrite: true, pattern: "merge_knead.tsv"

    scratch '/tmp'

    input:
        path(knead_count)

    output:
        file("merge_knead.tsv")

    script:
    """
        python ${params.project_dir}/script/merge_knead.py -i ${knead_count.join(' ')} -o .
    """ 
}


// ----------------------------------------------
//             Kraken
// ----------------------------------------------

process kraken {
    label 'maxi2'

    module 'kraken2'
   
    errorStrategy 'retry'

    maxRetries 5

    // Directive publishDir spécifique pour les fichiers .kreport
    publishDir path: "${params.output_dir}", mode: 'copy', overwrite: true, pattern: "${sample}_reads_minimizer.kreport"

    // Directive publishDir spécifique pour les autres fichiers
    // publishDir path: "${params.output_dir}/kraken/${sample}_kraken", mode: 'copy', overwrite: true, pattern: "${sample}_reads.krak"
    // publishDir path: "${params.output_dir}/kraken/${sample}_kraken", mode: 'copy', overwrite: true, pattern: "${sample}_classified_1.fastq"
    // publishDir path: "${params.output_dir}/kraken/${sample}_kraken", mode: 'copy', overwrite: true, pattern: "${sample}_classified_2.fastq"
    // publishDir path: "${params.output_dir}/kraken/${sample}_kraken", mode: 'copy', overwrite: true, pattern: "${sample}_unclassified_1.fastq"
    // publishDir path: "${params.output_dir}/kraken/${sample}_kraken", mode: 'copy', overwrite: true, pattern: "${sample}_unclassified_2.fastq"
    // publishDir path: "${params.output_dir}/kraken/${sample}_kraken", mode: 'copy', overwrite: true, pattern: ".command.log"

    input:
        tuple val(sample), file(reads)

    output:
        val(sample), emit: sampleName
        path("${sample}_reads.krak"), emit: krak
        path("${sample}_reads_minimizer.kreport"), emit: minimizer_report
        // file("${sample}_classified_1.fastq")
        // file("${sample}_classified_2.fastq")
        // file("${sample}_unclassified_1.fastq")
        // file("${sample}_unclassified_2.fastq")
        // file(".command.log")
        

    script:
    """
        kraken2 --db ${params.kraken_db} \
        --threads ${task.cpus} \
        --paired \
        --output ${sample}_reads.krak \
        --report ${sample}_reads_minimizer.kreport \
        --report-minimizer-data \
        --confidence 0.1 \
        ${reads} 
    """ 
}    



// ----------------------------------------------
//             Generate_original_kreport
// ----------------------------------------------
process generate_original_kreport {

    label 'midi'

    module 'python/3.7'

    errorStrategy 'retry'

    maxRetries 5


    // publishDir "${params.output_dir}/GenerateOriginalKreport", mode: 'copy', overwrite: true

    scratch '/tmp'

    input:
        val(sample)
        file(kreport)

    output:
        val(sample)
        path("${sample}_reads.kreport"), emit : kreport
    
    script:
    """
        python ${params.project_dir}/script/create_krakenfile.py ${kreport} ${sample}_reads.kreport
    """
}



// ----------------------------------------------
//             Filter_kreport_file
// ----------------------------------------------

process filter_kreport_file {
    label 'midi'

    module 'python/3.7'

    errorStrategy 'retry'

    maxRetries 5
    
   //  publishDir "${params.output_dir}/FilterKreportFile", mode: 'copy', overwrite: true
    
    scratch '/tmp'

    input:
        val(sample)
        file(kreport)

    output:
        val(sample)
        path("${sample}_filtered.kreport"), emit : kreport_filtered

    script:
    """
        python ${params.project_dir}/script/Minimizer_report_filtering.py -i ${kreport} -o ${sample}_filtered.kreport -t ${params.threshold}

    """
    }


// ----------------------------------------------
//             Filter_kreport_original_file
// ----------------------------------------------

process filter_kreport_original_file {
    label 'midi'

    module 'python/3.7'

    errorStrategy 'retry'

    maxRetries 5
    
    // publishDir "${params.output_dir}/FilterKreportOriginalFile", mode: 'copy', overwrite: true
    
    scratch '/tmp'

    input:
        val(sample)
        file(kreport_original)

    output:
        val(sample)
        path("${sample}_original_filtered.kreport"), emit : kreport_original_filtered

    script:
    """
        python ${params.project_dir}/script/Original_report_filtering.py -i ${kreport_original} -o ${sample}_original_filtered.kreport -t ${params.threshold}

    """
    }




// ----------------------------------------------
//             Kraken2biom
// ----------------------------------------------

process run_Kraken2biom {
    label 'midi'
    
    label 'biom'

    errorStrategy 'retry'

    maxRetries 5


    // publishDir "${params.output_dir}/kraken2biom", mode: 'copy', overwrite: true

    scratch '/tmp'

    input:
        file(inputFiles)

    output:
        path("all_sample_contig_table.biom"), emit: biom_data

    script:
    """
        kraken-biom ${inputFiles.join(' ')} --fmt json -o ./all_sample_contig_table.biom
    """ 
}



// ----------------------------------------------
//             Biom_tab
// ----------------------------------------------

process biom_tab {

    label 'midi'

    module 'r/4.2.3'

    errorStrategy 'retry'

    maxRetries 5


    publishDir "${params.output_dir}", mode: 'copy', overwrite: true

    input:
        file(biom_data)

    output:
        path("tax_table.csv"), emit: tax_table
        path("otu_table.csv"), emit: otu_table

    script:
    """
        Rscript ${params.project_dir}/script/biom.r -tab ${biom_data} 
    """ 
}


// -----------------------------------------------
//         Keep_bacteria_only
// -----------------------------------------------

process keep_bacteria_only {
    label 'midi'

    module 'r/4.2.3'

    errorStrategy 'retry'

    maxRetries 5

    publishDir "${params.output_dir}", mode: 'copy', overwrite: true

    scratch '/tmp'

    input:
        file(tax_table)
        file(otu_table)

    output:
        path("tax_table.csv"), emit: tax_table
        path("otu_table.csv"), emit: otu_table
    script:
    """
        Rscript ${params.project_dir}/script/keep_bacteria.r ${tax_table} ${otu_table}
    """ 
}



// -----------------------------------------------
//         MAIN
// -----------------------------------------------

// ext = "fastq,fastq.gz,fastq.bz2,fq,fq.gz,fq.bz2"

Channel
    .fromFilePairs("${params.fastq_dir}/${params.suffix}.{${params.ext}}")
    .ifEmpty { exit 1, "params.fastq_dir was empty - no input files supplied" }
    .set { reads_list } //stocke les pairs de fastq dans reads_list

workflow {
    kneaddata(reads_list, params.human_transcriptome, params.human_genome)
    merge_knead(kneaddata.out.tabCount.collect())
    kraken(kneaddata.out.paired)
    generate_original_kreport(kraken.out.sampleName, kraken.out.minimizer_report)
    filter_kreport_file(kraken.out.sampleName, kraken.out.minimizer_report)
    filter_kreport_original_file(kraken.out.sampleName, generate_original_kreport.out.kreport)
    run_Kraken2biom(generate_original_kreport.out.kreport.collect())
    biom_tab(run_Kraken2biom.out.biom_data)
    keep_bacteria_only(biom_tab.out.tax_table, biom_tab.out.otu_table)
}