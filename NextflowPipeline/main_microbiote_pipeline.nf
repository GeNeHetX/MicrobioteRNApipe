// ----------------------------------------------
//             Kneaddata
// ----------------------------------------------

process kneaddata {
    label 'maxi'

    errorStrategy 'retry'

    maxRetries 1
    publishDir "${params.output_dir}/kraken_results", mode: 'copy'
    //publishDir path: "${params.output_dir}/kneaddata", mode: 'copy', overwrite: true, pattern: "${sample}_kneaddata/*"
    //publishDir path: "${params.output_dir}/kneaddata", mode: 'copy', overwrite: true, pattern: "${sample}_kneaddata/fastqc/*"

    input:
        //récupération du nom du fichier pour chaque paire et du chemin du fichier fastq
        //val: accède à la valeur d'entrée par son nom dans le script de processus.
        //path : Gère la valeur d'entrée comme un chemin, en plaçant correctement le fichier dans le contexte d'exécution.
        tuple val(sample_id), path(reads)
        val(ref_transcriptome_human)
        val(ref_genome_human)
        // val(ref_silva16s)

    output:

        //tuple val(sample), file("${sample}_kneaddata/${sample}_Read_SILVA_128_LSUParc_SSUParc_ribosomal_RNA_bowtie2_paired_contam_1_{1,2}.fastq.gz")
        //file: Gère la valeur d'entrée en tant que fichier, en la mettant correctement en scène dans le contexte d'exécution.
        tuple val(sample_id), path("${sample_id}_kneaddata/${sample_id}_paired_{1,2}.fastq"), emit: paired
        tuple val(sample_id), path("${sample_id}_kneaddata/${sample_id}_unmatched_{1,2}.fastq"), emit: unmatched
        //tuple val(sample), file("${sample}_kneaddata/${sample}_paired_SILVA_128_LSUParc_SSUParc_ribosomal_RNA_bowtie2_contam_{1,2}.fastq"), emit: contam
        //tuple val(sample), file("${params.output_dir}/knead_count/${sample}_read_count_table.tsv")
        path("${sample_id}_kneaddata/${sample_id}_read_count_table.tsv"), emit: tabCount

        // path("${sample}_kneaddata/fastqc/${sample}{1,2}_fastqc.html")
        // path("${sample}_kneaddata/fastqc/${sample}_paired_{1,2}_fastqc.html")
        file("${sample_id}_kneaddata/${sample_id}.log")
        
        
    script:

    def adapter_path = "/shared/projects/microbiote_pdacrna/AGASH/TruSeq3-PE-2-GGGGG.fa"
    def id = sample_id instanceof List ? sample_id[0] : sample_id

    """
    kneaddata \\
        --verbose \\
        -i1 ${reads[0]} \\
        -i2 ${reads[1]} \\
        -db /shared/projects/microbiote_pdacrna/AGASH/human_genome \\
        -db /shared/projects/microbiote_pdacrna/AGASH/human_transcriptome \\
        --remove-intermediate-output \\
        --bowtie2-options="--sensitive" \\
        --trimmomatic-options "ILLUMINACLIP:/shared/projects/microbiote_pdacrna/AGASH/TruSeq3-PE-2-GGGGG.fa:2:30:10 HEADCROP:3 SLIDINGWINDOW:4:20 MINLEN:50" \\
        --bypass-trf \\
        --run-fastqc-start \\
        --run-fastqc-end \\
        --threads 10 \\
        --output-prefix ${id} \\
        -o ${id}_kneaddata

    kneaddata_read_count_table \\
        --input ${id}_kneaddata \\
        --output ${id}_kneaddata/${id}_read_count_table.tsv
    """
}


// ----------------------------------------------
//             Trimmomatic 
// ----------------------------------------------

process trimmomatic {
    tag "Trimmomatic sur ${sample_id}"
    label 'maxi'

    publishDir "${params.output_dir}/trimmomatic_results", mode: 'copy'

    input:
    tuple val(sample_id), path(reads)

    output:
    tuple val(sample_id), path("${sample_id}_trim_paired_1.fastq"), path("${sample_id}_trim_paired_2.fastq"), emit: paired
    tuple val(sample_id), path("${sample_id}_trim_unpaired_*.fastq"), emit: unpaired

    script:
    """
    java -jar /usr/share/java/trimmomatic.jar PE -threads $task.cpus \\
        ${reads[0]} ${reads[1]} \\
        ${sample_id}_trim_paired_1.fastq ${sample_id}_trim_unpaired_1.fastq \\
        ${sample_id}_trim_paired_2.fastq ${sample_id}_trim_unpaired_2.fastq \\
        ILLUMINACLIP:${params.adapter_path}:2:30:10 \\
        HEADCROP:3 \\
        SLIDINGWINDOW:4:20 \\
        MINLEN:50
    """
}
// ----------------------------------------------
//             Merge Knead
// ----------------------------------------------

process star_filter {
    label 'maxi'

    publishDir "${params.output_dir}/reads_non_humains", mode: 'copy'

    input:
        tuple val(sample), path(reads)
        path star_index 

    output:
        tuple val(sample), path("${sample}_unmapped_{1,2}.fastq.gz"), emit: clean_reads

    script:
    """
    STAR --runThreadN ${task.cpus} \
         --genomeDir ${star_index} \
         --readFilesIn ${reads[0]} ${reads[1]} \
         --outFileNamePrefix ${sample}_ \
         --outReadsUnmapped Fastx \
         --outSAMtype None

    mv ${sample}_Unmapped.out.mate1 ${sample}_unmapped_1.fastq
    mv ${sample}_Unmapped.out.mate2 ${sample}_unmapped_2.fastq

    gzip ${sample}_unmapped_1.fastq
    gzip ${sample}_unmapped_2.fastq
    """
}



// ----------------------------------------------
//             Kraken
// ----------------------------------------------
process kraken {
    label 'maxi2'
    publishDir { "${params.output_dir}/kraken_calibration/score_${params.confidence_score}" }, mode: 'copy'

    input:
    tuple val(sample), path(paired_reads), path(unpaired_reads)

    output:
    tuple val(sample), path("${sample}_${params.confidence_score}.kreport"), emit: report
    path "${sample}_${params.confidence_score}.krak", emit: minimizer_report

    script:
    """
    for f in ${paired_reads} ${unpaired_reads}; do
        if [[ "\$f" == *.gz ]]; then
            zcat "\$f" >> all_reads_for_kraken.fastq
        else
            cat "\$f" >> all_reads_for_kraken.fastq
        fi
    done

    if grep -q "^@" all_reads_for_kraken.fastq; then
        kraken2 --db ${params.kraken_db} \
                --threads ${task.cpus} \
                --confidence ${params.confidence_score} \
                --report-minimizer-data \
                --output ${sample}_${params.confidence_score}.krak \
                --report ${sample}_${params.confidence_score}.kreport \
                all_reads_for_kraken.fastq 
    else
        touch ${sample}_${params.confidence_score}.krak
        echo -e "100.00\t0\t0\tU\t0\tunclassified" > ${sample}_${params.confidence_score}.kreport
    fi
    """
}

// ----------------------------------------------
//             Generate_original_kreport
// --------------------------------------------kraken_--
process generate_original_kreport {
    label 'midi'
    module 'python/3.7'
    
    publishDir "${params.output_dir}/kreports_standards", mode: 'copy'

    input:
        val(sample)
        path(kreport) 

    output:
        val(sample), emit: sampleName 
        path("${sample}_reads.kreport"), emit: kreport
    
    script:
    """
    python ${params.script_dir}/create_krakenfile.py ${kreport} ${sample}_reads.kreport
    """
}



// ----------------------------------------------
//             Filter_kreport_file
// ----------------------------------------------

process filter_kreport_file {
    label 'midi'

    module 'python/3.7'

    errorStrategy 'retry'

    maxRetries 1
    
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
        python ${params.script_dir}/Minimizer_report_filtering.py -i ${kreport} -o ${sample}_filtered.kreport -t ${params.threshold}

    """
    }


// ----------------------------------------------
//             Filter_kreport_original_file
// ----------------------------------------------

process filter_kreport_original_file {
    label 'midi'

    module 'python/3.7'

    errorStrategy 'retry'

    maxRetries 1
    
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
        python ${params.script_dir}/Original_report_filtering.py -i ${kreport_original} -o ${sample}_original_filtered.kreport -t ${params.threshold}

    """
    }




// ----------------------------------------------
//             Kraken2biom
// ----------------------------------------------

process run_Kraken2biom {
    label 'midi'
    
    label 'biom'

    errorStrategy 'retry'

    maxRetries 1


    publishDir "${params.output_dir}/kraken2biom", mode: 'copy', overwrite: true

    scratch '/tmp'

    input:
        file(inputFiles)

    output:
        path("all_sample_contig_table.biom"), emit: biom_data

    script:
    """
        kraken-biom ${inputFiles} --max S -o ./all_sample_contig_table.biom
    """ 
}



// ----------------------------------------------
//             Biom_tab
// ----------------------------------------------

process biom_tab {

    label 'midi'

    module 'r/4.2.3'

    errorStrategy 'retry'

    maxRetries 1


    // publishDir "${params.output_dir}", mode: 'copy', overwrite: true

    input:
        file(biom_data)

    output:
        path("tax_table.csv"), emit: tax_table
        path("otu_table.csv"), emit: otu_table

    script:
    """
        Rscript ${params.script_dir}/biom.r -tab ${biom_data} 
    """ 
}


// -----------------------------------------------
//         Keep_bacteria_only
// -----------------------------------------------

process keep_bacteria_only {
    label 'midi'

    module 'r/4.2.3'

    errorStrategy 'retry'

    maxRetries 1

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
        Rscript ${params.script_dir}/keep_bacteria.r ${tax_table} ${otu_table}
    """ 
}



// -----------------------------------------------
//         MAIN
// -----------------------------------------------




workflow {
    reads_ch = Channel.fromPath("${params.krakendir}/*.fastq.gz", checkIfExists: true)

    reads_list = reads_ch
   .map { file ->
       def sample_id = file.simpleName.replaceAll(/_R[12]$/, "")
       return tuple(sample_id, file)
   }
   .groupTuple(size: 2)

    reads_list.view { sample_id, files -> "Échantillon détecté : $sample_id | Fichiers: $files" }

    trimmomatic(reads_list)

    ch_paired_ready = trimmomatic.out.paired
        .map { it -> tuple(it[0], [it[1], it[2]]) }

    ch_kraken_joined = ch_paired_ready.join(trimmomatic.out.unpaired)

    ch_kraken_ready = ch_kraken_joined.map { it ->
        return tuple(it[0], it[1], it[2])
    }


    kraken(ch_kraken_ready)
    run_Kraken2biom(kraken.out.report.collect())
    //biom_tab(run_Kraken2biom.out.biom_data)
    //keep_bacteria_only(biom_tab.out.tax_table, biom_tab.out.otu_table)
}