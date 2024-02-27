// ----------------------------------------------
//             Kneaddata
// ----------------------------------------------

process kneaddata {
    label 'maxi'


    errorStrategy 'retry'

    maxRetries 3

    publishDir path: "${params.output_dir}/kneaddata", mode: 'copy', overwrite: true, pattern: "${sample}_kneaddata/*"
    publishDir path: "${params.output_dir}/kneaddata", mode: 'copy', overwrite: true, pattern: "${sample}_kneaddata/fastqc/*"

    input:
        //récupération du nom du fichier pour chaque paire et du chemin du fichier fastq
        //val: accède à la valeur d'entrée par son nom dans le script de processus.
        //path : Gère la valeur d'entrée comme un chemin, en plaçant correctement le fichier dans le contexte d'exécution.
        tuple val(sample), path(reads)
        val(ref_transcriptome_human)
        val(ref_genome_human)

    output:

        //tuple val(sample), file("${sample}_kneaddata/${sample}_Read_SILVA_128_LSUParc_SSUParc_ribosomal_RNA_bowtie2_paired_contam_1_{1,2}.fastq.gz")
        //file: Gère la valeur d'entrée en tant que fichier, en la mettant correctement en scène dans le contexte d'exécution.
        tuple val(sample), path("${sample}_kneaddata/${sample}_paired_{1,2}.fastq"), emit: paired
        //tuple val(sample), file("${sample}_kneaddata/${sample}_paired_SILVA_128_LSUParc_SSUParc_ribosomal_RNA_bowtie2_contam_{1,2}.fastq"), emit: contam
        //tuple val(sample), file("${params.output_dir}/knead_count/${sample}_read_count_table.tsv")
        path("${sample}_kneaddata/${sample}_read_count_table.tsv"), emit: tabCount
        path("${sample}_kneaddata/fastqc/${sample}{1,2}_fastqc.html")
        path("${sample}_kneaddata/fastqc/${sample}_paired_{1,2}_fastqc.html")
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
        #--store-temp-output 

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

    maxRetries 3

    publishDir "${params.output_dir}/merge_knead", mode: 'copy', overwrite: true, pattern: "merge_knead.tsv"

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

    maxRetries 3

    // Directive publishDir spécifique pour les fichiers .kreport
    publishDir path: "${params.output_dir}/kraken/${sample}_kraken", mode: 'copy', overwrite: true, pattern: "${sample}_reads_minimizer.kreport"

    // Directive publishDir spécifique pour les autres fichiers
    publishDir path: "${params.output_dir}/kraken/${sample}_kraken", mode: 'copy', overwrite: true, pattern: "${sample}_reads.krak"
    publishDir path: "${params.output_dir}/kraken/${sample}_kraken", mode: 'copy', overwrite: true, pattern: "${sample}_classified_1.fastq"
    publishDir path: "${params.output_dir}/kraken/${sample}_kraken", mode: 'copy', overwrite: true, pattern: "${sample}_classified_2.fastq"
    publishDir path: "${params.output_dir}/kraken/${sample}_kraken", mode: 'copy', overwrite: true, pattern: "${sample}_unclassified_1.fastq"
    publishDir path: "${params.output_dir}/kraken/${sample}_kraken", mode: 'copy', overwrite: true, pattern: "${sample}_unclassified_2.fastq"
    //publishDir path: "${params.output_dir}/kraken/${sample}_kraken", mode: 'copy', overwrite: true, pattern: ".command.log"

    input:
        tuple val(sample), file(reads)

    output:
        val(sample), emit: sampleName
        path("${sample}_reads.krak"), emit: krak
        path("${sample}_reads_minimizer.kreport"), emit: minimizer_report
        file("${sample}_classified_1.fastq")
        file("${sample}_classified_2.fastq")
        file("${sample}_unclassified_1.fastq")
        file("${sample}_unclassified_2.fastq")
        file(".command.log")
        

    script:
    """
        kraken2 --db ${params.kraken_db} \
        --threads ${task.cpus} \
        --paired \
        --classified-out ${sample}_classified#.fastq \
        --unclassified-out ${sample}_unclassified#.fastq \
        --output ${sample}_reads.krak \
        --report ${sample}_reads_minimizer.kreport \
        --report-minimizer-data \
        ${reads}
    """ 
} // No confidence score is defined for kraken!     

// ----------------------------------------------
//             Generate_original_kreport
// ----------------------------------------------
process generate_original_kreport {

    label 'midi'

    module 'python/3.7'

    errorStrategy 'retry'

    maxRetries 3


    publishDir "${params.output_dir}/GenerateOriginalKreport", mode: 'copy', overwrite: true

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

    maxRetries 3
    
    
    publishDir "${params.output_dir}/FilterKreportFile", mode: 'copy', overwrite: true
    
    scratch '/tmp'

    input:
        val(sample)
        file(kreport)

    output:
        val(sample)
        path("${sample}_filtered.kreport"), emit : kreport_filtered

    script:
    """
        python ${params.project_dir}/script/Minimizer_report_filtering.py -i ${kreport} -o ${sample}_filtered.kreport

    """
    }


// ----------------------------------------------
//             K-mer correlation test
// ----------------------------------------------

process kct {
    label 'midi'

    module 'r/4.3.1'

    publishDir "${params.output_dir}/kct_plot", mode: 'copy', overwrite: true

    scratch '/tmp'

    input:
        file(kreport_filtered)

    output:
        path("scatterplot.png")
        path("correlation.csv"), emit: correlation

    script:
    """
        Rscript ${params.project_dir}/script/kct.R "${params.output_dir}/FilterKreportFile"
    """ 

}

// ----------------------------------------------
//             kCT filtering
// ----------------------------------------------

process kct_filtering {
    label 'midi'
    
    module 'python/3.7'

    errorStrategy 'retry'

    maxRetries 3

    publishDir "${params.output_dir}/kCT_filtering", mode: 'copy', overwrite: true

    scratch'/tmp'

    input:
        val(sample)
        file(correlation)
        file(kreport)

    output:
        val(sample)
        path("${sample}_filtered.kreport"), emit : kreport_filtered
        path("taxa_to_remove.csv")
    script:
        """
        python ${params.project_dir}/script/kct_filter.py ${correlation} "${params.output_dir}/GenerateOriginalKreport" "${sample}_filtered.kreport" 
        """
}

// ----------------------------------------------
//             Filter_kreport_original_file
// ----------------------------------------------

process filter_kreport_original_file {
    label 'midi'

    module 'python/3.7'

    errorStrategy 'retry'

    maxRetries 3
    
    
    publishDir "${params.output_dir}/FilterKreportOriginalFile", mode: 'copy', overwrite: true
    
    scratch '/tmp'

    input:
        val(sample)
        file(kreport_original)

    output:
        val(sample)
        path("${sample}_original_filtered.kreport"), emit : kreport_original_filtered

    script:
    """
        python ${params.project_dir}/script/Original_report_filtering.py -i ${kreport_original} -o ${sample}_original_filtered.kreport

    """
    }


// ----------------------------------------------
//             Kraken2Mpa
// ----------------------------------------------

process run_Kraken2Mpa {
    label 'midi'
    label 'krakentools'

    module 'krakentools'

    errorStrategy 'retry'

    maxRetries 3


    publishDir "${params.output_dir}/kraken2Mpa", mode: 'copy', overwrite: true

    scratch '/tmp'

    input:
        val(sample)
        file(kreport_original_filtered) 

    output:
        path("${sample}_mpa.tsv"), emit: mpa_report

    script:
    """
        kreport2mpa.py -r ${kreport_original_filtered} -o ${sample}_mpa.tsv --display-header --no-intermediate-ranks --read_count
    """ 
}

// ----------------------------------------------
//             KronaReport
// ----------------------------------------------
    
process run_KronaReport {
    label 'midi'
    label 'krona'

    module 'krona'

    errorStrategy 'retry'

    maxRetries 3


    publishDir "${params.output_dir}/krona", mode: 'copy', overwrite: true

    scratch '/tmp'
        
    input: 
        val(sample)
        file(kreport_original_filtered) 

    output:
        tuple val(sample), file("${sample}_krona.html")

    script:
    """
        ktImportTaxonomy -m 3 -t 5 ${kreport_original_filtered} -o ${sample}_krona.html -tax ${params.taxonomy}
    """    
}


// ----------------------------------------------
//             Kraken2biom
// ----------------------------------------------

process run_Kraken2biom {
    label 'midi'
    
    label 'biom'

    errorStrategy 'retry'

    maxRetries 3


    publishDir "${params.output_dir}/kraken2biom", mode: 'copy', overwrite: true

    scratch '/tmp'

    input:
        file(inputFiles)

    output:
        path("all_sample_contig_table.biom"), emit: biom_data

    script:
    """
        echo ${inputFiles}
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

    maxRetries 3


    publishDir "${params.output_dir}/biom_tab", mode: 'copy', overwrite: true

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
//         MPA combine Contig
// -----------------------------------------------
process run_MpaCombineContig {
    label 'midi'
    label 'krakentools'

    module 'krakentools'

    errorStrategy 'retry'

    maxRetries 3


    publishDir "${params.output_dir}/mpaCombineContig", mode: 'copy', overwrite: true

    scratch '/tmp'

    input:
        file(mpa_report)

    output:
        file("combine_mpa.tsv")

    script:
    """
        combine_mpa.py -i ${(mpa_report).join(' ')} -o combine_mpa.tsv
    """ 
}

// -----------------------------------------------
//         Analyse krak
// -----------------------------------------------
/*
process analyse_krak {
    label 'midi'

    module 'gcc/11.2.0'

    errorStrategy 'retry'

    maxRetries 3


    publishDir "${params.output_dir}/analyse_krak/", mode: 'copy', overwrite: true

    scratch '/tmp'

    input:
        val(sample)
        file(krak)

    output:
        file("${sample}_tabCount.tsv")

    script:
    """
        gcc ${params.project_dir}/script/analyse_krak.c -o analyse_krak
        ./analyse_krak ${krak} ${params.taxonomy}/categories.dmp ${sample} 
    """ 
}
*/


// -----------------------------------------------
//         Keep_bacteria_only
// -----------------------------------------------

process keep_bacteria_only {
    label 'midi'

    module 'r/4.2.3'

    errorStrategy 'retry'

    maxRetries 3

    publishDir "${params.output_dir}/biom_tab", mode: 'copy', overwrite: true

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

ext = "fastq,fastq.gz,fastq.bz2,fq,fq.gz,fq.bz2,_001.fastq.gz"

Channel
    .fromFilePairs("${params.fastq_dir}/*{1,2}.{${ext}}")
    .ifEmpty { exit 1, "params.fastq_dir was empty - no input files supplied" }
    .set { reads_list } //stocke les pairs de fastq dans reads_list

workflow {
    kneaddata(reads_list, params.human_transcriptome, params.human_genome)
    merge_knead(kneaddata.out.tabCount.collect())
    kraken(kneaddata.out.paired)
    generate_original_kreport(kraken.out.sampleName, kraken.out.minimizer_report)
    filter_kreport_file(kraken.out.sampleName, kraken.out.minimizer_report)
    kct(filter_kreport_file.out.kreport_filtered.collect())
    kct_filtering(kraken.out.sampleName, kct.out.correlation, generate_original_kreport.out.kreport)
    filter_kreport_original_file(kraken.out.sampleName, generate_original_kreport.out.kreport)
    run_Kraken2Mpa(kraken.out.sampleName, generate_original_kreport.out.kreport)
    run_KronaReport(kraken.out.sampleName, filter_kreport_original_file.out.kreport_original_filtered)
    run_MpaCombineContig(run_Kraken2Mpa.out.mpa_report.collect())
    run_Kraken2biom(kct_filtering.out.kreport_filtered.collect())
    biom_tab(run_Kraken2biom.out.biom_data)
    keep_bacteria_only(biom_tab.out.tax_table, biom_tab.out.otu_table)
    //analyse_krak(kraken.out.sampleName, kraken.out.krak)
}