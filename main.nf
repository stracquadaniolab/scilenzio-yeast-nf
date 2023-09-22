// enabling nextflow DSL v2
nextflow.enable.dsl=2

// process BUILD_REFERENCE_TRANSCRIPTOME {

//     publishDir "${params.resultsDir}/combined-genome", mode: 'copy', overwrite: true

//     input:
//         path(wt_reference)
//         path(syn_reference)
//         path(wt_gff)
//         path(syn_gff)

//     output:
//         path("wt-syn-chr11-ref.fasta"), emit: reference
//         path("wt-syn-chr11-ref.gff"), emit: annotation

//     """
//         # Replace fasta header of synthetic reference to "chr111" to indicate synthetic contig
//         # concatenate wt and synthetic reference (in that order)
//         cat ${wt_reference} <(sed 's/^>chr11/>chr111/' ${syn_reference}) > wt-syn-chr11-ref.fasta

//         # Change chromosome and gene names in synthetic GFF
//         # Replace "chr11" with "chr111" to indicate synthetic chromosome
//         sed 's/chr11/chr111/g' ${syn_gff} > tmp-syn.gff
//         # Add "x." prefix to synthetic GFF "ID=" column to indicate synthetic genes
//         awk -F'\t' -v OFS='\t' '{ sub(/ID=/, "ID=x.", $9); sub(/Name=/, "Name=x.", $9); print }' tmp-syn.gff > modified-syn-chr11.gff

//         # concatenate wt and synthetic annotation
//         cat ${wt_gff} modified-syn-chr11.gff > wt-syn-chr11-ref.gff
//     """

//     stub:
//     """
//         touch wt-syn-chr11-ref.fasta
//         touch wt-syn-chr11-ref.gff
//     """

// }

process TRIM_READS {

    tag "${sample}"

    publishDir "${params.resultsDir}/qc/reads", pattern: "*.json", mode: 'copy', overwrite: true
    publishDir "${params.resultsDir}/qc/reads", pattern: "*.html", mode: 'copy', overwrite: true

    input:
        tuple val(sample), path(read1), path(read2)
    output:
        tuple val(sample), path("${read1.simpleName}.trimmed.fastq.gz"), path("${read2.simpleName}.trimmed.fastq.gz"), emit: fastq
        path "${sample}.fastp.json", emit: qc
    
    script:
    """
        fastp -w ${task.cpus} \\
            ${params.fastp.args}  \\
            --in1 ${read1} \\
            --in2 ${read2} \\
            --out1 ${read1.simpleName}.trimmed.fastq.gz \\
            --out2 ${read2.simpleName}.trimmed.fastq.gz \\
            --json ${sample}.fastp.json \\
            --html ${sample}.fastp.html
    """

    stub:
    """
        touch ${read1.simpleName}.trimmed.fastq.gz
        touch ${read2.simpleName}.trimmed.fastq.gz
        touch ${sample}.fastp.json  
    """

}

process GENERATE_GENOME_INDEX {

    tag "genome_index"

    input:
        path(genome_fasta)
        path(annotation_gtf)

    output:
        path("genome-index")

    script:
    """
        STAR \\
            --runThreadN ${task.cpus} \\
            --runMode genomeGenerate \\
            --genomeDir genome-index \\
            --genomeFastaFiles ${genome_fasta} \\
            --sjdbGTFfile ${annotation_gtf} \\
            --sjdbOverhang ${params.star.sjdbOverhang} \\
            --genomeSAindexNbases ${params.star.genomeSAindexNbases}
    """

    stub:
    """
        mkdir genome-index
    """

}

process ALIGN_READS {

    publishDir "${params.resultsDir}/bams/", mode: 'copy', overwrite: true

    tag "${sample}"

    input:
        path(genome_index)
        tuple val(sample), path(trimmed_read1), path(trimmed_read2)

    output:
        tuple val(sample), path("${sample}.Aligned.sortedByCoord.out.bam")

    script:
    """
        STAR \\
            --runThreadN ${task.cpus} \\
            --genomeDir ${genome_index} \\
            --readFilesIn ${trimmed_read1} ${trimmed_read2} \\
            --readFilesCommand zcat \\
            --outFileNamePrefix ${sample}. \\
            --outSAMtype BAM SortedByCoordinate \\
            --alignIntronMax ${params.star.alignIntronMax} \\
            --limitBAMsortRAM ${params.star.limitBAMsortRAM} \\
            --outBAMsortingBinsN ${params.star.outBAMsortingBinsN}
    """

    stub: 
    """
        touch ${sample}.Aligned.sortedByCoord.out.bam
    """
    
}

process SALMON_INDEX {

    input: 
        path(transcriptome)
        path(genome)

    output: 
        path("transcriptome-index")

    script:
    """
        # extract names of genome targets
        grep '^>' < ${genome} | cut -d " " -f 1 > decoys.txt
        sed -i.bak -e 's/>//g' decoys.txt

        # concatenate transcriptome and genome reference file for index
        cat ${transcriptome} ${genome} > gentrome.fa.gz

        salmon index \\
            --threads ${task.cpus} \\
            --transcripts gentrome.fa.gz \\
            --decoys decoys.txt \\
            --index transcriptome-index \\
            --keepDuplicates \\
            ${params.salmon.index.args}
    """

    stub:
    """
        mkdir transcriptome-index
    """

}

process SALMON_QUANT {

    tag "${sample}"

    publishDir "${params.resultsDir}/salmon-quant/", mode: 'copy', overwrite: true

    input: 
        path(transcriptome_index)
        tuple val(sample), path(read1), path(read2)

    output:
        path("${sample}")

    script:
    """
        salmon quant \\
            --threads ${task.cpus} \\
            --libType ${params.salmon.quant.libtype} \\
            --index ${transcriptome_index} \\
            --mates1 ${read1} \\
            --mates2 ${read2} \\
            --output ${sample} \\
            ${params.salmon.quant.args}
    """

    stub:
    """
    mkdir ${sample}
    """
}

process SUMMARIZE_TO_GENE {

    publishDir "${params.resultsDir}/dataset/", mode: 'copy', overwrite: true

    input: 
        path(sample_sheet)
        path(annotation_gtf)
        path(salmon_results, stageAs: 'salmon-quant/quant*.sf')

    output:
        path 'summarized-experiment.rds'
    
    script:
    """
        quick-rnaseq-summarize-to-gene.R \\
            ${sample_sheet} \\
            -c ${params.summarize_to_gene.counts_from_abundance} \\
            -d ${params.summarize_to_gene.organism_db}
    """

    stub:
    """
    touch summarized-experiment.rds
    """
}

// TODO
process ANALYSIS_DGE {
    
    tag "${contrast1}-vs-${contrast2}"

    publishDir "${params.resultsDir}/analysis/", mode: 'copy', overwrite: true

    input: 
        path sefile
        tuple val(contrast1), val(contrast2)

    output:
        tuple path("dge-${contrast1}-vs-${contrast2}.csv"), val(contrast1), val(contrast2)
    
    script:
    """
        quick-rnaseq-dge.R \\
            ${sefile} \\
            dge-${contrast1}-vs-${contrast2}.csv \\
            --case ${contrast1} \\
            --control ${contrast2} \\
            -l ${params.dge.lfc_threshold} \\
            -f ${params.dge.fdr}
    """

    stub:
    """
    touch dge-${contrast1}-vs-${contrast2}.csv
    """
}

workflow {

    // sample channels
    samplesheet_file = file(params.samplesheet)
    samples_ch = channel.from(samplesheet_file)
                    .splitCsv(header: true)
                    .map{ record -> tuple(record.sample, file(record.read1), file(record.read2)) }

    // // build reference
    // BUILD_REFERENCE_TRANSCRIPTOME(file(params.reference.wt), file(params.reference.syn), file(params.annotation.wt), file(params.annotation.syn))

    // preprocess reads
    TRIM_READS(samples_ch)

    // perform alignment
    GENERATE_GENOME_INDEX(BUILD_REFERENCE_TRANSCRIPTOME.out.reference, BUILD_REFERENCE_TRANSCRIPTOME.out.annotation)
    ALIGN_READS(GENERATE_GENOME_INDEX.out, PREPROCESS_READS.out.trimmed_fastq)

    // create transcriptome index
    tx_ref_file = file(params.transcriptome.reference)
    tx_decoy_file = file(params.transcriptome.decoys)
    SALMON_INDEX(tx_ref_file, tx_decoy_file)

    // estimate transcript level abundance
    SALMON_QUANT(SALMON_INDEX.out, TRIM_READS.out.fastq)

    // summarise transcript-level abundance estimates to gene level
    SUMMARIZE_TO_GENE(samplesheet_file, SALMON_QUANT.out.collect())
    
}

workflow TEST {

    // sample channels
    samplesheet_file = file(params.samplesheet)
    samples_ch = channel.from(samplesheet_file)
                    .splitCsv(header: true)
                    .map{ record -> tuple(record.sample, file(record.read1), file(record.read2)) }

    // preprocess reads
    TRIM_READS(samples_ch)

}
