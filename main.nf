// enabling nextflow DSL v2
nextflow.enable.dsl=2

// process BUILD_REFERENCE_TRANSCRIPTOME {

//     publishDir "${params.resultsDir}/combined-genome", mode: 'copy', overwrite: true

//     input:
//         path(wt_reference)
//         path(syn_reference)
//         path(wt_annotation)
//         path(syn_annotation)

//     output:
//         path("wt-syn-chr11-ref.fasta"), emit: reference
//         path("wt-syn-chr11-ref.gff"), emit: annotation

//     """
//         # Replace fasta header of synthetic reference to "chr18" to indicate synthetic contig
//         # concatenate wt and synthetic reference (in that order)
//         cat ${wt_reference} <(sed 's/^>chr11/>chr18/' ${syn_reference}) > wt-syn-chr11-ref.fasta

//         # Change chromosome and gene names in synthetic GFF
//         # Replace "chr11" with "chr18" to indicate synthetic chromosome
//         sed 's/chr11/chr18/g' ${syn_annotation} > tmp-syn.gff
//         # Add "x." prefix to synthetic GFF "ID=" column to indicate synthetic genes
//         awk -F'\t' -v OFS='\t' '{ sub(/ID=/, "ID=x.", $9); sub(/Name=/, "Name=x.", $9); print }' tmp-syn.gff > modified-syn-chr11.gff

//         # concatenate wt and synthetic annotation
//         cat ${wt_annotation} modified-syn-chr11.gff > wt-syn-chr11-ref.gff
//     """

//     stub:
//     """
//         touch wt-syn-chr11-ref.fasta
//         touch wt-syn-chr11-ref.gff
//     """

// }

process TRIM_READS {

    tag "${sample}"

    maxForks 1

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

    maxForks 1

    input:
        path(genome_fasta)
        path(annotation_gff)

    output:
        path("genome-index")

    script:
    """
        STAR \\
            --runThreadN ${task.cpus} \\
            --runMode genomeGenerate \\
            --genomeDir genome-index \\
            --genomeFastaFiles ${genome_fasta} \\
            --sjdbGTFfile ${annotation_gff} \\
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

    maxForks 1

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


process GFFREAD_GET_WT_TRANSCRIPTOME {

    maxForks 1

    publishDir "${params.resultsDir}/wt-syn-transcriptome/", mode: 'copy', overwrite: true

    input:
        path(genome_fasta)
        path(annotation_gff)

    output:
        path("wt-syn-transcriptome.fa"), emit: transcriptome
        path("wt-syn-transcriptome.gff"), emit: annotation

    script:
    """
        gffread \\
            -g ${genome_fasta} \\
            -o wt-syn-transcriptome.gff \\
            -w wt-syn-transcriptome.fa \\
            -v \\
            ${annotation_gff}
    """

    stub:
    """
        touch wt-syn-transcriptome.fa
        touch wt-syn-transcriptome.gff
    """

}

process SALMON_INDEX {

    maxForks 1

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
            ${params.salmon.index.args}
    """

    stub:
    """
        mkdir transcriptome-index
    """

}

process SALMON_QUANT {

    tag "${sample}"

    maxForks 1

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

    publishDir "${params.resultsDir}/tximport/", mode: 'copy', overwrite: true

    input:
        path(sample_sheet)
        path(annotation_gff)
        path(salmon_results, stageAs: 'salmon-quant/quant*.sf')

    output:
        path("txi-summarized-experiment.rds")

    script:
    """
        summarize-to-gene.R \\
            ${sample_sheet} \\
            --gff ${annotation_gff} \\
            --quant-dir ${salmon_results} \\
            --counts-from-abundance ${params.summarize_to_gene.counts_from_abundance} \\
            --output txi-summarized-experiment.rds
    """

    stub:
    """
        touch txi-summarized-experiment.rds 
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
                    .map{ record -> tuple(record.Sample, file(record.read1), file(record.read2)) }

    // // build reference
    // BUILD_REFERENCE_TRANSCRIPTOME(file(params.reference.wt), file(params.reference.syn), file(params.annotation.wt), file(params.annotation.syn))

    // // preprocess reads
    // TRIM_READS(samples_ch)

    // perform alignment
    // GENERATE_GENOME_INDEX(file(params.genome.reference), file(params.genome.annotation))
    // ALIGN_READS(GENERATE_GENOME_INDEX.out, TRIM_READS.out.fastq)

    // generate transcriptome fasta
    GFFREAD_GET_WT_TRANSCRIPTOME(file(params.genome.reference), file(params.genome.annotation))

    // create transcriptome index
    SALMON_INDEX(GFFREAD_GET_WT_TRANSCRIPTOME.out.transcriptome, file(params.genome.reference))

    // // estimate transcript level abundance
    // SALMON_QUANT(SALMON_INDEX.out, TRIM_READS.out.fastq)

}

workflow ANALYSIS {

    // summarise transcript-level abundance estimates to gene level
    SUMMARIZE_TO_GENE(file(params.samplesheet), file(params.summarize_to_gene.annotation), file(params.summarize_to_gene.quant_dir))

}
