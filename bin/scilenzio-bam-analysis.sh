#!/bin/bash

# directory containing BAM files
bam_dir="results/bams"

# directory to save results to
output_dir="results/bam-analysis"

# iterate over each BAM file in directory
for bam_file in $bam_dir/*.bam; do
    # extract filename without extension
    filename=$(basename "$bam_file" .bam)
    
    # run samtools flagstat and save output to a text file
    samtools flagstat "$bam_file" > "${output_dir}/${filename}.bam-stats.txt"
done
