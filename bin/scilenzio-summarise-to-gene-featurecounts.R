#!/usr/bin/env Rscript

"scilenzio-summarise-to-gene-featurecounts.R

Usage:
    scilenzio-summarise-to-gene-featurecounts.R [--samplesheet=<samplesheet>] [--gff=<gff>] [--featurecounts-tpm-dir=<featurecounts_tpm_dir>] [--counts-from-abundance=<counts_model>] [--output=<output>] 
    
Options:
    --samplesheet=<samplesheet>                         Samplesheet containing sample ID information.
    --gff=<gff>                                         GFF file used to count transcripts.
    --featurecounts-tpm-dir=<featurecounts_tpm_dir>     Directory containing featureCounts files [default: .].
    --counts-from-abundance=<counts_model>              Generate estimated counts using abundance estimates either:
                                                        no, scaledTPM, lengthScaledTPM, dtuScaledTPM [default: no].
    --output=<output>                                   RDS output file [default: 
                                                        txi-featurecounts-summarized-experiment.rds].
    -h --help                                           Show this screen.
    --version                                           Show version.
" -> doc

# parsing command line arguments
suppressMessages(library(docopt))
arguments <- docopt(doc, version = "scilenzio-summarise-to-gene-featurecounts.R")

# load required packages
suppressMessages(library(readr)) # tximport uses readr package if available
suppressMessages(library(dplyr))
suppressMessages(library(tximport))
suppressMessages(library(GenomicFeatures)) # for converting gff file into tx2gene object

# read files
samplesheet <- read_csv(arguments$samplesheet)
gff <- rtracklayer::import(arguments$gff)

# prepare vector of file names and associated paths from featureCounts for tximport
samplesheet$files <- file.path(arguments$featurecounts_tpm_dir, 
                                 paste(samplesheet$Sample, 
                                       "featureCounts-tpm.txt", 
                                       sep = "."))
files <- samplesheet$files
names(files) <- samplesheet$Sample

# for debugging purposes only
print(arguments)
print(all(file.exists(files)))
print(files)

# create tx2gene object (TXNAME, GENEID)
# because our GFF file is formatted differently, we need to create
# our tx2gene object manually
tx2gene <- as.data.frame(gff) %>%
  dplyr::filter(., type == "gene") %>%
  dplyr::select(., Name, ID) %>%
  dplyr::rename(., TXNAME = Name) %>%
  dplyr::rename(., GENEID = ID)

# quantify transcripts and summarize counts to gene level
txi <- tximport(files,
                type = "none",
                countsFromAbundance = "no",
                tx2gene = tx2gene,
                geneIdCol = "Geneid", # N.B. refers to column in featureCounts file
                txIdCol = "transcript_id", # N.B. refers to column in featureCounts file
                abundanceCol = "TPM",
                countsCol = "NumReads",
                lengthCol = "Length",
                importer=function(x) read_tsv(x)) # specify the importer function on how to read lines

print(dim(txi$abundance)) # for debugging/sanity-check purposes

# save matrix to rds file
saveRDS(txi, file = arguments$output)