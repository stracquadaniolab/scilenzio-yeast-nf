#!/usr/bin/env Rscript

"scilenzio-summarize-to-gene.R

Usage:
    scilenzio-summarize-to-gene.R [<inputfile>] [--gff=<gff>] [--quant-dir=<quant_dir>] [--counts-from-abundance=<counts_model>] [--output=<output>] 
    
Options:
    --gff=<gff>                             GFF file used to count transcripts.
    --quant-dir=<quant_dir>                 Directory containing salmon results [default: .].
    --counts-from-abundance=<counts_model>  Generate estimated counts using abundance estimates either:
                                            no, scaledTPM, lengthScaledTPM, dtuScaledTPM [default: no].
    --output=<output>                       RDS output file [default: txi-summarized-experiment.rds].
    -h --help                               Show this screen.
    --version                               Show version.
" -> doc

# parsing command line arguments
suppressMessages(library(docopt))
arguments <- docopt(doc, version = "scilenzio-summarize-to-gene.R")

# load required packages
message(">>> Loading packages.")
suppressMessages(library(tidyverse)) # tximport uses readr package if available
suppressMessages(library(tximport))
suppressMessages(library(GenomicFeatures)) # for converting gff file into tx2gene object

# read metadata file
message(">>> Reading metadata file.")
metadata_file <- read_csv(arguments$inputfile)

# prepare vector of file names and associated paths (from salmon-quant) for tximport
message(">>> Preparing vector file names.")
metadata_file$files <- file.path(arguments$quant_dir, metadata_file$Sample, "quant.sf")
files <- metadata_file$files
names(files) <- metadata_file$Sample
message(">>> Done.")

# files <- arguments$quant_dir

# for debugging purposes only
print(arguments)
print(all(file.exists(files)))
print(files)

# obtain annotation
message(">>> Making txdb object.")
txdb <- makeTxDbFromGFF(arguments$gff)

# Associate transcripts with gene IDs for gene-level summarization
# NB: this must have specific column order: 1) transcript ID and 2) gene ID
message(">>> Getting transcript names.")
k <- keys(txdb, keytype = "TXNAME") # get transcript names
message(">>> Making tx2gene object.")
tx2gene <- select(txdb, keys = k, columns="GENEID", keytype = "TXNAME")

# quantify transcripts and summarize counts to gene level
message(">>> Creating txi object.")
txi <- tximport(files, type = "salmon", countsFromAbundance = arguments$counts_model, tx2gene = tx2gene)
print(dim(txi$abundance)) # for debugging/sanity-check purposes

# save matrix to rds file
message(">>> Writing to file.")
saveRDS(txi, file = arguments$output)

message(">>> Finished!")