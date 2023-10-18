#!/usr/bin/env Rscript

"scilenzio-featurecounts-tpm.R

Usage:
    scilenzio-featurecounts-tpm.R <input> [--gff=<gff>] [--output=<output>]
    
Options:
    --gff=<gff>                             GFF file.
    --output=<output>                       CSV featureCounts output file with additional RPK and TPM counts.
    -h --help                               Show this screen.
    --version                               Show version.
" -> doc

# parsing command line arguments
suppressMessages(library(docopt))
arguments <- docopt(doc, version = "scilenzio-featurecounts-tpm.R")

# load required packages
suppressMessages(library(dplyr))
suppressMessages(library(readr))

# load data
gff <- rtracklayer::import(arguments$gff)
counts <- read_delim(arguments$input, delim = "\t", skip = 1)

# calculate TPMs
counts <- counts %>%
  dplyr::rename(NumReads = 7) %>%
  dplyr::select(Geneid, Length, NumReads) %>%
  # calculate RPK
  mutate(RPK = NumReads/(Length/1000))
# calculate PM scaling factor
PM <- sum(counts$RPK)/1e6
# calculate TPM
counts <- counts %>%
  mutate(TPM = RPK / PM)

# create map of gene IDs to transcript IDs
trID.map <- as.data.frame(gff) %>%
  dplyr::filter(., type == "gene") %>%
  dplyr::select(., ID, Name)

# map featurecounts Geneids to GFF gene IDs
# this is required by tximport so it can correctly map to the 
# first column of its tx2gene object containing transcript IDs
# counts <- merge(counts, trID.map, by.x = "Geneid", by.y = "gene_id", all = TRUE)
counts <- left_join(counts, trID.map, by = c("Geneid" = "ID")) %>%
  dplyr::distinct(Geneid, .keep_all = TRUE) %>%
  dplyr::rename(., transcript_id = Name)

# write to file
write_tsv(counts, arguments$output) # need to save as TSV file for tximport