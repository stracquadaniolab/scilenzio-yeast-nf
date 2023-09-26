#!/usr/bin/env Rscript

"split-salmon.R

Usage:
    split-salmon.R [<inputfile>] [--quant-dir=<quant_dir>] [--gff=<gff>] [--output=<output>] 
    
Options:
    --gff=<gff>                             GFF file used to get genes in synthetic strain.
    --quant-dir=<quant_dir>                 Directory containing salmon results [default: .].
    --output=<output>                       Output directory [default: split-salmon-quant].
    -h --help                               Show this screen.
    --version                               Show version.
" -> doc

# parsing command line arguments
suppressMessages(library(docopt))
arguments <- docopt(doc, version = "split-salmon.R")

# load required packages
message(">>> Loading packages.")
suppressMessages(library(tidyverse))
suppressMessages(library(rtracklayer)) # to import GFF

# read metadata file
message(">>> Reading metadata file.")
metadata_file <- read_csv(arguments$inputfile)

# prepare vector of file names and associated paths (from salmon-quant) for tximport
message(">>> Preparing vector file names.")
metadata_file$files <- file.path(arguments$quant_dir, metadata_file$Sample, "quant.sf")
files <- metadata_file$files
names(files) <- metadata_file$Sample
message(">>> Done.")
print(files)

# read in synthetic-transcriptome GFF file and get list of gene IDs
gff <- import(arguments$gff)
syn_gene_id <- gff$ID
syn_gene_id <- na.omit(syn_gene_id)

# create output directory if it doesn't exist
output_dir <- arguments$output
if (!dir.exists(output_dir)) {
  dir.create(output_dir, recursive = TRUE)
}

message(">>> Getting synthetic quant.sf files.")
# filter quant.sf files to keep only genes present in synthetic GFF
for (sample_id in names(files)) {

    # load quant.sf file as a dataframe
    file <- files[[sample_id]]
    quant_df <- read.delim(file, header = TRUE, stringsAsFactors = FALSE, sep = "\t")
    print(head(quant_df))

    # filter quant.sf file to keep only genes present in synthetic gene list
    message("Filtering and cleaning quant.sf file.")
    filtered_quant_df <- quant_df %>%
      filter(Name %in% syn_gene_id)

    # remove duplicate synthetic genes (i.e. any rows with "x.transcript:" in Name column)
    filtered_quant_df <- filtered_quant_df[!grepl("x.transcript:", filtered_quant_df$Name), ]

    # remove any "x." or "transcript:" prefixes
    filtered_quant_df$Name <- sub("^(x\\.|transcript:)", "", filtered_quant_df$Name)

    # remove any "_BY4742" suffixes
    filtered_quant_df$Name <- sub("_BY4742$", "", filtered_quant_df$Name)

    # create subdirectory if it doesn't exist
    subdirectory <- file.path(output_dir, paste0("syn-", sample_id))
    if (!dir.exists(subdirectory)) {
      dir.create(subdirectory, recursive = TRUE)
    }

    # define output path
    output_path <- file.path(subdirectory, paste0("quant.sf"))
    print(output_path)

    # write as TSV
    message(">>> Writing to tsv.")
    write.table(filtered_quant_df, output_path, sep = "\t", row.names = FALSE)
    message(">>> TSV created.")

}

message(">>> Getting wildtype quant.sf files.")

# create consensus GFF file
message(">>> Creating consensus GFF file.")
gff <- as.data.frame(gff) %>%
  # remove duplicate synthetic genes (i.e. any rows with 
  # "x.transcript:" in ID column) and remove NAs in ID column
  filter(!grepl("x.transcript:", ID) & !is.na(ID)) %>%
  # strip any IDs containing "x." | "x.transcript:" | "transcript:" 
  # prefixes
  mutate(ID = gsub("^(x.transcript:|x.|transcript:)", "", ID)) %>%
  mutate(seqnames = gsub("chr111", "chr11", seqnames)) %>%
  # remove any "_BY4742" suffixes from gene names
  mutate(ID = sub("_BY4742$", "", ID))

for (sample_id in names(files)) {

    # load quant.sf file as a dataframe
    file <- files[[sample_id]]
    quant_df <- read.delim(file, header = TRUE, stringsAsFactors = FALSE, sep = "\t")
    print(head(quant_df))

    message("Filtering.")
    filtered_quant_df <- quant_df %>%
      # filter out synthetic genes - any rows with "x." | "x.transcript:" | 
      # "transcript:" in "Name" column
      filter(!grepl("^(x\\.|x.transcript:|transcript:)", Name)) %>%
      # remove any "_BY4742" suffixes from gene names
      mutate(Name = sub("_BY4742$", "", Name))

    # filter quant.sf file to keep only genes in consensus GFF file
    filtered_quant_df <- filtered_quant_df[filtered_quant_df$Name %in% gff$ID, ]

    # create subdirectory if it doesn't exist
    subdirectory <- file.path(output_dir, paste0("wt-", sample_id))
    if (!dir.exists(subdirectory)) {
      dir.create(subdirectory, recursive = TRUE)
    }

    # define output path
    output_path <- file.path(subdirectory, paste0("quant.sf"))
    print(output_path)

    # write as TSV
    message(">>> Writing to tsv.")
    write.table(filtered_quant_df, output_path, sep = "\t", row.names = FALSE)
    message(">>> TSV created.")
}

message(">>> Finished!")

# create "wt" GFF file


