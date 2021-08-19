#!/usr/bin/env Rscript

library(tidyverse)
library(biomaRt)
library(SummarizedExperiment)

transcripts_expected <- read.delim(as.character(snakemake@input$expected), header = T, row.names = 1, sep = "\t", stringsAsFactors = F)
transcripts_tpm <- read.delim(as.character(snakemake@input$tpm), header = T, row.names = 1, sep = "\t", stringsAsFactors = F)
sample_sheet <- read.delim(as.character(snakemake@input$sample_sheet), header = T, sep = "\t", stringsAsFactors = F)

fastqc <- read.delim(as.character(snakemake@input$fastqc), header = T, sep = "\t", stringsAsFactors = F) %>% 
  dplyr::mutate(Sample = gsub("_[12]$", "", Sample)) %>% 
  dplyr::group_by(Sample) %>% 
  dplyr::summarise(across(everything(), function(x) {paste(x, collapse = "|")}))
colnames(fastqc)[-1] <- paste0("fastqc.", colnames(fastqc)[-1])

# Join sample sheet and QC files together
final_coldata <- dplyr::left_join(sample_sheet, fastqc, by = c("sample" = "Sample"))

## Important - set parameter in config file as "Human" or "Mouse" to determine whether to get human or mouse Ensembl IDs
organism <- snakemake@params$organism

outfile <- as.character(snakemake@output)
outpath <- dirname(outfile)
if  (!dir.exists(outpath)) {
  dir.create(outpath)
}
# Given organism specified in config file, get isoform annotation for given matrix
get_transcript_anno <- function(organism, matrix) {
  # Determine whether to get datasets from human or mouse from Biomart
  if (organism == "Human") {
    ensembl <- useMart("ensembl", dataset = "hsapiens_gene_ensembl")
  } else if (organism == "Mouse") {
    ensembl <- useMart("ensembl", dataset = "mmusculus_gene_ensembl")
  }
  # If rownames of transcript counts matrix start with Ensembl ID, retrieve transcript symbols
  if (any(grepl("^ENST|^ENSMUST", rownames(matrix)))) {
    # If Ensembl IDs contain ".", retrieve Ensembl transcript ID with version number
    if (any(grepl("\\.", rownames(matrix)))) {
      ensembl_entrez_map <- getBM(attributes = c("ensembl_transcript_id_version", "ensembl_transcript_id", "external_transcript_name"),
                                  filter = "ensembl_transcript_id_version",
                                  value = rownames(matrix),
                                  mart = ensembl,
                                  useCache = FALSE) %>% 
        group_by(ensembl_transcript_id_version) %>% 
        slice_head(n = 1) %>% 
        dplyr::left_join(data.frame(ensembl_transcript_id_version = rownames(matrix)), by = c("ensembl_transcript_id_version"))
      # If Ensembl IDs do not contain ".", retrieve stable Ensembl transcript ID
    } else {
      ensembl_entrez_map <- getBM(attributes = c("ensembl_transcript_id", "ensembl_transcript_id", "external_transcript_name"),
                                  filter = "ensembl_transcript_id",
                                  value = rownames(matrix),
                                  mart = ensembl,
                                  useCache = FALSE) %>% 
        group_by(ensembl_transcript_id) %>% 
        slice_head(n = 1) %>% 
        dplyr::left_join(data.frame(ensembl_transcript_id = rownames(matrix)), by = c("ensembl_transcript_id"))
    }
  }
    # Check if rownames of matrix are UCSC IDs
    else if (any(grepl("^uc.", rownames(matrix)))) {
      ensembl_entrez_map <- getBM(attributes = c("ucsc", "ensembl_transcript_id", "external_transcript_name"),
                                  filter = "ucsc",
                                  value = rownames(matrix),
                                  mart = ensembl,
                                  useCache = FALSE) %>% 
        group_by(ensembl_transcript_id) %>% 
        slice_head(n = 1) %>% 
        dplyr::left_join(data.frame(ucsc = rownames(matrix)), by = c("ucsc"))
    }
    # If rownames are already transcript symbols, get Ensembl IDs
  else {
    ensembl_entrez_map <- getBM(attributes = c("ensembl_transcript_id", "external_transcript_name"),
                                filter = "external_transcript_name",
                                value = rownames(matrix),
                                mart = ensembl,
                                useCache = FALSE) %>% 
      dplyr::rename(transcript_name = 2) %>% 
      dplyr::left_join(data.frame(transcript_name = rownames(matrix)), by = c("transcript_name"))
  }
  return(ensembl_entrez_map)
}

# Get transcript annotation from expected counts
transcript_anno <- get_transcript_anno(organism, transcripts_expected)

filter_transcript_counts <- function(matrix) {
  # If rownames of transcript counts matrix start with Ensembl ID, retrieve transcript symbols
  if (any(grepl("^ENST|^ENSMUST", rownames(matrix)))) {
    # If Ensembl IDs contain ".", retrieve Ensembl transcript ID with version number
    if (any(grepl("\\.", rownames(matrix)))) {
      filtered_matrix <- matrix[transcript_anno$ensembl_transcript_id_version,]
      # If Ensembl IDs do not contain ".", retrieve stable Ensembl transcript ID
    } else {
      filtered_matrix <- matrix[transcript_anno$ensembl_transcript_id,]
    }
  }
  # Check if rownames of matrix are UCSC IDs
  else if (any(grepl("^uc.", rownames(matrix)))) {
    ucsc_transcripts <- transcript_anno$ucsc
    filtered_matrix <- matrix[ucsc_transcripts[ucsc_transcripts != ""],]
  }
  # If rownames are already transcript symbols, get Ensembl IDs
  else {
    transcripts <- transcript_anno$external_transcript_name
    filtered_matrix <- matrix[transcripts[transcripts != ""], ]
  }
  return(filtered_matrix)
}

filtered_transcripts_expected <- filter_transcript_counts(transcripts_expected)
filtered_transcripts_tpm <- filter_transcript_counts(transcripts_tpm)

# Create SummarizedExperiments
transcripts_se <- SummarizedExperiment::SummarizedExperiment(assays = list(transcripts_expected = filtered_transcripts_expected,
                                                                        transcripts_tpm = filtered_transcripts_tpm),
                                                          colData = final_coldata,
                                                          rowData = transcript_anno)

saveRDS(transcripts_se, file = as.character(snakemake@output))
