#!/usr/bin/env Rscript

library(tidyverse)
library(biomaRt)
library(SummarizedExperiment)
# library(HDF5Array)

genes_expected <- read.delim(snakemake@input$genes_expected, header = T, row.names = 1, sep = "\t", stringsAsFactors = F)
genes_tpm <- read.delim(snakemake@input$genes_tpm, header = T, row.names = 1, sep = "\t", stringsAsFactors = F)
transcripts_expected <- read.delim(snakemake@input$transcripts_expected, header = T, row.names = 1, sep = "\t", stringsAsFactors = F)
transcripts_tpm <- read.delim(snakemake@input$transcripts_tpm, header = T, row.names = 1, sep = "\t", stringsAsFactors = F)
sample_sheet <- read.delim(snakemake@input$sample_sheet, header = T, sep = "\t", stringsAsFactors = F)

# counts_dir <- "../../outs/20XX_11_11_test/RSEM/"
# genes_expected <- read.delim(paste0(counts_dir, "genes.expected_counts.tsv"), header = T, row.names = 1, sep = "\t", stringsAsFactors = F)
# genes_tpm <- read.delim(paste0(counts_dir, "genes.tpm_counts.tsv"), header = T, row.names = 1, sep = "\t", stringsAsFactors = F)
# transcripts_expected <- read.delim(paste0(counts_dir, "transcripts.expected_counts.tsv"), header = T, row.names = 1, sep = "\t", stringsAsFactors = F)
# transcripts_tpm <- read.delim(paste0(counts_dir, "transcripts.tpm_counts.tsv"), header = T, row.names = 1, sep = "\t", stringsAsFactors = F)
# 
# metadata_name <- "../../resources/sample_sheet.tsv"
# sample_sheet <- read.delim(metadata_name, header = T, sep = "\t", stringsAsFactors = F)

# Read in MultiQC files
# MultiQC_dir <- "../../outs/20XX_11_11_test/qc/multiqc_report.20XX_11_11_test_data/"
#fastqc <- read.delim(paste0(MultiQC_dir, "multiqc_fastqc.txt"), header = T, sep = "\t", stringsAsFactors = F) %>%
fastqc <- read.delim(snakemake@input$fastqc, header = T, sep = "\t", stringsAsFactors = F) %>% 
  dplyr::mutate(Sample = gsub("_[12]$", "", Sample)) %>% 
  dplyr::group_by(Sample) %>% 
  dplyr::summarise(across(everything(), function(x) {paste(x, collapse = "|")}))
colnames(fastqc)[-1] <- paste0("fastqc.", colnames(fastqc)[-1])
# star <- read.delim(paste0(MultiQC_dir, "multiqc_star.txt"), header = T, sep = "\t", stringsAsFactors = F)
star <- read.delim(snakemake@input$star, header = T, sep = "\t", stringsAsFactors = F)
colnames(star)[-1] <- paste0("star.", colnames(star)[-1])
# sortmerna <- read.delim(paste0(MultiQC_dir, "multiqc_sortmerna.txt"), header = T, sep = "\t", stringsAsFactors = F)
sortmerna <- read.delim(snakemake@input$sortmerna, header = T, sep = "\t", stringsAsFactors = F)
colnames(sortmerna)[-1] <- paste0("sortmerna.", colnames(sortmerna)[-1])
# rseqc <- read.delim(paste0(MultiQC_dir, "multiqc_rseqc_infer_experiment.txt"), header = T, sep = "\t", stringsAsFactors = F)
rseqc <- read.delim(snakemake@input$rseqc, header = T, sep = "\t", stringsAsFactors = F)
colnames(rseqc)[-1] <- paste0("rseqc.", colnames(rseqc)[-1])

# Join sample sheet and QC files together
final_coldata <- dplyr::left_join(sample_sheet, fastqc, by = c("sample" = "Sample")) %>% 
  dplyr::left_join(star, by = c("sample" = "Sample"), suffix = c(".fastqc", ".star")) %>% 
  dplyr::left_join(sortmerna, by = c("sample" = "Sample")) %>% 
  dplyr::left_join(rseqc, by = c("sample" = "Sample"))

## Important - set parameter in config file as "Human" or "Mouse" to determine whether to get human or mouse Ensembl IDs
organism <- snakemake@params$organism
# organism <- "Human"

outfile <- snakemake@output$genes_SummExp
outpath <- dirname(outfile)
if  (!dir.exists(outpath)) {
  dir.create(outpath)
}

# Given organism specified in config file, get gene annotation for given matrix
get_gene_anno <- function(organism, matrix) {
  # Determine whether to get datasets from human or mouse from Biomart
  if (organism == "Human") {
    ensembl <- useMart("ensembl", dataset = "hsapiens_gene_ensembl")
    organism_symbol <- "hgnc_symbol"
    } else if (organism == "Mouse") {
      ensembl <- useMart("ensembl", dataset = "mmusculus_gene_ensembl")
      organism_symbol <- "mgi_symbol"
      }
  # If rownames of gene counts matrix start with Ensembl ID, retrieve gene symbols
  if (any(grepl("^ENSG|^ENSMUSG", rownames(matrix)))) {
    # If Ensembl IDs contain ".", retrieve Ensembl gene ID with version number
    if (any(grepl("\\.", rownames(matrix)))) {
      ensembl_entrez_map <- getBM(attributes = c("ensembl_gene_id_version", "ensembl_gene_id",organism_symbol),
                                  filter = "ensembl_gene_id_version",
                                  value = rownames(matrix),
                                  mart = ensembl,
                                  useCache = FALSE) %>% 
        group_by(ensembl_gene_id_version) %>% 
        slice_head(n = 1) %>% 
        dplyr::left_join(data.frame(ensembl_gene_id_version = rownames(matrix)), by = c("ensembl_gene_id_version"))
      # If Ensembl IDs do not contain ".", retrieve stable Ensembl gene ID
      } else {
        ensembl_entrez_map <- getBM(attributes = c("ensembl_gene_id", organism_symbol),
                                    filter = "ensembl_gene_id",
                                    value = rownames(matrix),
                                    mart = ensembl,
                                    useCache = FALSE) %>% 
          group_by(ensembl_gene_id) %>% 
          slice_head(n = 1) %>% 
          dplyr::left_join(data.frame(ensembl_gene_id = rownames(matrix)), by = c("ensembl_gene_id"))
      }
  }
    # Check if rownames of matrix are UCSC IDs
    else if (any(grepl("^uc.", rownames(matrix)))) {
      ensembl_entrez_map <- getBM(attributes = c("ucsc","ensembl_gene_id", organism_symbol),
                                  filter = "ucsc",
                                  value = rownames(matrix),
                                  mart = ensembl,
                                  useCache = FALSE) %>% 
        group_by(ensembl_gene_id) %>% 
        slice_head(n = 1) %>% 
        dplyr::left_join(data.frame(ucsc = rownames(matrix)), by = c("ucsc"))
      }
    # If rownames are already gene symbols, get Ensembl IDs
  else {
      ensembl_entrez_map <- getBM(attributes = c("ensembl_gene_id", organism_symbol),
                                  filter = organism_symbol,
                                  value = rownames(matrix),
                                  mart = ensembl,
                                  useCache = FALSE) %>% 
        dplyr::rename(gene_symbol = 2) %>% 
        dplyr::left_join(data.frame(gene_symbol = rownames(matrix)), by = c("gene_symbol"))
      }
  return(ensembl_entrez_map)
  }

# Get gene annotation from expected counts
gene_anno <- get_gene_anno(organism, genes_expected)

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

# Given organism specified in config file, get gene annotation for given matrix
filter_gene_counts <- function(organism, matrix) {
  if (organism == "Human") {
    organism_symbol <- "hgnc_symbol"
  } else if (organism == "Mouse") {
    organism_symbol <- "mgi_symbol"
  }
  if (any(grepl("^ENSG|^ENSMUSG", rownames(matrix)))) {
    if (any(grepl("\\.", rownames(matrix)))) {
      filtered_matrix <- matrix[gene_anno$ensembl_gene_id_version,]
    } else {
      filtered_matrix <- matrix[gene_anno$ensembl_gene_id, ]    }
  }
  else if (any(grepl("^uc.", rownames(matrix)))) {
    filtered_matrix <- matrix[gene_anno$ucsc, ]
    }
  else {
    genes <- as.character(unlist(gene_anno[,organism_symbol]))
    filtered_matrix <- matrix[genes[genes != ""], ]
  }
  return(filtered_matrix)
}

filtered_genes_expected <- filter_gene_counts(organism, genes_expected)
filtered_genes_tpm <- filter_gene_counts(organism, genes_tpm)

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
genes_se <- SummarizedExperiment::SummarizedExperiment(assays = list(genes_expected = filtered_genes_expected,
                                                                     genes_tpm = filtered_genes_tpm),
                                                       colData = final_coldata,
                                                       rowData = gene_anno)
transcripts_se <- SummarizedExperiment::SummarizedExperiment(assays = list(transcripts_expected = filtered_transcripts_expected,
                                                                        transcripts_tpm = filtered_transcripts_tpm),
                                                          colData = final_coldata,
                                                          rowData = transcript_anno)

saveRDS(genes_se, file = as.character(snakemake@output$genes_SummExp))
saveRDS(transcripts_se, file = as.character(snakemake@output$transcripts_SummExp))
