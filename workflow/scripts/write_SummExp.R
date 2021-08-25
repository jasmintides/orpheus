#!/usr/bin/env Rscript

library(tidyverse)
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

outfile <- as.character(snakemake@output)
outpath <- dirname(outfile)
if  (!dir.exists(outpath)) {
  dir.create(outpath)
}

# Create SummarizedExperiments
transcripts_se <- SummarizedExperiment::SummarizedExperiment(
	assays = list(transcripts_expected = transcripts_expected,
		transcripts_tpm = transcripts_tpm),
	colData = final_coldata)

saveRDS(transcripts_se, file = as.character(snakemake@output))
