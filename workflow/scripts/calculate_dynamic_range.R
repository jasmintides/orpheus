#!/usr/bin/env Rscript

library(dupRadar)

determine_strand <- function(log){
  stranded <- 0
  if (log$V2[1] > 0.9 & log$V2[2] < 0.9) {
    stranded <- 1
  } else if (log$V2[1] < 0.9 & log$V2[2] > 0.9) {
    stranded <- 2
  }
  return(stranded)
}

strand_log <- read.delim(snakemake@input$strand_log, skip = 4, header = F, sep = ":", stringsAsFactors = F)
  
bamDuprm <- snakemake@input$bam
gtf <- snakemake@input$gtf
stranded <- determine_strand(strand_log)
paired <- as.logical(snakemake@params$is_paired_end)
threads <- snakemake@threads

dm <- analyzeDuprates(bamDuprm, gtf, stranded, paired, threads)
df <- data.frame(sample = snakemake@wildcards$sample,
                 dynrange.all = getDynamicRange(dm)$dynrange.all,
                 dynrange.duprm = getDynamicRange(dm)$dynrange.duprm,
                 stringsAsFactors = F)

write.table(df, as.character(snakemake@output), row.names = F, col.names = T, sep = "\t", quote = F)

##### for aggregating results together -- will revisit later
# for (wildcard in 1:length(snakemake@wildcards$sample)) {
#   bam = grep(wildcards[wildcard], snakemake@input$bam, value = T)
#   strand_log_filename = grep(wildcards[wildcard], snakemake@input$strand_log, value = T)
#   strand_log <- read.delim(strand_log_filename, skip = 4, header = F, sep = ":", stringsAsFactors = F)
#   
#   bamDuprm <- bam
#   gtf <- snakemake@input$gtf
#   stranded <- determine_strand(strand_log)
#   paired <- as.logical(snakemake@input$is_paired_end[wildcard])
#   threads <- snakemake@threads
#   
#   list_of_results[[df]] <- data.frame(sample = snakemake@wildcards$sample[wildcard],
#                                       dynrange.all = getDynamicRange(dm)$dynrange.all,
#                                       dynrange.duprm = getDynamicRange(dm)$dynrange.duprm,
#                                       stringsAsFactors = F)
#   }
# 
# dynamic_range_df <- do.call(rbind, list_of_results)
# rownames(dynamic_range_df) <- NULL