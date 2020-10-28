import glob, os
import pandas as pd

configfile: "config/config.test.yaml"

samples, = glob_wildcards(config['fastqs'] + '/' + '{sample}_1.fq.gz')
pairs = [1, 2]
ID = config['ID']

print(ID)
print(samples)

rule all:
	input:
		## STAR index ##
		directory("outs/{}/{}".format(config["ID"], config["ref"]["build"])),
		## BAM ##
		expand('outs/{ID}/star/{sample}/Aligned.toTranscriptome.out.bam', ID = ID, sample = samples),
		## QC report (fastqc, STAR) ##
		"outs/{}/qc/multiqc_report.{}.html".format(config["ID"], config["ID"]),
		## Isoform quantification ##
#		directory("outs/{}/RSEM_{}".format(config["ID"], config["ref"]["build"])),
#		directory("outs/{}/ref".format(config["ID"])),
		"outs/{}/ref/{}.transcripts.fa".format(config["ID"], config["ref"]["build"]),
		expand("outs/{ID}/RSEM/{sample}.isoforms.results", ID = ID, sample = samples),
		expand("outs/{ID}/RSEM/{sample}.genes.results", ID = ID, sample = samples)
		## Filtered VCFs ##
#		expand("outs/{ID}/calls/filtered/{sample}.vcf.gz", ID = ID, sample = samples)
		

### include rules ###
include: 'workflow/rules/align.smk'
include: 'workflow/rules/qc.smk'
#include: 'workflow/rules/call.smk'
include: 'workflow/rules/quant.smk'

#rule raw_counts:
#	input:
#		gtf = config['ref']['gtf'],
#		bams = expand('outs/STAR/bams/{sample}.Aligned.sortedByCoord.out.bam')
#	output:
#		'outs/counts/Pipeline.Counts.tsv'
#	threads: 
#		2
#	conda: 
#		'workflow/envs/raw_counts.yaml'
#	log: 
#		'logs/raw_counts.log'
#	script:
#		'workflow/scripts/create_counts.R'

