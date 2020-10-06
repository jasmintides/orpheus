import glob, os
import pandas as pd

configfile: "config/config.NCI-H1703_ADAR.yaml"

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
		expand('outs/{ID}/star/{sample}/Aligned.sortedByCoord.out.bam', ID = ID, sample = samples),
		## QC report (fastqc, STAR) ##
		"outs/{}/qc/multiqc_report.{}.html".format(config["ID"], config["ID"]),
		## Filtered VCFs ##
		expand("outs/{ID}/calls/filtered/{sample}.vcf.gz", ID = ID, sample = samples)


### include rules ###
include: 'workflow/rules/align.smk'
include: 'workflow/rules/qc.smk'
include: 'workflow/rules/call.smk'

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

