import glob, os
import pandas as pd

configfile: "config/config.2020_06_25.yaml"

samples, = glob_wildcards(config['fastqs'] + '/' + '{sample}_1.fq.gz')
pairs = [1, 2]
ID = config['ID']

rule all:
	input:
		directory("outs/{}/{}".format(config["ID"], config["ref"]["build"])),
		expand('outs/star/{sample}/Aligned.sortedByCoord.out.bam', sample = samples),
		"outs/{}/calls/all.filtered.vcf.gz".format(config["ID"]),
#		expand("outs/qc/{sample}_{pair}_fastqc.zip", sample = samples, pair = pairs),
		"outs/{}/multiqc_report.html".format(config["ID"]) 


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

