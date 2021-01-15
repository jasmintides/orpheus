import pandas as pd

configfile: "config/config.test.yaml"

samples = pd.read_table(config["samples"], dtype = str).set_index("sample", drop = False)
list_of_samples = samples["sample"].tolist()
ID = config['ID']
template_sample = list_of_samples[0]

rule all:
	input:
		## QC report (fastqc, STAR) ##
		"outs/{}/qc/multiqc_report.{}.html".format(config["ID"], config["ID"]),
		## Aggregated gene raw and TPM counts ##
		"outs/{}/star/raw_counts.tsv".format(config["ID"]),
		"outs/{}/RSEM/tpm_counts.tsv".format(config["ID"])

rule all_2:
	input:
		## QC report (fastqc, STAR) ##
		"outs/{}/qc/multiqc_report.{}.html".format(config["ID"], config["ID"]),
		## Aggregated gene raw and TPM counts ##
		"outs/{}/star/raw_counts.tsv".format(config["ID"]),
		"outs/{}/RSEM/tpm_counts.tsv".format(config["ID"]),
		## Per sample variant calls ##
		expand("outs/{ID}/calls/filtered/{sample}.vcf.gz", ID = ID, sample = list_of_samples)

#rule all:
#	input:
		## STAR index ##
#		directory("outs/{}/{}".format(config["ID"], config["ref"]["build"])),
		## BAM ##
#		expand('outs/{ID}/star/{sample}/Aligned.toTranscriptome.out.bam', ID = ID, sample = samples),
		## QC report (fastqc, STAR) ##
#		"outs/{}/qc/multiqc_report.{}.html".format(config["ID"], config["ID"]),
		## Isoform quantification ##
#		"outs/{}/ref/{}.transcripts.fa".format(config["ID"], config["ref"]["build"]),
#		expand("outs/{ID}/RSEM/{sample}.isoforms.results", ID = ID, sample = samples),
#		expand("outs/{ID}/RSEM/{sample}.genes.results", ID = ID, sample = samples),
#		"outs/{}/star/gene_ids.txt".format(config["ID"]),
#		expand("outs/{ID}/star/{sample}/{sample}.counts", ID = ID, sample = samples),
#		"outs/{}/star/raw_counts.tsv".format(config["ID"]),
#		"outs/{}/RSEM/tpm_counts.tsv".format(config["ID"])
		## Filtered VCFs ##
#		expand("outs/{ID}/calls/filtered/{sample}.vcf.gz", ID = ID, sample = samples)
		

### include rules ###
include: 'workflow/rules/qc.smk'
include: "workflow/rules/00_common.smk"
include: "workflow/rules/01_trim.smk"
include: 'workflow/rules/02_align.smk'
include: "workflow/rules/call.smk"
include: "workflow/rules/quant.smk"

