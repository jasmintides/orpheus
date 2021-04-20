import pandas as pd

configfile: "config/config.test.yaml"

samples = pd.read_table(config["samples"], dtype = str).set_index("sample", drop = False)
list_of_samples = samples["sample"].tolist()
ID = config['ID']
outpath = config['outpath']
template_sample = list_of_samples[0]
build = config["ref"]["build"]

chr_dict = {"chunk_01": "-L chr1", "chunk_02": "-L chr2", "chunk_03": "-L chr3",
	"chunk_04": "-L chr4", "chunk_05": "-L chr5", "chunk_06": "-L chr6",
	"chunk_07": "-L chr7", "chunk_08": "-L chr8", "chunk_09": "-L chr9",
	"chunk_10": "-L chr10", "chunk_11": "-L chr11", "chunk_12": "-L chr12 -L chr13",
	"chunk_13": "-L chr14 -L chr15", "chunk_14": "-L chr16 -L chr17 -L chr18",
	"chunk_15": "-L chr19 -L chr20 -L chr21 -L chr22", "chunk_16": "-L chrX -L chrY"}
chr_chunks = list(chr_dict.keys())
chr_params = list(chr_dict.values())

#if config["run_sra_pipeline"] is True:
#	subworkflow sra_pipeline:
#		workdir:
#			config["sra_workdir"]
#		snakefile:
#			config["sra_Snakefile"]
#		configfile:
#			config["sra_configfile"]
#
#	sample_sheet = "{}/{}/{}_orpheus_sample_sheet.tsv".format(config["sra_outpath"], 
#			config["sra_BioProject_ID"], 
#			config["sra_BioProject_ID"])

#	rule all_sra:
#		input:
#			sra_pipeline(sample_sheet),
#			"{}/{}/qc/multiqc_report.{}.html".format(outpath, ID, ID),
#			"{}/{}/star/raw_counts.tsv".format(outpath, ID),
#			"{}/{}/SummExp/{}.genes_SummExp.Rds".format(outpath, ID, ID),
#			"{}/{}/SummExp/{}.transcripts_SummExp.Rds".format(outpath, ID, ID)

#else:
#	samples = pd.read_table(config["samples"], dtype = str).set_index("sample", drop = False)
#	list_of_samples = samples["sample"].tolist()
#	template_sample = list_of_samples[0]

#	rule all_no_sra:
#		input:
#			"{}/{}/qc/multiqc_report.{}.html".format(outpath, ID, ID),
#			"{}/{}/star/raw_counts.tsv".format(outpath, ID),
#			"{}/{}/SummExp/{}.genes_SummExp.Rds".format(outpath, ID, ID),
#			"{}/{}/SummExp/{}.transcripts_SummExp.Rds".format(outpath, ID, ID),
#			expand("{outpath}/{ID}/rseqc/{sample}.dynamic_range.txt",
#				outpath = outpath, ID = ID, sample = list_of_samples)

rule all:
	input:
		## QC report (fastqc, STAR) ##
		"{}/{}/qc/multiqc_report.{}.html".format(outpath, ID, ID),
		## Aggregated gene raw and TPM counts ##
		"{}/{}/star/raw_counts.tsv".format(outpath, ID),
		"{}/{}/SummExp/{}.genes_SummExp.Rds".format(outpath, ID, ID),
		"{}/{}/SummExp/{}.transcripts_SummExp.Rds".format(outpath, ID, ID),
		expand("{outpath}/{ID}/rseqc/{sample}.dynamic_range.txt",
			outpath = outpath, ID = ID, sample = list_of_samples)

rule quant_call_variants:
	input:
		## QC report (fastqc, STAR) ##
		"{}/{}/qc/multiqc_report.{}.html".format(outpath, ID, ID),
		## Aggregated gene raw and TPM counts ##
		"{}/{}/star/raw_counts.tsv".format(outpath, ID),
		"{}/{}/RSEM/tpm_counts.tsv".format(outpath, ID),
		"{}/{}/RSEM/tpm_isoforms_counts.tsv".format(outpath, ID),
		## Per sample variant calls ##
		expand("{outpath}/{ID}/final/{sample}.vcf.gz", 
			outpath = outpath, ID = ID, sample = list_of_samples)

rule multiqc_only:
	input:
		"{}/{}/qc/multiqc_report.{}.html".format(outpath, ID, ID),

### include rules ###
include: "workflow/rules/00_common.smk"
include: "workflow/rules/01_trim.smk"
include: 'workflow/rules/02_align.smk'
include: "workflow/rules/call.smk"
include: "workflow/rules/quant.smk"
include: "workflow/rules/rRNA_filter.smk"
include: 'workflow/rules/qc.smk'
include: 'workflow/rules/calc_dynamic_range.smk'
include: 'workflow/rules/write_se.smk'
