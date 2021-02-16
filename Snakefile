import pandas as pd

configfile: "config/config.test.yaml"

samples = pd.read_table(config["samples"], dtype = str).set_index("sample", drop = False)
list_of_samples = samples["sample"].tolist()
ID = config['ID']
outpath = config['outpath']
template_sample = list_of_samples[0]

chr_dict = {"chunk_01": "-L chr1", "chunk_02": "-L chr2", "chunk_03": "-L chr3",
	"chunk_04": "-L chr4", "chunk_05": "-L chr5", "chunk_06": "-L chr6",
	"chunk_07": "-L chr7", "chunk_08": "-L chr8", "chunk_09": "-L chr9",
	"chunk_10": "-L chr10", "chunk_11": "-L chr11", "chunk_12": "-L chr12 -L chr13",
	"chunk_13": "-L chr14 -L chr15", "chunk_14": "-L chr16 -L chr17 -L chr18",
	"chunk_15": "-L chr19 -L chr20 -L chr21 -L chr22", "chunk_16": "-L chrX -L chrY"}
chr_chunks = list(chr_dict.keys())
chr_params = list(chr_dict.values())

rule all:
	input:
		## QC report (fastqc, STAR) ##
		"{}/{}/qc/multiqc_report.{}.html".format(outpath, ID, ID),
		## Aggregated gene raw and TPM counts ##
		"{}/{}/star/raw_counts.tsv".format(outpath, ID),
		"{}/{}/RSEM/genes.expected_counts.tsv".format(outpath, ID),
		"{}/{}/RSEM/genes.tpm_counts.tsv".format(outpath, ID),
		"{}/{}/RSEM/isoforms.expected_counts.tsv".format(outpath, ID),
		"{}/{}/RSEM/isoforms.tpm_counts.tsv".format(outpath, ID)

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

rule all_2:
	input:
		expand("{outpath}/{ID}/sortmerna/{sample}.clean.fastq",
			outpath = outpath, ID = ID, sample = list_of_samples)


### include rules ###
include: "workflow/rules/00_common.smk"
include: "workflow/rules/01_trim.smk"
include: 'workflow/rules/02_align.smk'
include: "workflow/rules/call.smk"
include: "workflow/rules/quant.smk"
include: "workflow/rules/rRNA_filter.smk"
include: 'workflow/rules/qc.smk'
