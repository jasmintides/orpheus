def get_fq1_qc(wildcards):
	return samples.loc[(wildcards.sample), ["fq1"]].dropna()

def get_fq2_qc(wildcards):
	if not is_single_end(wildcards.sample):
		return samples.loc[(wildcards.sample), ["fq2"]].dropna()

rule fastqc_fq1:
	input:
		get_fq1_qc
	output:
		html="{outpath}/{ID}/qc/{sample}_1.html",
		zip=temp("{outpath}/{ID}/qc/{sample}_1_fastqc.zip")
	wrapper:
		"0.31.1/bio/fastqc"

rule fastqc_fq2:
	input:
		get_fq2_qc
	output:
		html=temp("{outpath}/{ID}/qc/{sample}_2.html"),
		zip=temp("{outpath}/{ID}/qc/{sample}_2_fastqc.zip")
	wrapper:
		"0.31.1/bio/fastqc"

rule multiqc:
	input:
		expand("{outpath}/{ID}/qc/{sample}_1_fastqc.zip", 
			outpath = outpath, ID = ID, sample = list_of_samples),
		expand("{outpath}/{ID}/qc/{sample}_2_fastqc.zip", 
			outpath = outpath, ID = ID,
			sample = [f for f in list_of_samples if is_paired_end(f)]),
		expand("{outpath}/{ID}/star/{sample}/Log.final.out", 
			outpath = outpath, ID = ID, sample = list_of_samples),
		expand("{outpath}/{ID}/sortmerna/{sample}.rRNA.log", 
			outpath = outpath, ID = ID, sample = list_of_samples),
		expand("{outpath}/{ID}/rseqc/{sample}.infer_experiment.txt",
			outpath = outpath, ID = ID, sample = list_of_samples)
	output:
		"{}/{}/qc/multiqc_report.{}.html".format(outpath, ID, ID)
	wrapper:
		"0.51.3/bio/multiqc"
