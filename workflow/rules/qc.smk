def get_fq1_qc(wildcards):
	return samples.loc[(wildcards.sample), ["fq1"]].dropna()

def get_fq2_qc(wildcards):
	return samples.loc[(wildcards.sample), ["fq2"]].dropna()

rule fastqc_fq1:
	input:
		get_fq1_qc
	output:
		html="outs/{ID}/qc/{sample}_1.html",
		zip=temp("outs/{ID}/qc/{sample}_1_fastqc.zip")
	wrapper:
		"0.31.1/bio/fastqc"

rule fastqc_fq2:
	input:
		get_fq2_qc
	output:
		html=temp("outs/{ID}/qc/{sample}_2.html"),
		zip=temp("outs/{ID}/qc/{sample}_2_fastqc.zip")
	wrapper:
		"0.31.1/bio/fastqc"

rule multiqc:
	input:
		expand("outs/{ID}/qc/{sample}_{group}_fastqc.zip", group = [1, 2], sample = list_of_samples, ID = ID),
		expand("outs/{ID}/star/{sample}/Log.final.out", sample = list_of_samples, ID = ID)
	output:
		"outs/{}/qc/multiqc_report.{}.html".format(config["ID"], config["ID"])
	wrapper:
		"0.51.3/bio/multiqc"
