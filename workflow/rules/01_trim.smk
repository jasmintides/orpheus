def get_fastq(wildcards):
	return samples.loc[(wildcards.sample), ["fq1", "fq2"]].dropna()

def get_fq1(wildcards):
	return {'r1': samples.loc[(wildcards.sample), ["fq1"]].dropna().values[0]}

def get_fq2(wildcards):
	return {'r2': samples.loc[(wildcards.sample), ["fq2"]].dropna().values[0]}

rule trimmomatic_pe:
	input:
		unpack(get_fq1),
		unpack(get_fq2)
	output:
		r1 = temp("outs/{ID}/trimmed/{sample}_1.fq.gz"),
		r2 = temp("outs/{ID}/trimmed/{sample}_2.fq.gz"),
		r1_unpaired = temp("outs/{ID}/trimmed/{sample}_1.unpaired.fq.gz"),
		r2_unpaired = temp("outs/{ID}/trimmed/{sample}_2.unpaired.fq.gz")
	params:
		trimmer = ["TRAILING:3"],
		extra = "",
		compression_level = "-9"
	threads:
		8
	wrapper:
		"0.67.0/bio/trimmomatic/pe"

#rule trimmomatic:
#	input:
#		get_fastq
#	output:
#		temp("outs/{ID}/trimmed/{sample}.fq.gz")
#	params:
#		trimmer = ["TRAILING:3"],
#		extra = "",
#		compression_level = "-9"
#	threads:
#		8
#	wrapper:
#		"0.67.0/bio/trimmomatic/se"
