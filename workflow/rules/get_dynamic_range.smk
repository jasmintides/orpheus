def get_fastq(wildcards):
        return samples.loc[(wildcards.sample), ["fq1", "fq2"]].dropna()

def get_fq1(wildcards):
        return {'r1': samples.loc[(wildcards.sample), ["fq1"]].dropna().values[0]}

def get_fq2(wildcards):
        return {'r2': samples.loc[(wildcards.sample), ["fq2"]].dropna().values[0]}

def get_trimmed(wildcards):
	if not is_single_end(wildcards.sample):
		return {'fq1': expand("outs/{ID}/trimmed/{sample}_{group}.fq.gz", group = 1, **wildcards),
			'fq2': expand("outs/{ID}/trimmed/{sample}_{group}.fq.gz", group = 2, **wildcards)}
	else:
		return "outs/{ID}/{sample}.fq.gz".format(**wildcards)

rule merge_paired_reads:
	input:
		unpack(get_fq1),
		unpack(get_fq2)
	output:
		"outs/{ID}/sortmerna/{sample}/out/other.fasta"
	params:
		workdir = "/data/exploratory/Users/jeff.alvarez/orpheus"
	conda:
		"../envs/sortmerna.yaml"
	log:
		"logs/{ID}/sortmerna/{sample}.log"
	threads:
		8
	shell:
		"sortmerna -ref {params.workdir}/resources/smr_v4.3_default_db.fasta "
		"-reads {params.workdir}/{input.r1} -reads {params.workdir}/{input.r2} "
		"-workdir {params.workdir}/outs/{ID}/sortmerna/{wildcards.sample}/ "
		"-fastx -paired_in -v -threads {threads}"
