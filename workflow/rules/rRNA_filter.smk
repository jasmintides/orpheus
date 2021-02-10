def get_fastq(wildcards):
        return samples.loc[(wildcards.sample), ["fq1", "fq2"]].dropna()

def get_fq1(wildcards):
        return {'r1': samples.loc[(wildcards.sample), ["fq1"]].dropna().values[0]}

def get_fq2(wildcards):
        return {'r2': samples.loc[(wildcards.sample), ["fq2"]].dropna().values[0]}

rule merge_paired_reads:
	input:
		unpack(get_fq1),
		unpack(get_fq2)
	output:
		"outs/{ID}/sortmerna/{sample}.fastq"
	conda:
		"../envs/sortmerna.yaml"
	log:
		"logs/{ID}/sortmerna/{sample}.log"
	shell:
#		"if [[ {input.r1} =~ \.gz ]] ; then gunzip -cd {input.r1} > 
		"bash merge-paired-reads.sh {input.r1} {input.r2} outs/{ID}/sortmerna/{wildcards.sample}.fastq"

rule sortmerna:
	input:
		fastq = "outs/{ID}/sortmerna/{sample}.fastq"
	params:
		workdir = "/data/exploratory/Users/jeff.alvarez/orpheus"
	output:
		clean = "outs/{ID}/sortmerna/{sample}.clean.fastq"
	conda:
		"../envs/sortmerna.yaml"
	shell:
		"sortmerna --ref resources/silva-arc-16s-id95.fasta,resources/silva-arc-16s-db "
		"--reads {input.fastq} --aligned outs/{ID}/sortmerna/{wildcards.sample}.rRNA "
		"--paired_in --fastx --other outs/{ID}/sortmerna/{wildcards.sample}.clean --log -v" 
