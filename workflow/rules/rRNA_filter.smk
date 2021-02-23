rule unzip_paired_reads:
	input:
		unpack(get_fastq)
	output:
		r1 = temp("{outpath}/{ID}/sortmerna/paired/{sample}_1.fastq"),
		r2 = temp("{outpath}/{ID}/sortmerna/paired/{sample}_2.fastq")
	run:
		if input.r1.endswith(".gz") and is_paired_end(wildcards.sample):
			shell("gunzip -cd {input.r1} > {output.r1}")
		else:
			shell("ln -sr {input.r1} {output.r1}")
		if input.r2.endswith(".gz") and is_paired_end(wildcards.sample):
			shell("gunzip -cd {input.r2} > {output.r2}")
		else:
			shell("ln -sr {input.r2} {output.r2}")

rule merge_paired_reads:
	input:
		r1 = "{outpath}/{ID}/sortmerna/paired/{sample}_1.fastq",
		r2 = "{outpath}/{ID}/sortmerna/paired/{sample}_2.fastq"
	output:
		merged = "{outpath}/{ID}/sortmerna/merged/{sample}.fastq"
	conda:
		"../envs/sortmerna.yaml"
	shell:
		"bash merge-paired-reads.sh {input.r1} {input.r2} {output.merged}"
	
rule unzip_reads:
	input:
		unpack(get_fastq)
	output:
		merged = "{outpath}/{ID}/sortmerna/single_end/{sample}.fastq"
	run:
		if input.r1.endswith(".gz") and not is_paired_end(wildcards.sample):
			shell("gunzip -cd {input.r1} > {output.merged}")
		else:
			shell("ln -sr {input.r1} {output.merged}")

def get_merged_reads(wildcards):
	if not is_single_end(wildcards.sample):
		# paired-end
		return{'fastq': expand("{outpath}/{ID}/sortmerna/merged/{sample}.fastq", **wildcards)}
	else:
		# single-end
		return{'fastq': expand("{outpath}/{ID}/sortmerna/single_end/{sample}.fastq", **wildcards)}

rule sortmerna:
	input:
		unpack(get_merged_reads)
	params:
		is_single_end = lambda wildcards: is_single_end(wildcards.sample),
		rRNA_db_path = "/data/exploratory/genome/sortmerna"
	output:
		log = "{outpath}/{ID}/sortmerna/{sample}.rRNA.log",
		clean_fq = temp("{outpath}/{ID}/sortmerna/{sample}.clean.fastq"),
		rRNA_fq = temp("{outpath}/{ID}/sortmerna/{sample}.rRNA.fastq")
	conda:
		"../envs/sortmerna.yaml"
	shell:
		"is_single_end={params.is_single_end} ; if [[ $is_single_end == False ]]; then "
		"sortmerna --ref {params.rRNA_db_path}/silva-arc-16s-id95.fasta,"
		"{params.rRNA_db_path}/silva-arc-16s-db --reads {input.fastq} "
		"--aligned {wildcards.outpath}/{wildcards.ID}/sortmerna/{wildcards.sample}.rRNA "
		"--fastx --paired_in --other {outpath}/{ID}/sortmerna/{wildcards.sample}.clean "
		"--log -v ; elif [[ $is_single_end == True ]]; then "
		"sortmerna --ref {params.rRNA_db_path}/silva-arc-16s-id95.fasta,"
		"{params.rRNA_db_path}/silva-arc-16s-db --reads {input.fastq} "
		"--aligned {wildcards.outpath}/{wildcards.ID}/sortmerna/{wildcards.sample}.rRNA --fastx "
		"--other {wildcards.outpath}/{wildcards.ID}/sortmerna/{wildcards.sample}.clean --log -v ; fi"
