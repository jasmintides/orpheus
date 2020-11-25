rule star_index_new:
	input:
		fasta = config["ref"]["fa"],
		gtf = config["ref"]["gtf"],
	threads:
		8
	params:
		extra = "",
		build = config["ref"]["build"],
		analysis = config["ID"]
	output:
		directory("outs/{}/{}".format(config["ID"], config["ref"]["build"])),
	benchmark:
		"benchmarks/align/00_star_index.txt"
	log:
		"logs/star_index_{}.log".format(config["ref"]["build"])
	wrapper:
		"0.59.1/bio/star/index"

def get_trimmed(wildcards):
	if not is_single_end(wildcards.sample):
		### paired-end sample ###
		#checkpoint_f1 = checkpoints.trimmomatic_pe.get(**wildcards).output.r1
		#checkpoint_f2 = checkpoints.trimmomatic_pe.get(**wildcards).output.r2
		#return {'fq1': checkpoint_f1[0],
		#	'fq2': checkpoint_f2[0]}
		return {'fq1': expand("outs/{ID}/trimmed/{sample}_{group}.fq.gz", group = 1, **wildcards),
			'fq2': expand("outs/{ID}/trimmed/{sample}_{group}.fq.gz", group = 2, **wildcards)}
	# single end sample
	else:
		return "outs/{ID}/{sample}.fq.gz".format(**wildcards)	

rule star_pe_multi:
#checkpoint star_pe_multi:
	input:
		directory("outs/{}/{}".format(config["ID"], config["ref"]["build"])),
		unpack(get_trimmed)
	benchmark:
		"benchmarks/{ID}/align/01_star_align.{sample}.txt"
	log:
		"logs/{ID}/star/{sample}/{sample}.log"
	params:
		index = "outs/{}/{}".format(config["ID"], config["ref"]["build"]),
		extra = "--twopassMode Basic --outSAMtype BAM SortedByCoordinate "
			"--quantMode TranscriptomeSAM GeneCounts"
	output:
		temp("outs/{ID}/star/{sample}/ReadsPerGene.out.tab"),
		temp("outs/{ID}/star/{sample}/Aligned.sortedByCoord.out.bam"),
		temp("outs/{ID}/star/{sample}/Aligned.toTranscriptome.out.bam"),
		log = temp("outs/{ID}/star/{sample}/Log.final.out")
	threads:
		8
	wrapper:
		"0.59.1/bio/star/align"

def get_template(wildcards):
	checkpoint_output = checkpoints.star_pe_multi.get(ID = wildcards.ID, sample = list_of_samples[0]).output[0]
	#template = checkpoint_output[0]
	return {'template': checkpoint_output[0]}

rule create_gene_ids_star:
	input:
#		unpack(get_template)
		expand("outs/{ID}/star/{sample}/ReadsPerGene.out.tab", ID = config["ID"], sample = list_of_samples[0])
	params:
		ID = config["ID"]
	output:
		#temp("outs/{}/star/gene_ids.txt".format(config["ID"]))
		temp("outs/{ID}/star/gene_ids.txt")
	shell:
		"tail -n +5 {input} | cut -f1 | sed '1i \n' > "
		"outs/{params.ID}/star/gene_ids.txt"
