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
		# paired-end sample
		return {'fq1': expand("outs/{ID}/trimmed/{sample}_{group}.fq.gz", group = 1, **wildcards),
			'fq2': expand("outs/{ID}/trimmed/{sample}_{group}.fq.gz", group = 2, **wildcards)}
	# single end sample
	else:
		return "outs/{ID}/{sample}.fq.gz".format(**wildcards)	

rule star_pe_multi:
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
		temp("outs/{ID}/star/{sample}/Aligned.sortedByCoord.out.bam"),
		temp("outs/{ID}/star/{sample}/Aligned.toTranscriptome.out.bam"),
		log = temp("outs/{ID}/star/{sample}/Log.final.out")
	threads:
		8
	wrapper:
		"0.59.1/bio/star/align"
