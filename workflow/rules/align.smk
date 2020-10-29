rule star_index_new:
	input:
		fasta = config["ref"]["fa"],
		gtf = config["ref"]["gtf"],
	threads:
		4
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

rule star_pe_multi:
	input:
		directory("outs/{}/{}".format(config["ID"], config["ref"]["build"])),
		fq1 = [config["fastqs"] + "{sample}_1.fq.gz"],
		fq2 = [config["fastqs"] + "{sample}_2.fq.gz"]
	benchmark:
		"benchmarks/{ID}/align/01_star_align.{sample}.txt"
	log:
		"logs/{ID}/star/{sample}/{sample}.log"
	params:
		index = "outs/{}/{}".format(config["ID"], config["ref"]["build"]),
		extra = "--outSAMtype BAM SortedByCoordinate "
			"--quantMode TranscriptomeSAM GeneCounts"
	output:
	#	"outs/star/{sample}/Aligned.sortedByCoord.out.bam"
		"outs/{ID}/star/{sample}/Aligned.toTranscriptome.out.bam",
		temp("outs/{ID}/star/{sample}/Log.final.out")
	threads:
		4
	wrapper:
		"0.59.1/bio/star/align"
