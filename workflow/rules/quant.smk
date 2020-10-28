rule rsem_prepare_reference:
	input:
		fasta = config["ref"]["fa"],
		gtf = config["ref"]["gtf"]
	params:
		ID = config["ID"],
		build = config["ref"]["build"]
	output:
		"outs/{}/ref/{}.transcripts.fa".format(config["ID"], config["ref"]["build"])
#	benchmark:
#		"benchmarks/quant/00_rsem_prepare_reference.txt"
#	log:
#		"logs/rsem/reference.log"
	threads:
		4
	conda:
		"../envs/quant.yaml"
	shell:
		"rsem-prepare-reference --num-threads {threads} --gtf {input.gtf} "
		"{input.fasta} outs/{params.ID}/ref/{params.build}"

rule rsem_calculate_expression:
	input:
		bam = "outs/{ID}/star/{sample}/Aligned.toTranscriptome.out.bam"
#		ref = directory("outs/{}/RSEM/{}".format(config["ID"], config["ref"]["build"]))
	params:
		ID = config["ID"],
		build = config["ref"]["build"]
	output:
		"outs/{ID}/RSEM/{sample}.isoforms.results",
		"outs/{ID}/RSEM/{sample}.genes.results"
	benchmark:
		"benchmarks/{ID}/quant/01_rsem.{sample}.txt"
	log:
		"logs/{ID}/rsem/01_rsem.{sample}.log"
	threads:
		4
	conda:
		"../envs/quant.yaml"
	shell:
		"rsem-calculate-expression --num-threads {threads} "
		"--fragment-length-max 1000 --no-bam-output --paired-end --estimate-rspd --forward-prob 0.0 "
		"--bam {input.bam} outs/{params.ID}/ref/{params.build} "
		"outs/{ID}/RSEM/{wildcards.sample}"

rule star_raw_counts:
	input:
		counts = expand("outs/{ID}/star/{sample}/ReadsPerGene.out.tab")

