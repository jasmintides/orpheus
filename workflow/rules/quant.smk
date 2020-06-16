rule rsem_prepare_reference:
	input:
		fasta = config["ref"]["fa"],
		gtf = config["ref"]["gtf"],
		build = config["ref"]["build"],
		ID = config["ID"]
	output:
		directory("outs/{}/RSEM/{}".format(config["ID"], config["ref"]["build"]))
	benchmark:
		"benchmarks/quant/00_rsem_prepare_reference.txt"
	log:
		"logs/rsem/reference.log"
	threads:
		12
	conda:
		"workflow/envs/quant.yaml"
	shell:
		"rsem-prepare-reference --num-threads {threads} --gtf {input.gtf} "
		"{input.fasta} outs/{input.ID}/RSEM/{input.build}"

rule rsem_calculate_expression:
	input:
		bam = "outs/star/{sample}/Aligned.sortedByCoord.out.bam",
		ref = directory("outs/{}/RSEM/{}".format(config["ID"], config["ref"]["build"]))
	output:
		directory("outs/RSEM/{sample}")
	benchmark:
		"benchmarks/quant/01_rsem.{sample}.txt"
	log:
		"logs/rsem/01_rsem.{sample}.log"
	threads:
		12
	conda:
		"workflow/envs/quant.yaml"
	shell:
		"rsem-calculate-expression --num-threads {threads} "
		"--fragment-length-max 1000 --no-bam-output --paired-end --estimate-rspd --forward-prob 0.0 "
		"--bam {input.bam} {input.ref} outs/RSEM/{sample}"
