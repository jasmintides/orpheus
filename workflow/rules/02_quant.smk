rule star_index_new:
	input: 
		fasta = config["ref"]["fa"], 
		gtf = config["ref"]["gtf"]
	params:
		extra = ""
	output: temp(directory("{}/{}/STAR/{}".format(outpath, ID, build)))
	log: "{}/{}/STAR/log/{}_index.txt".format(outpath, ID, build)
	benchmark: "{}/{}/STAR/benchmark/{}_index.txt".format(outpath, ID, build)
	wrapper: "0.59.1/bio/star/index"

rule star_pe_multi:
	input:
		lambda wildcards: get_star_index(wildcards),
		unpack(fastq_to_aligner)
	params:
		index = "{}/{}/{}".format(outpath, ID, config["ref"]["build"]),
		extra = "--twopassMode Basic --quantMode TranscriptomeSAM"
	output:
		temp("{outpath}/{ID}/STAR/{sample}/Aligned.toTranscriptome.out.bam")
	log: "{outpath}/{ID}/STAR/log/{sample}.txt"
	benchmark: "{outpath}/{ID}/STAR/benchmark/{sample}.txt"
	threads: 8
	wrapper: "0.59.1/bio/star/align"

rule rsem_prepare_reference:
	input: fasta = config["ref"]["fa"], gtf = config["ref"]["gtf"]
	params: outpath = config["outpath"], ID = config["ID"], 
		build = config["ref"]["build"]
	output: temp(directory("{}/{}/RSEM/{}".format(outpath, ID, build)))
	log: "{}/{}/RSEM/log/{}_index.txt".format(outpath, ID, build)
	benchmark: "{}/{}/RSEM/benchmark/{}_index.txt".format(outpath, ID, build)
	conda: "../envs/RSEM.yaml"
	shell:
		"rm -rf {params.outpath}/{params.ID}/RSEM/{params.build} ; "
		"mkdir {params.outpath}/{params.ID}/RSEM/{params.build} ; "
		"rsem-prepare-reference --num-threads {threads} --gtf {input.gtf} "
		"{input.fasta} {params.outpath}/{params.ID}/RSEM/{params.build}"

rule rsem_calculate_expression:
	input:
		bam = "{outpath}/{ID}/STAR/{sample}/Aligned.toTranscriptome.out.bam",	
		ref = directory("{}/{}/RSEM/{}".format(outpath, ID, build))
	params:
		is_single_end = lambda wildcards: is_single_end(wildcards.sample),
		build = config["ref"]["build"]
	output: temp("{outpath}/{ID}/RSEM/{sample}.isoforms.results")
	log: "{outpath}/{ID}/RSEM/log/{sample}.txt"
	benchmark: "{outpath}/{ID}/RSEM/benchmark/{sample}.txt"
	conda:
		"../envs/RSEM.yaml"
	shell:
		"is_single_end={params.is_single_end} ; if [[ $is_single_end == False ]]; then "
		"rsem-calculate-expression --num-threads {threads} "
		"--fragment-length-max 1000 --no-bam-output --paired-end "
		"--bam {input.bam} {input.ref}/{params.build} "
		"{wildcards.outpath}/{wildcards.ID}/RSEM/{wildcards.sample} ; "
		"elif [[ $is_single_end == True ]]; then "
		"rsem-calculate-expression --num-threads {threads} "
		"--fragment-length-max 1000 --no-bam-output "
		"--bam {input.bam} {input.ref}/{params.build} "
		"{wildcards.outpath}/{wildcards.ID}/RSEM/{wildcards.sample} ; fi"

rule kallisto_index:
	input: fasta = config["ref"]["fa"]
	output: temp("{}/{}/kallisto/{}.idx".format(outpath, ID, build))
	conda: "../envs/kallisto.yaml"
	log: "{}/{}/kallisto/log/{}_index.log".format(outpath, ID, build)
	benchmark: "{}/{}/kallisto/benchmark/{}_index.log".format(outpath, ID, build)
	shell: "kallisto index -i {output} {input}"

rule kallisto_quant:
	input:
		unpack(fastq_to_aligner),
		unpack(get_kallisto_index)
	params:
		is_single_end = lambda wildcards: is_single_end(wildcards.sample),
		outdir = "{outpath}/{ID}/kallisto/{sample}",
		gtf = config["ref"]["gtf"]
	output:
		counts_h5 = temp("{outpath}/{ID}/kallisto/{sample}/abundance.h5"),
		counts_tsv = temp("{outpath}/{ID}/kallisto/{sample}/abundance.tsv"),
		log = temp("{outpath}/{ID}/kallisto/{sample}/run_info.json")
	conda: "../envs/kallisto.yaml"
	log: "{outpath}/{ID}/kallisto/log/{sample}.log"
	benchmark: "{outpath}/{ID}/kallisto/benchmark/{sample}.log"
	resources: mem_mb = 4000
	shell:
		"is_single_end={params.is_single_end} ; if [[ $is_single_end == False ]]; then "
		"kallisto quant -i {input.index} -o {params.outdir} -b 100 "
		"-g {params.gtf} {input[0]} {input[1]} ; "
		"elif [[ $is_single_end == True ]]; then "
		"kallisto quant -i {input.index} -o {params.outdir} -b 100 --single "
		"-l 180 -s 20 -g {params.gtf} {input[0]} ; fi"
