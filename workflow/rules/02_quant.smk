rule star_index_new:
	input: 
		fasta = config["ref"]["fa"], 
		gtf = config["ref"]["gtf"]
	threads: 8
	params:
		extra = ""
	output: directory("{}/{}/STAR/{}".format(outpath, ID, build))
	benchmark: "benchmarks/align/00_star_index.txt"
	log: "logs/star_index_{}.log".format(build)
	wrapper: "0.59.1/bio/star/index"

rule star_pe_multi:
	input:
		lambda wildcards: get_star_index(wildcards),
		unpack(fastq_to_aligner)
	params:
		index = "{}/{}/{}".format(outpath, ID, config["ref"]["build"]),
		extra = "--twopassMode Basic --outSAMtype BAM SortedByCoordinate "
		"--quantMode TranscriptomeSAM"
	output:
		"{outpath}/{ID}/STAR/{sample}/Aligned.sortedByCoord.out.bam",
		temp("{outpath}/{ID}/STAR/{sample}/Aligned.toTranscriptome.out.bam")
	threads: 8
	wrapper: "0.59.1/bio/star/align"

rule rsem_prepare_reference:
	input: fasta = config["ref"]["fa"], gtf = config["ref"]["gtf"]
	params: outpath = config["outpath"], ID = config["ID"], 
		build = config["ref"]["build"]
	output: directory("{}/{}/RSEM/{}".format(outpath, ID, build))
	conda: "../envs/quant.yaml"
	threads: 4
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
	output:
		temp("{outpath}/{ID}/RSEM/{sample}.isoforms.results")
	threads:
		4
	conda:
		"../envs/quant.yaml"
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
	output: "{}/{}/kallisto/{}.idx".format(outpath, ID, build)
	conda: "../envs/kallisto.yaml"
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
		counts_h5 = "{outpath}/{ID}/kallisto/{sample}/abundance.h5",
		counts_tsv = "{outpath}/{ID}/kallisto/{sample}/abundance.tsv",
		log = "{outpath}/{ID}/kallisto/{sample}/run_info.json"
	conda: "../envs/kallisto.yaml"
	log: "{outpath}/{ID}/kallisto/log/{sample}.log"
	benchmark: "{outpath}/{ID}/kallisto/benchmark/{sample}.log"
	threads: 16
	shell:
		"is_single_end={params.is_single_end} ; if [[ $is_single_end == False ]]; then "
		"kallisto quant -i {input.index} -o {params.outdir} -b 100 "
		"-g {params.gtf} {input.fq1} {input.fq2} ; "
		"elif [[ $is_single_end == False ]]; then "
		"kallisto quant -i {input.index} -o {output[0]} -b 100 --single "
		"-l 180 -s 20 -g {params.gtf} {input.fq1} ; fi"
