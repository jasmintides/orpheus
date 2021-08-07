rule star_index_new:
	input: fasta = config["ref"]["fa"], gtf = config["ref"]["gtf"]
	threads: 8
	params: extra = "", build = config["ref"]["build"], analysis = config["ID"]
	output:
		directory("{}/{}/STAR/{}".format(outpath, ID, config["ref"]["build"])),
	benchmark: "benchmarks/align/00_star_index.txt"
	log: "logs/star_index_{}.log".format(config["ref"]["build"])
	wrapper: "0.59.1/bio/star/index"

rule star_pe_multi:
	input:
		directory("{}/{}/STAR/{}".format(outpath, ID, config["ref"]["build"])),
		unpack(get_trimmed)
	benchmark:
		"{outpath}/{ID}/benchmarks/STAR/01_star_align.{sample}.txt"
	log: "{outpath}/{ID}/logs/STAR/{sample}/{sample}.log"
	params:
		index = "{}/{}/{}".format(outpath, ID, config["ref"]["build"]),
		extra = "--twopassMode Basic --outSAMtype BAM SortedByCoordinate "
			"--quantMode TranscriptomeSAM"
	output:
		"{outpath}/{ID}/STAR/{sample}/Aligned.sortedByCoord.out.bam",
		"{outpath}/{ID}/STAR/{sample}/Aligned.toTranscriptome.out.bam",
	log: temp("{outpath}/{ID}/STAR/{sample}/Log.final.out")
	threads: 8
	wrapper: "0.59.1/bio/star/align"

rule rsem_prepare_reference:
	input: fasta = config["ref"]["fa"], gtf = config["ref"]["gtf"]
	params: outpath = config["outpath"], ID = config["ID"], 
			build = config["ref"]["build"]
	output: directory("{}/{}/RSEM/{}".format(outpath, ID, build))
	conda: "../envs/quant.yaml"
	threads: 4
	shell: "rm -rf {params.outpath}/{params.ID}/{params.build} && "
	"mkdir {params.outpath}/{params.ID}/{params.build} && "
	"rsem-prepare-reference --num-threads {threads} --gtf {input.gtf} "
	"{input.fasta} {params.outpath}/{params.ID}/{params.build}"

rule rsem_calculate_expression:
	input:
		bam = "{outpath}/{ID}/STAR/{sample}/Aligned.toTranscriptome.out.bam",	
		ref = directory("{}/{}/RSEM/{}".format(outpath, ID, build))
	params:
		is_single_end = lambda wildcards: is_single_end(wildcards.sample),
		build = config["ref"]["build"]
	output:
		temp("{outpath}/{ID}/RSEM/{sample}.isoforms.results"),
		temp("{outpath}/{ID}/RSEM/{sample}.genes.results")
	threads: 4
	conda: "../envs/quant.yaml"
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
