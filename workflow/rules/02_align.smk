rule star_index_new:
	input: fasta = config["ref"]["fa"], gtf = config["ref"]["gtf"]
	threads: 8
	params: extra = "", build = config["ref"]["build"], analysis = config["ID"], aligner = config["aligner"]
	output:
		directory("{}/{}/STAR/{}".format(outpath, ID, aligner, config["ref"]["build"])),
	benchmark: "benchmarks/align/00_star_index.txt"
	log: "logs/star_index_{}.log".format(config["ref"]["build"])
	wrapper: "0.59.1/bio/star/index"

rule star_pe_multi:
	input:
		directory("{}/{}/{}/{}".format(outpath, ID, aligner, config["ref"]["build"])),
		unpack(fastq_to_aligner)
	benchmark:
		"{outpath}/{ID}/benchmarks/{aligner}/01_star_align.{sample}.txt"
	log: "{outpath}/{ID}/logs/{aligner}/{sample}/{sample}.log"
	params:
		index = "{}/{}/{}".format(outpath, ID, config["ref"]["build"]),
		extra = "--twopassMode Basic --outSAMtype BAM SortedByCoordinate "
		"--quantMode TranscriptomeSAM"
	output:
		"{outpath}/{ID}/{aligner}/{sample}/Aligned.sortedByCoord.out.bam",
		"{outpath}/{ID}/{aligner}/{sample}/Aligned.toTranscriptome.out.bam"
	log: temp("{outpath}/{ID}/{aligner}/{sample}/Log.final.out")
	threads: 8
	wrapper: "0.59.1/bio/star/align"

rule kallisto_index:
	input: fasta = config["ref"]["fa"]
	output: "{}/{}/kallisto/{}.idx".format(outpath, ID, build)
	conda: "../envs/kallisto.yaml"
	shell: "kallisto index -i {output} {input}"

rule kallisto_quant:
	input:
		unpack(fastq_to_aligner),
		#index = "{}/{}/kallisto/{}.idx".format(outpath, ID, build),
		unpack(get_kallisto_index)
	params:
		is_single_end = lambda wildcards: is_single_end(wildcards.sample),
		outdir = "{outpath}/{ID}/kallisto/{sample}"
	output:
		counts_h5 = "{outpath}/{ID}/kallisto/{sample}/abundance.h5",
		counts_tsv = "{outpath}/{ID}/kallisto/{sample}/abundance.tsv",
		log = "{outpath}/{ID}/kallisto/{sample}/run_info.json"
	conda: "../envs/kallisto.yaml"
	log: "{outpath}/{ID}/kallisto/log/{sample}.log"
	shell:
		"is_single_end={params.is_single_end} ; if [[ $is_single_end == False ]]; then "
		"kallisto quant -i {input.index} -o {params.outdir} -b 100 "
		"{input.fq1} {input.fq2} ; "
		"elif [[ $is_single_end == False ]]; then "
		"kallisto quant -i {input.index} -o {output[0]} -b 100 --single "
		"-l 180 -s 20 {input.fq1} ; fi"
