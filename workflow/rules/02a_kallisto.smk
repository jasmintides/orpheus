rule kallisto_index:
	input: fasta = config["ref"]["fa"]
	output: "{}/{}/kallisto/{}".format(outpath, ID, build)
	conda: "../envs/kallisto.yaml"
	shell: "kallisto index -i {output} {input}"

rule kallisto_quant:
	input:
		unpack(get_trimmed),
		index = "{}/{}/kallisto/{}".format(outpath, ID, build)
	params:
		is_single_end = lambda wildcards: is_single_end(wildcards.sample),
	output:
		outdir = directory("{outpath}/{ID}/kallisto/{sample}"),
		counts_h5 = "{outpath}/{ID}/kallisto/{sample}/abundance.h5",
		counts_tsv = "{outpath}/{ID}/kallisto/{sample}/abundance.tsv",
		log = "{outpath}/{ID}/kallisto/{sample}/run_info.json"
	conda: "../envs/kallisto.yaml"
	log: "{outpath}/{ID}/kallisto/log/{sample}.log"
	shell:
		"is_single_end={params.is_single_end} ; if [[ $is_single_end == False ]]; then "
		"kallisto quant -i {input.index} -o {output[0]} -b 100 "
		"{input.fq1} {input.fq2} ; "
		"elif [[ $is_single_end == False ]]; then "
		"kallisto quant -i {input.index} -o {output[0]} -b 100 --single "
		"-l 180 -s 20 {input.fq1} ; fi"
