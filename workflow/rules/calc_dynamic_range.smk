rule mark_duplicates_again:
	input:
		"{outpath}/{ID}/star/{sample}/Aligned.sortedByCoord.out.bam"
	log:
		"{outpath}/{ID}/logs/picard/dedup/{sample}.log"
	output:
		bam = "{outpath}/{ID}/picard/{sample}.dupMarked.bam",
		metrics = "{outpath}/{ID}/picard/{sample}.dupMarked.metrics.txt"
	params:
		mem = "-Xmx8g"
	wrapper:
		"0.57.0/bio/picard/markduplicates"

rule samtools_index:
	input:
		"{outpath}/{ID}/picard/{sample}.dupMarked.bam",
	log:
		"{outpath}/{ID}/logs/samtools/{sample}.index.log"
	output:
		temp("{outpath}/{ID}/picard/{sample}.dupMarked.bam.bai")
	params:
		""
	wrapper:
		"0.72.0/bio/samtools/index"

rule gtf2bed:
	input:
		gtf = config["ref"]["gtf"]
	log:
		"{}/logs/{}/bedops/{}.log".format(outpath, ID, build)
	output:
		bed = "{}/{}/bedops/{}.gtf2bed.bed".format(outpath, ID, build)
	conda:
		"../envs/bedops.yaml"
	shell:
		"gtf2bed < {input.gtf} > {output.bed}"

rule infer_experiment:
	input:
		bam = "{outpath}/{ID}/picard/{sample}.dupMarked.bam",
		bai = "{outpath}/{ID}/picard/{sample}.dupMarked.bam.bai",
		bed = "{}/{}/bedops/{}.gtf2bed.bed".format(outpath, ID, build)
	output:
		"{outpath}/{ID}/rseqc/{sample}.txt"
	log:
		"{outpath}/logs/{ID}/rseqc/{sample}.log"
	conda:
		"../envs/RSeQC.yaml"
	shell:
		"infer_experiment.py -r {input.bed} -i {input.bam} > {output}"

rule calc_dynamic_range:
	input:
		bam = "{outpath}/{ID}/picard/{sample}.dupMarked.bam",
		strand_log = "{outpath}/{ID}/rseqc/{sample}.txt",
		gtf = config["ref"]["gtf"]
	params:
		is_paired_end = lambda wildcards: is_paired_end(wildcards.sample)
	output:
		"{outpath}/{ID}/rseqc/{sample}.dynamic_range.txt"
	threads:
		4
	conda:
		"../envs/dupRadar.yaml"
	script:
		"../scripts/calculate_dynamic_range.R"
