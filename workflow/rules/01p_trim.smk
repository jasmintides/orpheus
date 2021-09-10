rule trimmomatic_pe:
	input: unpack(get_fastq)
	output:
		r1 = "{outpath}/{ID}/trimmed/{sample}_1.fq.gz",
		r2 = "{outpath}/{ID}/trimmed/{sample}_2.fq.gz",
		r1_unpaired = temp("{outpath}/{ID}/trimmed/{sample}_1.unpaired.fq.gz"),
		r2_unpaired = temp("{outpath}/{ID}/trimmed/{sample}_2.unpaired.fq.gz")
	params: trimmer = ["TRAILING:3"], extra = "", compression_level = "-9"
	log: "{outpath}/{ID}/trimmed/log/{sample}.txt"
	benchmark: "{outpath}/{ID}/trimmed/benchmark/{sample}.txt"
	threads: 32
	wrapper: "0.67.0/bio/trimmomatic/pe"

rule trimmomatic:
	input: lambda wildcards: get_fastq(wildcards.sample)
	output: temp("{outpath}/{ID}/trimmed/{sample}.fq.gz")
	params: trimmer = ["TRAILING:3"], extra = "", compression_level = "-9"
	benchmark: "{outpath}/{ID}/trimmed/log/{sample}.txt"
	threads: 32
	wrapper: "0.67.0/bio/trimmomatic/se"
