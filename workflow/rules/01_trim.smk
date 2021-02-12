rule trimmomatic_pe:
	input:
		unpack(get_fastq)
	output:
		r1 = temp("{outpath}/{ID}/trimmed/{sample}_1.fq.gz"),
		r2 = temp("{outpath}/{ID}/trimmed/{sample}_2.fq.gz"),
		r1_unpaired = temp("{outpath}/{ID}/trimmed/{sample}_1.unpaired.fq.gz"),
		r2_unpaired = temp("{outpath}/{ID}/trimmed/{sample}_2.unpaired.fq.gz")
	params:
		trimmer = ["TRAILING:3"],
		extra = "",
		compression_level = "-9"
	threads:
		8
	wrapper:
		"0.67.0/bio/trimmomatic/pe"

rule trimmomatic:
	input:
		unpack(get_fastq)
	output:
		temp("{outpath}/{ID}/trimmed/{sample}.fq.gz")
	params:
		trimmer = ["TRAILING:3"],
		extra = "",
		compression_level = "-9"
	threads:
		8
	wrapper:
		"0.67.0/bio/trimmomatic/se"
