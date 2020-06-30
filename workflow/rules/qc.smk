rule fastqc:
	input:
		config['fastqs'] + '/{sample}_{pair}.fq.gz' 
	output:
		html="outs/qc/{sample}_{pair}.html",
		zip="outs/qc/{sample}_{pair}_fastqc.zip"
	wrapper:
		"0.31.1/bio/fastqc"

rule multiqc:
	input:
		expand("outs/qc/{sample}_{pair}_fastqc.zip", sample = samples, pair = pairs),
		expand("outs/star/{sample}/Log.final.out", sample = samples, pair = pairs)
	output:
		"outs/{}/multiqc_report.html".format(config["ID"])
	wrapper:
		"0.51.3/bio/multiqc"
