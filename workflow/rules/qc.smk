rule fastqc:
	input:
		config['fastqs'] + '{sample}_{pair}.fq.gz' 
	output:
		html="outs/{ID}/qc/{sample}_{pair}.html",
		zip=temp("outs/{ID}/qc/{sample}_{pair}_fastqc.zip")
	wrapper:
		"0.31.1/bio/fastqc"

rule multiqc:
	input:
		expand("outs/{ID}/qc/{sample}_{pair}_fastqc.zip", ID = ID, sample = samples, pair = pairs),
		expand("outs/{ID}/star/{sample}/Log.final.out", ID = ID, sample = samples, pair = pairs)
	output:
		"outs/{}/qc/multiqc_report.{}.html".format(config["ID"], config["ID"])
	wrapper:
		"0.51.3/bio/multiqc"
