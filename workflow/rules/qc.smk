rule fastqc_fq1:
	input:
		get_fq1_qc
	output:
		html="{outpath}/{ID}/qc/{sample}_1.html",
		zip=temp("{outpath}/{ID}/qc/{sample}_1_fastqc.zip")
	wrapper:
		"0.31.1/bio/fastqc"

rule fastqc_fq2:
	input:
		get_fq2_qc
	output:
		html=temp("{outpath}/{ID}/qc/{sample}_2.html"),
		zip=temp("{outpath}/{ID}/qc/{sample}_2_fastqc.zip")
	wrapper:
		"0.31.1/bio/fastqc"

rule multiqc:
	input:
		expand("{outpath}/{ID}/qc/{sample}_1_fastqc.zip", 
			outpath = outpath, ID = ID, sample = list_of_samples),
		expand("{outpath}/{ID}/qc/{sample}_2_fastqc.zip", 
			outpath = outpath, ID = ID,
			sample = [f for f in list_of_samples if is_paired_end(f)]),
		expand("{outpath}/{ID}/star/{sample}/Log.final.out", 
			outpath = outpath, ID = ID, sample = list_of_samples),
		expand("{outpath}/{ID}/sortmerna/{sample}.rRNA.log", 
			outpath = outpath, ID = ID, sample = list_of_samples),
		expand("{outpath}/{ID}/rseqc/{sample}.txt",
			outpath = outpath, ID = ID, sample = list_of_samples)
	output:
		report = "{}/{}/qc/multiqc_report.{}.html".format(outpath, ID, ID),
		fastqc = "{}/{}/qc/multiqc_report.{}_data/multiqc_fastqc.txt".format(outpath, ID, ID),
		star = "{}/{}/qc/multiqc_report.{}_data/multiqc_star.txt".format(outpath, ID, ID),
		sortmerna = "{}/{}/qc/multiqc_report.{}_data/multiqc_sortmerna.txt".format(outpath, ID, ID),
		rseqc = "{}/{}/qc/multiqc_report.{}_data/multiqc_rseqc_infer_experiment.txt".format(outpath, ID, ID)
	wrapper:
		"0.51.3/bio/multiqc"
