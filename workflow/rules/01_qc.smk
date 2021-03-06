rule fastqc_fq1:
	input: get_fq1_qc
	output:
		html="{outpath}/{ID}/qc/{sample}_1.html",
		zip=temp("{outpath}/{ID}/qc/{sample}_1_fastqc.zip")
	log: "{outpath}/{ID}/qc/log/{sample}_1.log"
	benchmark: "{outpath}/{ID}/qc/benchmark/{sample}_1.log"
	wrapper:
		"0.31.1/bio/fastqc"

rule fastqc_fq2:
	input: get_fq2_qc
	output:
		html=temp("{outpath}/{ID}/qc/{sample}_2.html"),
		zip=temp("{outpath}/{ID}/qc/{sample}_2_fastqc.zip")
	log: "{outpath}/{ID}/qc/log/{sample}_2.log"
	benchmark: "{outpath}/{ID}/qc/benchmark/{sample}_2.log"
	wrapper: "0.31.1/bio/fastqc"

rule multiqc:
	input:
		expand("{outpath}/{ID}/qc/{sample}_1_fastqc.zip", 
			outpath = outpath, ID = ID, sample = list_of_samples),
		expand("{outpath}/{ID}/qc/{sample}_2_fastqc.zip", 
			outpath = outpath, ID = ID,
			sample = [f for f in list_of_samples if is_paired_end(f)])
	output:
		report = "{}/{}/qc/multiqc_report.{}.html".format(outpath, ID, ID),
		fastqc = temp("{}/{}/qc/multiqc_report.{}_data/multiqc_fastqc.txt".format(outpath, ID, ID))
	log: "{}/{}/qc/log/{}_multiqc.log".format(outpath, ID, ID)
	benchmark: "{}/{}/qc/benchmark/{}_multiqc.log".format(outpath, ID, ID)
	wrapper: "0.51.3/bio/multiqc"
