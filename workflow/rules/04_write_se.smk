rule write_se:
	input:
		**get_full_counts(aligner),
		sample_sheet = config["samples"],
		fastqc = "{}/{}/qc/multiqc_report.{}_data/multiqc_fastqc.txt".format(outpath, ID, ID)
	output:
		"{outpath}/{ID}/SummExp/{ID}.SummExp.Rds"
	conda:
		"../envs/SummExp.yaml"
	script:
		"../scripts/write_SummExp.R"
