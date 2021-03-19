rule write_se:
	input:
		genes_expected = "{}/{}/RSEM/genes.expected_counts.tsv".format(outpath, ID),
		genes_tpm = "{}/{}/RSEM/genes.tpm_counts.tsv".format(outpath, ID),
		isoforms_expected = "{}/{}/RSEM/isoforms.expected_counts.tsv".format(outpath, ID),
		isoforms_tpm = "{}/{}/RSEM/isoforms.tpm_counts.tsv".format(outpath, ID),
		sample_sheet = config["samples"],
		fastqc = "{}/{}/qc/multiqc_report.{}_data/multiqc_fastqc.txt".format(outpath, ID, ID),
		star = "{}/{}/qc/multiqc_report.{}_data/multiqc_star.txt".format(outpath, ID, ID),
		sortmerna = "{}/{}/qc/multiqc_report.{}_data/multiqc_sortmerna.txt".format(outpath, ID, ID),
		rseqc = "{}/{}/qc/multiqc_report.{}_data/multiqc_rseqc_infer_experiment.txt".format(outpath, ID, ID)
	params:
		organism = config["organism"]
	output:
		genes_SummExp = "{}/{}/SummExp/{}.genes_SummExp.Rds".format(outpath, ID, ID),
		transcripts_SummExp = "{}/{}/SummExp/{}.transcripts_SummExp.Rds".format(outpath, ID, ID)
	conda:
		"../envs/SummExp.yaml"
	script:
		"../scripts/write_SummExp.R"
