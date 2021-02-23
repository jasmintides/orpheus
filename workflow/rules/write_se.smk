rule write_se:
	input:
		genes_expected = "{}/{}/RSEM/genes.expected_counts.tsv".format(outpath, ID),
		genes_tpm = "{}/{}/RSEM/genes.tpm_counts.tsv".format(outpath, ID).
		isoforms_expected = "{}/{}/RSEM/isoforms.expected_counts.tsv".format(outpath, ID),
		isoforms_tpm = "{}/{}/RSEM/isoforms.tpm_counts.tsv".format(outpath, ID)
	output:
		"{}/{}/RSEM/{}.dynamic_range.txt".format(outpath, ID, ID)
	conda:
		"../envs/summarizedexperiment.yaml"
	script:
		"../scripts/write_SummarizedExperiment.R"
