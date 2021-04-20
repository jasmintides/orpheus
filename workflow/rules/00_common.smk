def get_samples(sample_sheet):
	if config["run_sra_pipeline"] is True:
		samples = pd.read_table(sample_sheet, dtype = str).set_index("sample", drop = False)
		list_of_samples = samples["sample"].tolist()
		template_sample = list_of_samples[0]
		return {'samples': list_of_samples,
			'template_sample': template_sample}
	else:
		samples = pd.read_table(config["samples"], dtype = str).set_index("sample", drop = False)
		list_of_samples = samples["sample"].tolist()
		template_sample = list_of_samples[0]
		return {'samples': list_of_samples,
			'template_sample': template_sample}

def is_single_end(sample):
	return pd.isnull(samples.loc[(sample), "fq2"])

def is_paired_end(sample):
	return pd.notnull(samples.loc[(sample), "fq2"])

def get_fastq(wildcards):
	if not is_single_end(wildcards.sample):
		# paired-end
		return {'r1': samples.loc[(wildcards.sample), ["fq1"]].dropna().values[0],
			'r2': samples.loc[(wildcards.sample), ["fq2"]].dropna().values[0]}
	else:
		# single-end
		return {'r1': samples.loc[(wildcards.sample), ["fq1"]].dropna().values[0]}
