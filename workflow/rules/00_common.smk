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
