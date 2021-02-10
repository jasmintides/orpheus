def is_single_end(sample):
	return pd.isnull(samples.loc[(sample), "fq2"])

def is_paired_end(sample):
	return pd.notnull(samples.loc[(sample), "fq2"])
