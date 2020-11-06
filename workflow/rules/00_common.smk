def is_single_end(sample):
	return pd.isnull(samples.loc[(sample), "fq2"])
