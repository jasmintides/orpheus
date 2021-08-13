def is_single_end(sample):
	return pd.isnull(samples.loc[(sample), "fq2"])

def is_paired_end(sample):
	return pd.notnull(samples.loc[(sample), "fq2"])

def get_fastq(wildcards):
	if not is_single_end(wildcards.sample):
		return {'r1': samples.loc[(wildcards.sample), ["fq1"]].dropna().values[0],
			'r2': samples.loc[(wildcards.sample), ["fq2"]].dropna().values[0]}
	else:
		return {'r1': samples.loc[(wildcards.sample), ["fq1"]].dropna().values[0]}

def fastq_to_aligner(wildcards):
	if config["skip_trimming"] is True:
		if not is_single_end(wildcards.sample):
			return {'fq1': samples.loc[(wildcards.sample), ["fq1"]].dropna().values[0],
				'fq2': samples.loc[(wildcards.sample), ["fq2"]].dropna().values[0]}
		else:
			return {'fq1': samples.loc[(wildcards.sample), ["fq1"]].dropna().values[0]}

	elif config["skip_trimming"] is False:
		if not is_single_end(wildcards.sample):
			return {'fq1': expand("{outpath}/{ID}/trimmed/{sample}_{group}.fq.gz", 
				group = 1, **wildcards),
				'fq2': expand("{outpath}/{ID}/trimmed/{sample}_{group}.fq.gz", 
					group = 2, **wildcards)}
		else:
			return {'fq1': "{outpath}/{ID}/trimmed/{sample}.fq.gz".format(**wildcards)}

def get_transcript_counts(aligner):
	transcript_ids = list()
	if aligner in ["Kallisto", "KALLISTO", "kallisto"]:
		string = "{outpath}/{ID}/kallisto/{sample}/abundance.tsv"
		transcript_ids.append(string)
	elif config["aligner"] in ["STAR", "Star", "star"]:
		string = "{outpath}/{ID}/RSEM/{sample}.isoforms.results"
		transcript_ids.append(string)
	return transcript_ids

def name_id_file(aligner):
	id_file = list()
	if aligner in ["Kallisto", "KALLISTO", "kallisto"]:
		string = "counts/kallisto/transcripts_ids.txt".format(outpath, ID)
		id_file.append(string)
	elif config["aligner"] in ["STAR", "Star", "star"]:
		string = "{}/{}/RSEM/transcripts_ids.txt".format(outpath, ID)
		id_file.append(string)
	return id_file

def get_quant_outs(aligner):
	quant_outs = list()
	if aligner in ["Kallisto", "KALLISTO", "kallisto"]:
		string = "{outpath}/{ID}/kallisto/{sample}/abundance.tsv",
		quant_outs.append(string)
	elif config["aligner"] in ["STAR", "Star", "star"]:
		string = "{outpath}/{ID}/RSEM/{sample}.isoforms.results"
		quant_outs.append(string)
	return quant_outs

def name_counts_file(aligner):
	counts_file = list()
	if aligner in ["Kallisto", "KALLISTO", "kallisto"]:
		string = "{outpath}/{ID}/counts/kallisto/{sample}.transcripts.expected.counts"
		counts_file.append(string)
	elif config["aligner"] in ["STAR", "Star", "star"]:
		string = "{outpath}/{ID}/counts/RSEM/{sample}.transcripts.expected.counts"
		counts_file.append(string)
	return counts_file

def get_fq1_qc(wildcards):
	return samples.loc[(wildcards.sample), ["fq1"]].dropna()

def get_fq2_qc(wildcards):
	if not is_single_end(wildcards.sample):
		return samples.loc[(wildcards.sample), ["fq2"]].dropna()
