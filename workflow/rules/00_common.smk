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
	my_dict = dict()
	if config["skip_trimming"]:
		my_dict["fq1"] = samples.loc[wildcards.sample, ["fq1"]].dropna().values[0]
		if not is_single_end(wildcards.sample):
			my_dict["fq2"] = samples.loc[wildcards.sample, ["fq2"]].dropna().values[0]
	else:
		my_dict["fq1"] = "{outpath}/{ID}/trimmed/{sample}_1.fq.gz"
		if not is_single_end(wildcards.sample):
			my_dict["fq2"] = "{outpath}/{ID}/trimmed/{sample}_2.fq.gz"
	return my_dict

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
		string = "{}/{}/kallisto/template/transcripts_ids.txt".format(outpath, ID)
		id_file.append(string)
	elif config["aligner"] in ["STAR", "Star", "star"]:
		string = "{}/{}/RSEM/template/transcripts_ids.txt".format(outpath, ID)
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
		string = "{outpath}/{ID}/kallisto/{sample}.transcripts.expected.counts"
		counts_file.append(string)
	elif config["aligner"] in ["STAR", "Star", "star"]:
		string = "{outpath}/{ID}/RSEM/{sample}.transcripts.expected.counts"
		counts_file.append(string)
	return counts_file

def name_all_counts_files(wildcards):
	counts_dict = dict()
	if aligner in ["Kallisto", "KALLISTO", "kallisto"]:
		counts_dict["expected"] = "{outpath}/{ID}/kallisto/{sample}.transcripts.expected.counts"
		counts_dict["tpm"] = "{outpath}/{ID}/kallisto/{sample}.transcripts.tpm.counts"
	elif config["aligner"] in ["STAR", "Star", "star"]:
		counts_dict["expected"] = "{outpath}/{ID}/kallisto/{sample}.transcripts.expected.counts"
		counts_dict["tpm"] = "{outpath}/{ID}/kallisto/{sample}.transcripts.tpm.counts"
	return counts_dict

def get_kallisto_index(wildcards):
	my_dict = dict()
	if config["kallisto_index"] != "":
		my_dict["index"] = config["kallisto_index"]
	else:
		my_dict["index"] = "{}/{}/kallisto/{}.idx".format(outpath, ID, build)
	return my_dict

def get_star_index(wildcards):
	my_list = list()
	if config["star_index"] != "":
		my_list.append(config["star_index"])
	else:
		my_list.append("{}/{}/STAR/{}".format(outpath, ID, build))
	return my_list

def get_full_counts(aligner):
	full_dict = dict()
	if aligner in ["Kallisto", "KALLISTO", "kallisto"]:
		full_dict["expected"] = "{outpath}/{ID}/kallisto/transcripts.expected.counts.tsv"
		full_dict["tpm"] = "{outpath}/{ID}/kallisto/transcripts.tpm.counts.tsv"
	elif config["aligner"] in ["STAR", "Star", "star"]:
		full_dict["expected"] = "{outpath}/{ID}/RSEM/transcripts.expected.counts.tsv"
		full_dict["tpm"] = "{outpath}/{ID}/RSEM/transcripts.tpm.counts.tsv"
	return full_dict

def get_fq1_qc(wildcards):
	return samples.loc[(wildcards.sample), ["fq1"]].dropna()

def get_fq2_qc(wildcards):
	if not is_single_end(wildcards.sample):
		return samples.loc[(wildcards.sample), ["fq2"]].dropna()
