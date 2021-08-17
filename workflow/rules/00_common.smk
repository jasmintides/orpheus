def is_single_end(sample):
	return pd.isnull(samples.loc[(sample), "fq2"])

def get_fastq(sample):
	fastq_dict = dict()
	fastq_dict["r1"] = samples[(samples["sample"] == sample), "fq1" & samples["fq1" ==]
#	if not is_single_end(sample):
#		fastq_dict["r2"] = samples.loc[samples.index.intersection(sample, "fq2")]
	return fastq_dict

def get_trimmed(w):
	trimmed_dict = dict()
	string = "{}/{}/trimmed/{}_{}.fq.gz"
	trimmed_dict["r1"] = string.format(w.outpath, w.ID, w.sample, 1)
	if not is_single_end(w.sample):
		trimmed_dict["r2"] = string.format(w.outpath, w.ID, w.sample, 2)
	return trimmed_dict

def fastq_to_aligner(wildcards):
	if wildcards.skip_trimming:
		x = get_fastq(wildcards.sample)
	elif not skip_trimming:
		x = wildcards.get_trimmed(wildcards)
	return x

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
		string = "{outpath}/{ID}/kallisto/{sample}.abundance.tsv",
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

def get_fq1_qc(wildcards):
	return samples.loc[(wildcards.sample), ["fq1"]].dropna()

def get_fq2_qc(wildcards):
	if not is_single_end(wildcards.sample):
		return samples.loc[(wildcards.sample), ["fq2"]].dropna()
