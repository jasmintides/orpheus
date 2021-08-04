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

def get_trimmed(wildcards):
	if not is_single_end(wildcards.sample):
		return {'fq1': expand("{outpath}/{ID}/trimmed/{sample}_{group}.fq.gz", 
			group = 1, **wildcards),
			'fq2': expand("{outpath}/{ID}/trimmed/{sample}_{group}.fq.gz", 
			group = 2, **wildcards)}
	else:
		return {'fq1': "{outpath}/{ID}/trimmed/{sample}.fq.gz".format(**wildcards)}

myoutput = list()

if config["aligner"] in ["Kallisto", "KALLISTO", "kallisto"]:
	myoutput.append("{outpath}/{ID}/kallisto/outs/{sample}/abundance.tsv")
elif config["aligner"] in ["STAR", "Star", "star"]:
	myoutput.append("{outpath}/{ID}/RSEM/{sample}.transcripts.expected.counts")

transcript_ids = list()

def get_transcript_ids(aligner):
	if aligner in ["Kallisto", "KALLISTO", "kallisto"]:
		return {'template':
		template = expand("{outpath}/{ID}/RSEM/{sample}.isoforms.results", 
