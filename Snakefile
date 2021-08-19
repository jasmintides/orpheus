from snakemake.utils import min_version
min_version("6.0")

import pandas as pd
import re
configfile: "config/config.test.yaml"
#wildcard_constraints: sample = '\w'

samples = pd.read_table(config["samples"], dtype = str).set_index("sample", drop = False)
list_of_samples = samples["sample"].tolist()
template_sample = list_of_samples[0]

ID = config["ID"]
outpath = config["outpath"]
build = config["ref"]["build"]
aligner = config["aligner"]
skip_trimming = config["skip_trimming"]

expected_counts = list()
tpm_counts = list()

if aligner in ["Kallisto", "KALLISTO", "kallisto"]:
	expected_counts.append("{outpath}/{ID}/kallisto/transcripts.expected.counts.tsv")
	tpm_counts.append("{outpath}/{ID}/kallisto/transcripts.tpm.counts.tsv")
elif aligner in ["STAR", "Star", "star"]:
	expected_counts.append("{outpath}/{ID}/RSEM/transcripts.expected.counts")
	tpm_counts.append("{outpath}/{ID}/RSEM/transcripts.tpm.counts")

rule all:
	input:
		expand(expected_counts, outpath = outpath, ID = ID),
		expand(tpm_counts, outpath = outpath, ID = ID)

include: "workflow/rules/00_common.smk"
include: "workflow/rules/01_trim.smk"
include: "workflow/rules/02_quant.smk"
include: "workflow/rules/03_aggregate.smk"
