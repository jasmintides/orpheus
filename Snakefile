from snakemake.utils import min_version
min_version("6.0")

import pandas as pd
import re

configfile: "config/config.test.yaml"
wildcard_constraints: sample = '[^_\W]+'

samples = pd.read_table(config["samples"], dtype = str)
samples.set_index("sample", drop = False)
samples.set_index("fq1", drop = False)
list_of_samples = samples["sample"].tolist()
template_sample = list_of_samples[0]


ID = config['ID']
outpath = config['outpath']
build = config["ref"]["build"]
aligner = config["aligner"]
samplez = "sample_1"
skip_trimming = config["skip_trimming"]

myoutput = list()

if config["aligner"] in ["Kallisto", "KALLISTO", "kallisto"]:
#	myoutput.append("{outpath}/{ID}/kallisto/transcripts.expected.counts.tsv")
#	myoutput.append("{outpath}/{ID}/kallisto/{sample}/abundance.tsv")
	myoutput.append("{outpath}/{ID}/trimmed/{sample}_{groups}.fq.gz")
elif config["aligner"] in ["STAR", "Star", "star"]:
	myoutput.append("{outpath}/{ID}/RSEM/transcripts.expected.counts")

rule all:
	input:
		expand(myoutput, outpath = outpath, ID = ID, sample = template_sample, groups = [1, 2])

include: "workflow/rules/00_common.smk"
include: "workflow/rules/01_trim.smk"
#include: "workflow/rules/02_align.smk"
#include: "workflow/rules/03_quant-prime.smk"
