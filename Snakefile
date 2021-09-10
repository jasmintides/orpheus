from snakemake.utils import min_version
min_version("6.0")

import pandas as pd
configfile: "config/config.test.yaml"

samples = pd.read_table(config["samples"], dtype = str).set_index("sample", drop = False)
list_of_samples = samples["sample"].tolist()
template_sample = list_of_samples[0]

ID = config["ID"]
outpath = config["outpath"]
build = config["ref"]["build"]
aligner = config["aligner"]
skip_trimming = config["skip_trimming"]

if aligner.lower() == "kallisto":
	if "kallisto_index" not in list(config.keys()):
		config["kallisto_index"] = ""
if aligner.lower() == "star":
	if "STAR_index" not in list(config.keys()):
		config["star_index"] = ""

rule all:
	input:
		"{}/{}/SummExp/{}.SummExp.Rds".format(outpath, ID, ID),
		"{}/{}/qc/multiqc_report.{}.html".format(outpath, ID, ID)

include: "workflow/rules/00_common.smk"
include: "workflow/rules/01_qc.smk"
include: "workflow/rules/01p_trim.smk"
include: "workflow/rules/02_quant.smk"
include: "workflow/rules/03_aggregate.smk"
include: "workflow/rules/04_write_se.smk"
