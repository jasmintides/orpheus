rule rsem_prepare_reference:
	input:
		fasta = config["ref"]["fa"],
		gtf = config["ref"]["gtf"]
	params:
		ID = config["ID"],
		build = config["ref"]["build"]
	output:
		#temp("outs/{}/ref/{}.transcripts.fa".format(config["ID"], config["ref"]["build"])),
		directory("outs/{}/ref".format(config["ID"]))
	conda:
		"../envs/quant.yaml"
	threads:
		4
	shell:
		"rm -rf outs/{params.ID}/ref && "
		"mkdir outs/{params.ID}/ref && "
		"rsem-prepare-reference --num-threads {threads} --gtf {input.gtf} "
		"{input.fasta} outs/{params.ID}/ref/{params.build}"

rule rsem_calculate_expression_pe:
	input:
		bam = "outs/{ID}/star/{sample}/Aligned.toTranscriptome.out.bam",	
		ref = directory("outs/{}/ref".format(config["ID"]))
	params:
		is_single_end = lambda wildcards: is_single_end(wildcards.sample),
		ID = config["ID"],
		build = config["ref"]["build"]
	output:
		"outs/{ID}/RSEM/{sample}.isoforms.results",
		"outs/{ID}/RSEM/{sample}.genes.results"
	benchmark:
		"benchmarks/{ID}/quant/34_rsem.{sample}.txt"
	log:
		"logs/{ID}/rsem/01_rsem.{sample}.log"
	threads:
		5
	conda:
		"../envs/quant.yaml"
	shell:
		"is_single_end={params.is_single_end} ; if [[ $is_single_end == False ]]; then "
		"rsem-calculate-expression --num-threads {threads} "
		"--fragment-length-max 1000 --no-bam-output --paired-end "
		"--bam {input.bam} {input.ref}/{params.build} "
		"outs/{ID}/RSEM/{wildcards.sample} ; "
		"elif [[ $is_single_end == True ]]; then "
		"rsem-calculate-expression --num-threads {threads} "
		"--fragment-length-max 1000 --no-bam-output "
		"--bam {input.bam} {input.ref}/{params.build} "
		"outs/{ID}/RSEM/{wildcards.sample} ; fi"

rule create_raw_counts_star:
	input:
		ind_counts = "outs/{ID}/star/{sample}/ReadsPerGene.out.tab"
	output:
		temp("outs/{ID}/star/{sample}/{sample}.counts")
	shell:
		"tail -n +5 {input.ind_counts} | cut -f4 | sed '1i {wildcards.sample}' > "
		"outs/{wildcards.ID}/star/{wildcards.sample}/{wildcards.sample}.counts"

rule create_raw_counts_table:
	input:
		gene_ids = "outs/{}/star/gene_ids.txt".format(config["ID"]),
		ind_counts = expand('outs/{ID}/star/{sample}/{sample}.counts', ID = ID, sample = list_of_samples)
	params:
		ID = config["ID"]
	output:
		counts = "outs/{ID}/star/raw_counts.tsv"
	shell:
		"paste outs/{params.ID}/star/gene_ids.txt {input.ind_counts} > {output.counts}"

###### RSEM (TPM - gene) ######

rule create_gene_ids_rsem:
	input:
		template = expand("outs/{ID}/RSEM/{sample}.genes.results", ID = config["ID"], sample = list_of_samples[0])
	params:
		ID = config["ID"]
	output:
		temp("outs/{}/RSEM/gene_ids.txt".format(config["ID"]))
	shell:
		"tail -n +2 {input.template} | cut -f1 | sed '1i \n' > "
		"{output}"

rule create_tpm_counts:
	input:
		ind_counts = "outs/{ID}/RSEM/{sample}.genes.results"
	output:
		temp("outs/{ID}/RSEM/{sample}.genes.counts")
	shell:
		"tail -n +2 {input.ind_counts} | cut -f6 | sed '1i {wildcards.sample}' > "
		"{output}"
#		"outs/{wildcards.ID}/RSEM/{wildcards.sample}.genes.counts"

rule aggregate_tpm_counts_table:
	input:
		gene_ids = "outs/{}/RSEM/gene_ids.txt".format(config["ID"]),
		ind_counts = expand('outs/{ID}/RSEM/{sample}.genes.counts', ID = ID, sample = list_of_samples)
	params:
		ID = config["ID"]
	output:
		counts = "outs/{ID}/RSEM/tpm_counts.tsv"
	shell:
		"paste outs/{params.ID}/RSEM/gene_ids.txt {input.ind_counts} > {output.counts}"

### RSEM (TPM - transcript) ###

rule create_isoforms_ids_rsem:
	input:
		template = expand("outs/{ID}/RSEM/{sample}.isoforms.results", ID = config["ID"], sample = list_of_samples[0])
	params:
		ID = config["ID"]
	output:
		temp("outs/{}/RSEM/isoforms_ids.txt".format(config["ID"]))
	shell:
		"tail -n +2 {input.template} | cut -f1 | sed '1i \n' > "
		"{output}"
#		"outs/{params.ID}/RSEM/isoforms_ids.txt"

rule create_tpm_isoforms_counts:
	input:
		ind_counts = "outs/{ID}/RSEM/{sample}.isoforms.results"
	output:
		temp("outs/{ID}/RSEM/{sample}.isoforms.counts")
	shell:
		"tail -n +2 {input.ind_counts} | cut -f6 | sed '1i {wildcards.sample}' > "
		"{output}"
#		"outs/{wildcards.ID}/RSEM/{wildcards.sample}.isoforms.counts"

rule aggregate_tpm_isoforms_counts:
	input:
		gene_ids = "outs/{}/RSEM/isoforms_ids.txt".format(config["ID"]),
		ind_counts = expand('outs/{ID}/RSEM/{sample}.isoforms.counts', ID = ID, sample = list_of_samples)
	params:
		ID = config["ID"]
	output:
		counts = "outs/{ID}/RSEM/tpm_isoforms_counts.tsv"
	shell:
		"paste outs/{params.ID}/RSEM/isoforms_ids.txt {input.ind_counts} > {output.counts}"
