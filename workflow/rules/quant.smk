rule rsem_prepare_reference:
	input:
		fasta = config["ref"]["fa"],
		gtf = config["ref"]["gtf"]
	params:
		ID = config["ID"],
		build = config["ref"]["build"]
	output:
		temp("outs/{}/ref/{}.transcripts.fa".format(config["ID"], config["ref"]["build"]))
	conda:
		"../envs/quant.yaml"
	threads:
		4
	shell:
		"rsem-prepare-reference --num-threads {threads} --gtf {input.gtf} "
		"{input.fasta} outs/{params.ID}/ref/{params.build}"

rule rsem_calculate_expression:
	input:
		bam = "outs/{ID}/star/{sample}/Aligned.toTranscriptome.out.bam"
	params:
		ID = config["ID"],
		build = config["ref"]["build"]
	output:
		temp("outs/{ID}/RSEM/{sample}.isoforms.results"),
		temp("outs/{ID}/RSEM/{sample}.genes.results")
	benchmark:
		"benchmarks/{ID}/quant/01_rsem.{sample}.txt"
	log:
		"logs/{ID}/rsem/01_rsem.{sample}.log"
	threads:
		4
	conda:
		"../envs/quant.yaml"
	shell:
		"rsem-calculate-expression --num-threads {threads} "
		"--fragment-length-max 1000 --no-bam-output --paired-end "
		"--bam {input.bam} outs/{params.ID}/ref/{params.build} "
		"outs/{ID}/RSEM/{wildcards.sample}"

rule create_gene_ids_star:
	input:
		template = "outs/{}/star/{}/ReadsPerGene.out.tab".format(config["ID"], list_of_samples[0])
	params:
		ID = config["ID"]
	output:
		temp("outs/{}/star/gene_ids.txt".format(config["ID"]))
	shell:
		"tail -n +5 {input.template} | cut -f1 | sed '1i \n' > "
		"outs/{params.ID}/star/gene_ids.txt"

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

###### RSEM (TPM) ######

rule create_gene_ids_rsem:
	input:
		template = "outs/{}/RSEM/{}.genes.results".format(config["ID"], list_of_samples[0])
	params:
		ID = config["ID"]
	output:
		temp("outs/{}/RSEM/gene_ids.txt".format(config["ID"]))
	shell:
		"tail -n +2 {input.template} | cut -f1 | sed '1i \n' > "
		"outs/{params.ID}/RSEM/gene_ids.txt"

rule create_tpm_counts:
	input:
		ind_counts = "outs/{ID}/RSEM/{sample}.genes.results"
	output:
		temp("outs/{ID}/RSEM/{sample}.counts")
	shell:
		"tail -n +2 {input.ind_counts} | cut -f6 | sed '1i {wildcards.sample}' > "
		"outs/{wildcards.ID}/RSEM/{wildcards.sample}.counts"

rule aggregate_tpm_counts_table:
	input:
		gene_ids = "outs/{}/RSEM/gene_ids.txt".format(config["ID"]),
		ind_counts = expand('outs/{ID}/RSEM/{sample}.counts', ID = ID, sample = list_of_samples)
	params:
		ID = config["ID"]
	output:
		counts = "outs/{ID}/RSEM/tpm_counts.tsv"
	shell:
		"paste outs/{params.ID}/RSEM/gene_ids.txt {input.ind_counts} > {output.counts}"
		
