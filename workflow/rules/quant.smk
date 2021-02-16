rule rsem_prepare_reference:
	input:
		fasta = config["ref"]["fa"],
		gtf = config["ref"]["gtf"]
	params:
		outpath = config["outpath"],
		ID = config["ID"],
		build = config["ref"]["build"]
	output:
		directory("{}/{}/ref".format(outpath, ID))
	conda:
		"../envs/quant.yaml"
	threads:
		4
	shell:
		"rm -rf {params.outpath}/{params.ID}/ref && "
		"mkdir {params.outpath}/{params.ID}/ref && "
		"rsem-prepare-reference --num-threads {threads} --gtf {input.gtf} "
		"{input.fasta} {params.outpath}/{params.ID}/ref/{params.build}"

rule rsem_calculate_expression:
	input:
		bam = "{outpath}/{ID}/star/{sample}/Aligned.toTranscriptome.out.bam",	
		ref = directory("{}/{}/ref".format(outpath, ID))
	params:
		is_single_end = lambda wildcards: is_single_end(wildcards.sample),
		build = config["ref"]["build"]
	output:
		"{outpath}/{ID}/RSEM/{sample}.isoforms.results",
		"{outpath}/{ID}/RSEM/{sample}.genes.results"
	threads:
		4
	conda:
		"../envs/quant.yaml"
	shell:
		"is_single_end={params.is_single_end} ; if [[ $is_single_end == False ]]; then "
		"rsem-calculate-expression --num-threads {threads} "
		"--fragment-length-max 1000 --no-bam-output --paired-end "
		"--bam {input.bam} {input.ref}/{params.build} "
		"{wildcards.outpath}/{wildcards.ID}/RSEM/{wildcards.sample} ; "
		"elif [[ $is_single_end == True ]]; then "
		"rsem-calculate-expression --num-threads {threads} "
		"--fragment-length-max 1000 --no-bam-output "
		"--bam {input.bam} {input.ref}/{params.build} "
		"{wildcards.outpath}/{wildcards.ID}/RSEM/{wildcards.sample} ; fi"

rule create_raw_counts_star:
	input:
		ind_counts = "{outpath}/{ID}/star/{sample}/ReadsPerGene.out.tab"
	output:
		temp("{outpath}/{ID}/star/{sample}/{sample}.counts")
	shell:
		"tail -n +5 {input.ind_counts} | cut -f4 | sed '1i {wildcards.sample}' > "
		"{wildcards.outpath}/{wildcards.ID}/star/{wildcards.sample}/{wildcards.sample}.counts"

rule create_raw_counts_table:
	input:
		gene_ids = "{}/{}/star/gene_ids.txt".format(outpath, ID),
		ind_counts = expand('{outpath}/{ID}/star/{sample}/{sample}.counts', 
				outpath = outpath, ID = ID, sample = list_of_samples)
	output:
		counts = "{outpath}/{ID}/star/raw_counts.tsv"
	shell:
		"paste {wildcards.outpath}/{wildcards.ID}/star/gene_ids.txt {input.ind_counts} > {output.counts}"

###### RSEM (expected counts, TPM - gene) ######

rule create_gene_ids_rsem:
	input:
		template = expand("{outpath}/{ID}/RSEM/{sample}.genes.results", 
				outpath = outpath, ID = config["ID"], sample = list_of_samples[0])
	output:
		temp("{}/{}/RSEM/gene_ids.txt".format(outpath, ID))
	shell:
		"tail -n +2 {input.template} | cut -f1 | sed '1i \n' > "
		"{output}"

rule create_gene_expected_counts:
	input:
		ind_counts = "{outpath}/{ID}/RSEM/{sample}.genes.results"
	output:
		temp("{outpath}/{ID}/RSEM/{sample}.genes.expected.counts")
	shell:
		"tail -n +2 {input.ind_counts} | cut -f5 | sed '1i {wildcards.sample}' > "
		"{output}"

rule create_gene_tpm_counts:
	input:
		ind_counts = "{outpath}/{ID}/RSEM/{sample}.genes.results"
	output:
		temp("{outpath}/{ID}/RSEM/{sample}.genes.tpm.counts")
	shell:
		"tail -n +2 {input.ind_counts} | cut -f6 | sed '1i {wildcards.sample}' > "
		"{output}"

rule aggregate_gene_expected_counts:
	input:
		gene_ids = "{}/{}/RSEM/gene_ids.txt".format(outpath, ID),
		ind_counts = expand('{outpath}/{ID}/RSEM/{sample}.genes.expected.counts', 
				outpath = outpath, ID = ID, sample = list_of_samples)
	output:
		counts = "{outpath}/{ID}/RSEM/genes.expected_counts.tsv"
	shell:
		"paste {wildcards.outpath}/{wildcards.ID}/RSEM/gene_ids.txt "
		"{input.ind_counts} > {output.counts}"

rule aggregate_gene_tpm_counts:
	input:
		gene_ids = "{}/{}/RSEM/gene_ids.txt".format(outpath, ID),
		ind_counts = expand('{outpath}/{ID}/RSEM/{sample}.genes.tpm.counts', 
				outpath = outpath, ID = ID, sample = list_of_samples)
	output:
		counts = "{outpath}/{ID}/RSEM/genes.tpm_counts.tsv"
	shell:
		"paste {wildcards.outpath}/{wildcards.ID}/RSEM/gene_ids.txt "
		"{input.ind_counts} > {output.counts}"

### RSEM (TPM - transcript) ###

rule create_isoform_ids_rsem:
	input:
		template = expand("{outpath}/{ID}/RSEM/{sample}.isoforms.results", 
				outpath = outpath, ID = ID, sample = list_of_samples[0])
	output:
		temp("{}/{}/RSEM/isoforms_ids.txt".format(outpath, ID))
	shell:
		"tail -n +2 {input.template} | cut -f1 | sed '1i \n' > "
		"{output}"

rule create_isoform_expected_counts:
	input:
		ind_counts = "{outpath}/{ID}/RSEM/{sample}.isoforms.results"
	output:
		temp("{outpath}/{ID}/RSEM/{sample}.isoforms.expected.counts")
	shell:
		"tail -n +2 {input.ind_counts} | cut -f5 | sed '1i {wildcards.sample}' > "
		"{output}"

rule create_isoform_tpm_counts:
	input:
		ind_counts = "{outpath}/{ID}/RSEM/{sample}.isoforms.results"
	output:
		temp("{outpath}/{ID}/RSEM/{sample}.isoforms.tpm.counts")
	shell:
		"tail -n +2 {input.ind_counts} | cut -f6 | sed '1i {wildcards.sample}' > "
		"{output}"

rule aggregate_isoform_expected_counts:
	input:
		gene_ids = "{}/{}/RSEM/isoforms_ids.txt".format(outpath, ID),
		ind_counts = expand('{outpath}/{ID}/RSEM/{sample}.isoforms.expected.counts', 
				outpath = outpath, ID = ID, sample = list_of_samples)
	output:
		counts = "{outpath}/{ID}/RSEM/isoforms.expected_counts.tsv"
	shell:
		"paste {wildcards.outpath}/{wildcards.ID}/RSEM/isoforms_ids.txt "
		"{input.ind_counts} > {output.counts}"

rule aggregate_isoform_tpm_counts:
	input:
		gene_ids = "{}/{}/RSEM/isoforms_ids.txt".format(outpath, ID),
		ind_counts = expand('{outpath}/{ID}/RSEM/{sample}.isoforms.tpm.counts', 
				outpath = outpath, ID = ID, sample = list_of_samples)
	output:
		counts = "{outpath}/{ID}/RSEM/isoforms.tpm_counts.tsv"
	shell:
		"paste {wildcards.outpath}/{wildcards.ID}/RSEM/isoforms_ids.txt "
		"{input.ind_counts} > {output.counts}"
