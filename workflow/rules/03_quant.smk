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
		temp("{outpath}/{ID}/RSEM/{sample}.isoforms.results"),
		temp("{outpath}/{ID}/RSEM/{sample}.genes.results")
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
		"paste {input.gene_ids} {input.ind_counts} > {output.counts}"

rule aggregate_gene_tpm_counts:
	input:
		gene_ids = "{}/{}/RSEM/gene_ids.txt".format(outpath, ID),
		ind_counts = expand('{outpath}/{ID}/RSEM/{sample}.genes.tpm.counts', 
				outpath = outpath, ID = ID, sample = list_of_samples)
	output:
		counts = "{outpath}/{ID}/RSEM/genes.tpm_counts.tsv"
	shell:
		"paste {input.gene_ids} {input.ind_counts} > {output.counts}"

### RSEM (TPM - transcript) ###

rule create_transcript_ids_rsem:
	input:
		template = expand("{outpath}/{ID}/RSEM/{sample}.isoforms.results", 
				outpath = outpath, ID = ID, sample = list_of_samples[0])
	output:
		temp("{}/{}/RSEM/transcripts_ids.txt".format(outpath, ID))
	shell:
		"tail -n +2 {input.template} | cut -f1 | sed '1i \n' > "
		"{output}"

rule create_transcript_expected_counts:
	input:
		ind_counts = "{outpath}/{ID}/RSEM/{sample}.isoforms.results"
	output:
		temp("{outpath}/{ID}/RSEM/{sample}.transcripts.expected.counts")
	shell:
		"tail -n +2 {input.ind_counts} | cut -f5 | sed '1i {wildcards.sample}' > "
		"{output}"

rule create_transcript_tpm_counts:
	input:
		ind_counts = "{outpath}/{ID}/RSEM/{sample}.isoforms.results"
	output:
		temp("{outpath}/{ID}/RSEM/{sample}.transcripts.tpm.counts")
	shell:
		"tail -n +2 {input.ind_counts} | cut -f6 | sed '1i {wildcards.sample}' > "
		"{output}"

rule aggregate_transcript_expected_counts:
	input:
		transcripts_ids = "{}/{}/RSEM/transcripts_ids.txt".format(outpath, ID),
		ind_counts = expand('{outpath}/{ID}/RSEM/{sample}.transcripts.expected.counts', 
				outpath = outpath, ID = ID, sample = list_of_samples)
	output:
		counts = "{outpath}/{ID}/RSEM/transcripts.expected_counts.tsv"
	shell:
		"paste {input.transcripts_ids} {input.ind_counts} > {output.counts}"

rule aggregate_transcript_tpm_counts:
	input:
		transcripts_ids = "{}/{}/RSEM/transcripts_ids.txt".format(outpath, ID),
		ind_counts = expand('{outpath}/{ID}/RSEM/{sample}.transcripts.tpm.counts', 
				outpath = outpath, ID = ID, sample = list_of_samples)
	output:
		counts = "{outpath}/{ID}/RSEM/transcripts.tpm_counts.tsv"
	shell:
		"paste {input.transcripts_ids} {input.ind_counts} > {output.counts}"
