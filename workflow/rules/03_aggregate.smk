rule create_transcript_ids:
	input: expand(*get_transcript_counts(aligner), outpath = outpath, ID = ID, sample = template_sample)
	output: *name_id_file(aligner)
	shell: "tail -n +2 {input} | cut -f1 | sed '1i \n' > {output}"  

rule create_transcript_counts:
	input: *get_quant_outs(aligner)
	output: **name_all_counts_files(aligner)
	run:
	if aligner.lower() == "kallisto":
		shell("tail -n +2 {input} | cut -f4 | "
			"sed '1i {wildcards.sample}' > {output.expected} ; "
			"tail -n +2 {input} | cut -f5 | "
			"sed '1i {wildcards.sample}' > {output.tpm}")
	elif aligner.lower() == "star":
		shell("tail -n +2 {input} | cut -f5 | "
			"sed '1i {wildcards.sample}' > {output.expected} ; "
			"tail -n +2 {input} | cut -f6 | "
			"sed '1i {wildcards.sample}' > {output.tpm}")

rule aggregate_transcript_counts:
	input:
		*name_id_file(aligner),
		expected = expand("{outpath}/{ID}/kallisto/{sample}.transcripts.expected.counts",
			outpath = outpath, ID = ID, sample = list_of_samples),
		tpm = expand("{outpath}/{ID}/kallisto/{sample}.transcripts.tpm.counts",
			outpath = outpath, ID = ID, sample = list_of_samples)
	output: 
		expected_counts = expected_counts,
		tpm_counts = tpm_counts
	shell:
		"paste {input[0]} {input.expected} > {output.expected_counts} ; "
		"paste {input[0]} {input.tpm} > {output.tpm_counts}"
