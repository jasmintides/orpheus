rule create_transcript_ids:
	input: expand(*get_transcript_counts(aligner), outpath = outpath, ID = ID, sample = template_sample)
	output: temp(*name_id_file(aligner))
	shell: "tail -n +2 {input} | cut -f1 | sed '1i \n' > {output}"  

rule create_transcript_counts:
	input: *get_quant_outs(aligner)
	output: 
		expected_file = temp(*name_expected_file(aligner)),
		tpm_file = temp(*name_tpm_file(aligner))
	run:
		if aligner.lower() == "kallisto":
			shell("tail -n +2 {input} | cut -f4 | "
				"sed '1i {wildcards.sample}' > {output.expected_file} ; "
				"tail -n +2 {input} | cut -f5 | "
				"sed '1i {wildcards.sample}' > {output.tpm_file}")
		elif aligner.lower() == "star":
			shell("tail -n +2 {input} | cut -f5 | "
				"sed '1i {wildcards.sample}' > {output.expected_file} ; "
				"tail -n +2 {input} | cut -f6 | "
				"sed '1i {wildcards.sample}' > {output.tpm_file}")

rule aggregate_transcript_counts:
	input:
		*name_id_file(aligner),
		all_expected = expand(name_expected_file(aligner),
			outpath = outpath, ID = ID, sample = list_of_samples),
		all_tpm = expand(name_tpm_file(aligner),
			outpath = outpath, ID = ID, sample = list_of_samples)
	output: 
		**get_full_counts(aligner)
	shell:
		"paste {input[0]} {input.all_expected} > {output.expected} ; "
		"paste {input[0]} {input.all_tpm} > {output.tpm}"
