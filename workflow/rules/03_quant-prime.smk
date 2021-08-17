rule create_transcript_ids:
	input: expand(*get_transcript_counts(aligner), outpath = outpath, ID = ID, sample = template_sample)
	output: *name_id_file(aligner)
	shell: "tail -n +2 {input} | cut -f1 | sed '1i \n' > {output}"  

rule create_transcript_counts:
	input: *get_quant_outs(aligner)
	output: *name_counts_file(aligner)
	run:
		if aligner in ["kallisto", "KALLISTO", "Kallisto"]:
			shell("tail -n +2 {input} | cut -f4 | "
				"sed '1i {wildcards.sample}' > {output}")
		elif aligner in ["STAR", "star", "Star"]:
			shell("tail -n +2 {input.ind_counts} | cut -f5 | "
				"sed '1i {wildcards.sample}' > {output}")

#rule aggregate_transcript_counts:
#	input:
#		*name_id_file(aligner),
#		expand(*name_counts_file(aligner), 
#			outpath = outpath, ID = ID, sample = list_of_samples)
#	output: myoutput
#	shell:
#		"paste {input[0]} {input[1]} > {output}"
