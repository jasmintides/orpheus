rule rsem_prepare_reference:
    input: fasta = config["ref"]["fa"], gtf = config["ref"]["gtf"]
    params: outpath = config["outpath"], ID = config["ID"], 
            build = config["ref"]["build"], aligner = config["aligner"]
    output: directory("{}/{}/{}/ref".format(outpath, ID, aligner))
    conda: "../envs/quant.yaml"
    threads: 4
    shell: "rm -rf {params.outpath}/{params.ID}/ref && "
    "mkdir {params.outpath}/{params.ID}/ref && "
    "rsem-prepare-reference --num-threads {threads} --gtf {input.gtf} "
    "{input.fasta} {params.outpath}/{params.ID}/ref/{params.build}"

rule rsem_calculate_expression:
    input:
        bam = "{outpath}/{ID}/{aligner}/{sample}/Aligned.toTranscriptome.out.bam",	
        ref = directory("{}/{}/{}/ref".format(outpath, ID, aligner))
    params:
        is_single_end = lambda wildcards: is_single_end(wildcards.sample),
        build = config["ref"]["build"]
    output:
        temp("{outpath}/{ID}/{aligner}/{sample}.isoforms.results"),
        temp("{outpath}/{ID}/{aligner}/{sample}.genes.results")
    threads: 4
    conda: "../envs/quant.yaml"
    shell:
        "is_single_end={params.is_single_end} ; if [[ $is_single_end == False ]]; then "
        "rsem-calculate-expression --num-threads {threads} "
        "--fragment-length-max 1000 --no-bam-output --paired-end "
        "--bam {input.bam} {input.ref}/{params.build} "
        "{wildcards.outpath}/{wildcards.ID}/{wildcards.aligner}/{wildcards.sample} ; "
        "elif [[ $is_single_end == True ]]; then "
        "rsem-calculate-expression --num-threads {threads} "
        "--fragment-length-max 1000 --no-bam-output "
        "--bam {input.bam} {input.ref}/{params.build} "
        "{wildcards.outpath}/{wildcards.ID}/{wildcards.aligner}/{wildcards.sample} ; fi"

rule create_gene_ids_rsem:
    input: get_transcript_ids
    output: temp("{}/{}/{}/gene_ids.txt".format(outpath, ID, aligner))
    shell:
        "tail -n +2 {input.template} | cut -f1 | sed '1i \n' > "
        "{output}"

rule create_gene_expected_counts:
    input:
        ind_counts = "{outpath}/{ID}/{aligner}/{sample}.genes.results"
    output:
        temp("{outpath}/{ID}/{aligner}/{sample}.genes.expected.counts")
    shell:
        "tail -n +2 {input.ind_counts} | cut -f5 | sed '1i {wildcards.sample}' > "
        "{output}"

rule create_gene_tpm_counts:
    input:
        ind_counts = "{outpath}/{ID}/{aligner}/{sample}.genes.results"
    output:
        temp("{outpath}/{ID}/{aligner}/{sample}.genes.tpm.counts")
    shell:
        "tail -n +2 {input.ind_counts} | cut -f6 | sed '1i {wildcards.sample}' > "
        "{output}"

rule aggregate_gene_expected_counts:
    input: gene_ids = "{}/{}/{}/gene_ids.txt".format(outpath, ID, aligner),
            ind_counts = expand('{outpath}/{ID}/{aligner}/{sample}.genes.expected.counts', 
            outpath = outpath, ID = ID, aligner = aligner, sample = list_of_samples)
    output: counts = "{outpath}/{ID}/{aligner}/genes.expected_counts.tsv"
    shell:
        "paste {input.gene_ids} {input.ind_counts} > {output.counts}"

rule aggregate_gene_tpm_counts:
    input:
        gene_ids = "{}/{}/{}/gene_ids.txt".format(outpath, ID, aligner),
        ind_counts = expand('{outpath}/{ID}/{aligner}/{sample}.genes.tpm.counts', 
            outpath = outpath, ID = ID, aligner = aligner, sample = list_of_samples)
    output:
        counts = "{outpath}/{ID}/{aligner}/genes.tpm_counts.tsv"
    shell:
        "paste {input.gene_ids} {input.ind_counts} > {output.counts}"

rule create_transcript_ids:
    input: unpack(get_transcript_ids)
    output: temp("{}/{}/counts/{}/transcripts_ids.txt".format(outpath, ID, aligner))
    shell: "tail -n +2 {input.transcript_ids} | cut -f1 | sed '1i \n' > {output}"  

rule create_transcript_counts:
    input: unpack(get_transcript_counts)
    output: 
        expected = temp("{outpath}/{ID}/counts/{aligner}/{sample}.transcripts.expected.counts"),
        tpm = temp("{outpath}/{ID}/counts/{aligner}/{sample}.transcripts.tpm.counts")
    run:
        if wildcards.aligner in ["kallisto", "KALLISTO", "Kallisto"]:
            shell("tail -n +2 {input.ind_counts} | cut -f4 | "
                "sed '1i {wildcards.sample}' > {output.expected} ; "
                "tail -n +2 {input.ind_counts} | cut -f5 | "
                "sed '1i {wildcards.sample}' > {output.tpm}")
        elif wildcards.aligner in ["STAR", "star", "Star"]:
            shell("tail -n +2 {input.ind_counts} | cut -f5 | "
                "sed '1i {wildcards.sample}' > {output.expected} ; "
                "tail -n +2 {input.ind_counts} | cut -f6 | "
                "sed '1i {wildcards.sample}' > {output.tpm}")

rule aggregate_transcript_counts:
    input:
        transcripts_ids = "{}/{}/counts/{}/transcripts_ids.txt".format(outpath, ID, aligner),
        all_expected = expand(
            '{outpath}/{ID}/counts/{aligner}/{all}.transcripts.expected.counts', 
            outpath = outpath, ID = ID, aligner = aligner, all = list_of_samples
            ),
        all_tpm = expand(
            '{outpath}/{ID}/counts/{aligner}/{all}.transcripts.tpm.counts', 
            outpath = outpath, ID = ID, aligner = aligner, all = list_of_samples
            )
    output: 
        expected_counts = "{outpath}/{ID}/counts/{aligner}/transcripts.expected.counts.tsv",
        tpm_counts = "{outpath}/{ID}/counts/{aligner}/transcripts.tpm.counts.tsv"
    shell:
        "paste {input.transcripts_ids} {input.all_expected} > {output.expected_counts} ; "
        "paste {input.transcripts_ids} {input.all_tpm} > {output.tpm_counts}"
