rule replace_rg:
	input:
		"outs/{ID}/star/{sample}/Aligned.sortedByCoord.out.bam"
	output:
		temp("outs/{ID}/star/{sample}/Aligned.sortedByCoord.out.rgAligned.bam")
	benchmark:
		"benchmarks/{ID}/call/00_replace_rg.{sample}.txt"
	log:
		"logs/{ID}/picard/replace_rg/{sample}.log"
	params:
		"RGID={sample} RGLB={sample} RGPL={sample} RGPU={sample} RGSM={sample} "
		"VALIDATION_STRINGENCY=LENIENT"
	wrapper:
		"0.57.0/bio/picard/addorreplacereadgroups"

rule mark_duplicates:
	input:
		"outs/{ID}/star/{sample}/Aligned.sortedByCoord.out.rgAligned.bam"
	output:
		bam = "outs/{ID}/star/{sample}/Aligned.sortedByCoord.out.markedAligned.bam",
		metrics = "outs/{ID}/star/{sample}/metrics.txt"
	benchmark:
		"benchmarks/{ID}/call/01_mark_duplicates.{sample}.txt"
	log:
		"logs/{ID}/picard/dedup/{sample}.log"
	params:
		mem = "-Xmx8g"
	wrapper:
		"0.57.0/bio/picard/markduplicates"

rule split_n_cigar_reads:
	input:
		bam = "outs/{ID}/star/{sample}/Aligned.sortedByCoord.out.markedAligned.bam",
		ref = config['ref']['fa']
	output:
		temp("outs/{ID}/split/{sample}.bam")
	benchmark:
		"benchmarks/{ID}/call/02_split_n_cigar_reads.{sample}.txt"
	log:
		"logs/{ID}/gatk/splitNCIGARreads/{sample}.log"
	params:
		extra = "--tmp-dir outs/{ID}/star/{sample}",
		java_opts = "-Xmx4g"
	wrapper:
		"0.57.0/bio/gatk/splitncigarreads"

rule gatk_bqsr:
	input:
		bam = "outs/{ID}/split/{sample}.bam",
		ref = config['ref']['fa'],
		known = config["ref"]["known_sites"]
	output:
		bam = temp("outs/{ID}/recal/{sample}.bam")
	benchmark:
		"benchmarks/{ID}/call/03_gatk_bqsr.{sample}.txt"
	log:
		"logs/{ID}/gatk/bqsr/{sample}.log"
	params:
		extra = "-DF NotDuplicateReadFilter --tmp-dir outs/{ID}/split",
		java_opts = "-Xmx4g"
	wrapper:
		"0.57.0/bio/gatk/baserecalibrator"

rule haplotype_caller:
	input:
		bam = "outs/{ID}/recal/{sample}.bam",
		ref = config['ref']['fa']
	output:
		gvcf = temp("outs/{ID}/calls/{sample}.{chr_chunks}.g.vcf.gz")
	benchmark:
		"benchmarks/{ID}/call/04_haplotype_caller.{sample}.{chr_chunks}.txt"
	log:
		"logs/{ID}/gatk/haplotypecaller/{sample}.{chr_chunks}.log"
	threads:
		4
	params:
		chr_intervals = lambda wildcards: chr_dict[wildcards.chr_chunks],
		extra = "--dont-use-soft-clipped-bases true -DF NotDuplicateReadFilter "
			"--minimum-mapping-quality 0 --base-quality-score-threshold 10 -mbq 13",
		java_opts = "-Xmx8g"
	conda:
		"../envs/gatk4.yaml"
	shell:
		"gatk --java-options {params.java_opts} HaplotypeCaller {params.extra} "
		"{params.chr_intervals} -R {input.ref} -I {input.bam} -ERC GVCF -O {output.gvcf}"

rule combine_gvcfs:
	input:
		gvcfs = expand("outs/{ID}/calls/{{sample}}.{chr_chunks}.g.vcf.gz",
				ID = ID, chr_chunks = chr_chunks),
		ref = config['ref']['fa']
	output:
		gvcf = temp("outs/{ID}/combined/{sample}.g.vcf.gz")
	benchmark:
		"benchmarks/{ID}/call/05_combine_gvcfs.{sample}.txt"
	log:
		"logs/{ID}/gatk/combine_gvcfs/{sample}.log"
	params:
		extra = "",
		java_opts = ""
	wrapper:
		"0.58.0/bio/gatk/combinegvcfs"

rule genotype_gvcfs:
	input:
		gvcf = "outs/{ID}/combined/{sample}.g.vcf.gz",
		ref = config['ref']['fa']
	output:
		vcf = temp("outs/{ID}/calls/unfiltered/{sample}.unfiltered.vcf.gz")
	benchmark:
		"benchmarks/{ID}/call/06_genotype_gvcfs.{sample}.txt"
	log:
		"logs/{ID}/gatk/genotypegvcfs.{sample}.log"
	params:
		extra = "-stand-call-conf 0.0",
		java_opts = "",
	wrapper:
		"0.58.0/bio/gatk/genotypegvcfs"

rule gatk_filter:
	input:
		vcf = "outs/{ID}/calls/unfiltered/{sample}.unfiltered.vcf.gz",
		ref = config["ref"]["fa"]
	output:
		vcf = temp("outs/{ID}/calls/filtered/{sample}.vcf.gz")
	benchmark:
		"benchmarks/{ID}/call/07_gatk_filter.{sample}.txt"
	log:
		"logs/{ID}/gatk/filter/snvs.{sample}.log"
	params:
		filters = {"FS": "FS > 30.0", "QD": "QD < 2.0", "DP": "DP < 20"},
		extra = "-window 35 -cluster 3",
		java_opts = "",
	wrapper:
		"0.59.1/bio/gatk/variantfiltration"

rule snpeff:
	input:
		calls = "outs/{ID}/calls/filtered/{sample}.vcf.gz"
	output:
		calls = temp("outs/{ID}/annotated/{sample}.vcf")
	benchmark:
		"benchmarks/{ID}/call/08_snpeff.{sample}.txt"
	log:
		"logs/{ID}/snpeff/{sample}.log"
	params:
		extra = "-Xmx4g -no-downstream -no-intergenic -no-intron -no-upstream",
		reference = "GRCh37.75"
	conda:
		"../envs/snpeff.yaml"
	shell:
		"snpEff {params.extra} {params.reference} "
		"{input.calls} > {output.calls}"

rule bcftools_annotate:
        input:
                calls = "outs/{ID}/annotated/{sample}.vcf",
                bed = "/data/exploratory/Users/jeff.alvarez/pipeline_ins/Alu.RepeatMasker.hg19.ID.bed",
                header = "/data/exploratory/Users/jeff.alvarez/pipeline_ins/Alu.RepeatMasker.hg19.ID.txt"
        output:
                vcf = "outs/{ID}/final/{sample}.vcf.gz"
        benchmark:
                "benchmarks/{ID}/call/{sample}.09_bcftools_annotate.txt"
        log:
                "logs/{ID}/09_bcftools_annotate/{sample}.log"
        params:
                columns = "CHROM,FROM,TO,ALU_NAME,ALU_ID,STRAND"
        conda:
                "../envs/bcftools.yaml"
        shell:
                """bcftools annotate -a {input.bed} -h {input.header} -c {params.columns} {input.calls} | """
                """bcftools filter -i "REF == 'A' & ALT == 'G' | REF == 'T' & ALT == 'C'" | """
		"""bcftools filter -e "FORMAT/AD[:1] < 5" | bgzip -c > {output.vcf}"""
