rule replace_rg:
	input:
		'outs/{ID}/star/{sample}/Aligned.sortedByCoord.out.bam'
	output:
		temp("outs/{ID}/star/{sample}/Aligned.sortedByCoord.out.rgAligned.bam")
	benchmark:
		"benchmarks/{ID}/call/00_replace_rg.{sample}.txt"
	log:
		"logs/{ID}/picard/replace_rg/{sample}.log"
	params:
		"RGID={sample} RGLB={sample} RGPL={sample} RGPU={sample} RGSM={sample} SO=coordinate "
		"VALIDATION_STRINGENCY=SILENT"
	wrapper:
		"0.57.0/bio/picard/addorreplacereadgroups"

rule mark_duplicates:
	input:
		"outs/{ID}/star/{sample}/Aligned.sortedByCoord.out.rgAligned.bam"
	output:
		bam = temp("outs/{ID}/star/{sample}/Aligned.sortedByCoord.out.markedAligned.bam"),
		metrics = "outs/{ID}/star/{sample}/metrics.txt"
	benchmark:
		"benchmarks/{ID}/call/01_mark_duplicates.{sample}.txt"
	log:
		"logs/{ID}/picard/dedup/{sample}.log"
	params:
		""
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
		java_opts = ""
	wrapper:
		"0.57.0/bio/gatk/splitncigarreads"

rule gatk_bqsr:
	input:
		bam = "outs/{ID}/split/{sample}.bam",
		ref = config['ref']['fa'],
		known = config["known_sites"]
	output:
		bam = temp("outs/{ID}/recal/{sample}.bam")
	benchmark:
		"benchmarks/{ID}/call/03_gatk_bqsr.{sample}.txt"
	log:
		"logs/{ID}/gatk/bqsr/{sample}.log"
	params:
		extra = "-DF NotDuplicateReadFilter --tmp-dir outs/split",
		java_opts = ""
	wrapper:
		"0.57.0/bio/gatk/baserecalibrator"

rule haplotype_caller:
	input:
		bam = "outs/{ID}/recal/{sample}.bam",
		ref = config['ref']['fa']
	output:
		gvcf = temp("outs/{ID}/calls/{sample}.g.vcf.gz")
	benchmark:
		"benchmarks/{ID}/call/04_haplotype_caller.{sample}.txt"
	log:
		"logs/{ID}/gatk/haplotypecaller/{sample}.log"
	threads:
		4
	params:
		extra = "--dont-use-soft-clipped-bases true -DF NotDuplicateReadFilter -RF MappingQualityReadFilter --minimum-mapping-quality 0 --base-quality-score-threshold 13 -mbq 13 -L /data/exploratory/Users/jeff.alvarez/adar/editing_analysis/A549/input/data/Alu.filtered.bed --tmp-dir outs/recal",
		java_opts = ""
	wrapper:
		"0.57.0/bio/gatk/haplotypecaller"

rule genotype_gvcfs:
	input:
		gvcf = "outs/{ID}/calls/{sample}.g.vcf.gz",
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
		vcf = "outs/{ID}/calls/filtered/{sample}.vcf.gz"
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
