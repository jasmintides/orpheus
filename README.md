<h1>orpheus: RNA-seq quantification</h1>

<p align="center">
<img src="img/Orpheus_schema.svg">
</p>

**orpheus** is a Snakemake workflow that performs quality control, alignment, and 
quantification on raw sequencing data. Processing big data with limited 
computational resources presents a challenge to time-sensitive exploratory
data science. The Data Science team has developed a modular, scalable pipeline 
that employs a suite of bioinformatic tasks into a concise workflow that 
processes several samples at once while being efficient in both runtime and
computational resource. The result is a standardized data output that is 
analysis-ready and easily communicated to other team members and departments. 

In-depth explanation of the workflow steps, file structure, benchmarking, 
rationale for software and techniques chosen, and other information 
may be found on the **orpheus** entry at the 
[Data Science DokuWiki](https://hpc.agios.local/dokuwiki/doku.php?id=orpheus_rna-seq).