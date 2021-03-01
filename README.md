<h1>Orpheus</h1>
Orpheus (Omicsoft-inspired RNA-seq pipeline) is an RNA-seq workflow that 
performs quality control, alignment, expression quantification and variant 
calling (based on
[GATK Best Practices](https://gatk.broadinstitute.org/hc/en-us/articles/360035531192-RNAseq-short-variant-discovery-SNPs-Indels-)). 
Orpheus was written with the Python-based workflow manager 
[Snakemake](https://snakemake.readthedocs.io/en/stable/), which is installed in 
a Conda environment in the Agios HPC. Users may test that the pipeline is 
working properly with the test data included in the repo:

```
module load conda
conda activate snakemake

git clone https://git.agios.local/Jeff.Alvarez/orpheus.git
cd /path/to/orpheus
snakemake --use-conda --cores 8
```

<h2>Overview</h2>
The git repository for Orpheus is laid out in the following
structure:

```
├── Snakefile
├── workflow
│   ├── envs
│   └── rules
│   └── scripts
├── config
│   └── config.local.yaml
│   └── config.docker.yaml
├── outs
├── benchmarks
├── rules
├── resources
├── .gitignore
├── README.md
└── LICENSE.md
```

The pipeline code for Orpheus is stored in the <code>workflow</code> directory. 
<code>workflow</code> contains subdirectories for <code>rules</code> written 
for each module of the pipeline: trimming, QC, alignment, quantification, and 
an optional step  for variant calling. Although most rules are run with 
Snakemake wrapper scripts, Conda environments for individual rules are stored 
in <code>envs</code>.

The <code>config</code> directory contains configuration forms with which 
users define analysis metadata and paths to reference data to be used in the 
Orpheus pipeline. There exists configuration files 
<code>config/config.local.yaml</code> and <code>config/config.docker.yaml</code>,
which are to be used with local and Docker deployments of the pipeline 
respectively. They are run with a test dataset included in the repo, 
located in the <code>data</code> directory. Configuration files are further
described in [Step 2](#step-2-configure-orpheus).

The central <code>Snakefile</code> marks the entrypoint of the pipeline by 
specifying a central rule for the expected output of the entire pipeline:

* SummarizedExperiment R object with RSEM gene counts results
* SummarizedExperiment R object with RSEM transcript counts results
* A MultiQC report aggregating FASTQC and STAR logs from all samples
</li>

To do this, the <code>Snakefile</code> parses user-defined info from the config 
file and include the appropriate modules from the <code>workflow</code> 
directory to generate the expected output defined in the 
<code>Snakefile</code>.

All output files generated in the workflow are stored in the <code>outs</code> 
directory, which contains subdirectories named after the analysis ID name 
specified in the config file of a given analysis. From the directory root, the 
raw counts, TPM counts, and MultiQC report are found in
<code>outs/{analysis_ID}/SummExp/{analysis_ID}.genes_SummExp.Rds</code>, 
<code>outs/{analysis_ID}/SummExp/{analysis_ID}.transcripts_SummExp.Rds</code>, 
<code>outs/{analysis_ID}/qc/multiqc_report.{analysis_ID}.html</code> 
respectively. Benchmarking data, which summarizes running time and memory 
usage, and system logs for each rule in the pipeline are also stored in a 
similar manner in <code>benchmarks</code> and <code>logs</code> 
respectively.

<h2>Usage</h2>
<h3>Step 1: Get Orpheus</h3>
[Clone](https://help.github.com/en/github/creating-cloning-and-archiving-repositories/cloning-a-repository)
the repo to your path of choice on HPC or a local machine:

```
git clone https://git.agios.local/Jeff.Alvarez/orpheus.git
```

Orpheus is also available as a Docker image, which may be run by following the 
directions at [Step 4](#step-4-execute workflow).

<h3>Step 2: Create sample sheet</h3>
Input fastqs and metadata are specified with a sample sheet. Each entry of the
sample sheet describes a unique sample ID, paths to corresponding fastq paths,
and metadata columns describing the experimental setup of the analysis. An 
example <code>data/sample_sheet.tsv</code> is included with Orpheus for analyses
run with the test dataset in the same directory. An empty template to write
in your own data is available at <code>data/template.tsv</code>. Alternatively,
you may create your own sample sheet in Microsoft Excel as long as the columns
are labeled exactly as shown below and the file is saved as "Text (Tab
delimited) (*.txt)"

![alt text](img/example_sample_sheet.png)

<h3>Step 3: Configure Orpheus</h3>
Analysis names and reference data are specified with a config file. The layout 
is based on the config file used for the 
[Array Studio RNA-seq pipeline](https://git.agios.local/Mark.Fletcher/array_studio_RNAseq_pipeline).
The example shown below is based on the config file 
<code>config/config.test.yaml</code>, which was made to run the test
data:

```
ID: 20XX_11_11_test
samples: /path/to/sample_sheet.tsv
organism: "Human"
ref:
        fa: /path/to/Human_B37.3.fasta
        gtf: /path/to/Human_B37.3.gtf
        build: Human_B37.3
        known_sites: /path/to/dbsnp_138.b37.vcf.gz
outpath: outs
```
The config file takes the following values:
* <b>ID</b>: Filename of output directory where analysis files are written.
* <b>samples</b>: Absolute path to sample sheet with input fastqs and metadata.
* <b>organism</b>: "Human" or "Mouse" - used for gene annotation of counts 
  data.
* <b>outpath</b>: Default is set to write output to <code>outs</code> directory 
  at top of Orpheus directory as a relative path, otherwise can be set to 
  absolute path of another location.
* <b>ref</b>: Reference data to be used in alignment and variant calling.
     - <b>fa</b>: Absolute path to fasta file (.fasta).
     - <b>gtf</b>: Absolute path to gene annotation file (.gtf).
     - <b>build</b>: Name of reference genome--output files will use this prefix.
     - <b>known_sites</b>: Absolute path to known variants file (.vcf.gz).

<h3>Step 4: Execute workflow</h3>
<h4>Local</h4>
Orpheus may be run on the HPC via a system-wide installation of Conda. At 
Agios, Conda may be loaded as a module then used to activate the virtual 
environment with Snakemake pre-installed:

```
module load conda
conda activate snakemake
```

Orpheus is run via the Snakemake command line interface. By default, Orpheus is
set to run using the toy data located in <code>resources</code>. Users can 
set the path to the config file to run their analysis as observed below:
```
git clone https://git.agios.local/Jeff.Alvarez/orpheus.git
cd /path/to/orpheus
snakemake --use-conda --cores 8 --configfile /path/to/configfile.yaml
```

<h4>[TO-DO: Update Docker]</h4>
The workflow may also be deployed as a Docker image, where a conda environment
is set up with Snakemake and dependencies installed. When run, the conda
environment is activated then the Snakemake directory, input FASTQ files,
and reference data are mounted as <code>analysis</code>, <code>input</code>, and <code>ref</code> respectively:

```
docker run -it --rm \
    --user "$(id -u):$(id -g)" \
    -v /data/exploratory/Users/jeff.alvarez/orpheus:/home/user/analysis \
    -v /data/exploratory/Users/jeff.alvarez/orpheus/data/samples/single:/home/user/input \
    -v /data/exploratory/Users/jeff.alvarez/orpheus/data:/home/user/ref \
    jeff.alvarez/orpheus:latest /bin/bash -c \
    "conda run -n snakemake \
    snakemake -j 8 --keep-remote --use-conda \
    --directory /home/user/analysis \
    --configfile /home/user/analysis/config/config.docker.yaml \
    -s /home/user/analysis/Snakefile"
```