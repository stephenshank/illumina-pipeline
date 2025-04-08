# Bioinformatics pipeline for Illumina viral deep sequencing

WARNING: This repository is a work in progress.

This `README.md` is intended to be a quickstart overview. For a deeper understanding of this pipeline, please see [our documentation](./DOCUMENTATION.md).

## Usage

Usage instructions assume that you've successfully followed the [installation instructions](#installation), and read about and adhere to the [conventions](#conventions) used by this software. Further, it assumes a basic understanding of command line interfaces, as well as [conda](https://docs.conda.io/en/latest/) and [Snakemake](https://snakemake.readthedocs.io/en/stable/).

### Quick start
Suppose you have several FASTQs downloaded from BaseSpace that you'd like to analyze for a project called `MyAnalysis`.

Get a copy of the code and set up your environment for analysis:

```
git clone https://github.com/monclalab/illumina-pipeline MyAnalysis
cd MyAnalysis
conda activate mlip
```

 Copy the file `config.yml.template` to `config.yml`. Make appropriate edits, which will likely involve adjusting only `virus` (what's used as the reference sequence) and the `data_root_directory` (where data from BaseSpace was downloaded) for a first run.

### Situate data

To move data out of BaseSpace's download folder and into this pipeline, a text file of BaseSpace/sequencing experiment IDs, one per line, must be created at a known path we'll call `/path/to/basespaceIDs.txt`. Each ID corresponds to exactly one FASTQ dataset (which is a pair of FASTQ files, as only paired-end reads are supported at present).

The first step looks at these IDs and builds a metadata spreadsheet at `data/metadata.tsv`:
```
python mlip/dataflow.py preprocess -f /path/to/basespaceIDs.txt
```

The pipeline will do its best to assign sample IDs to each sequencing experiment ID. The user should inspect that these sample IDs were correctly assigned, assign any that are missing, and assign each replicate.

Once the metadata spreadsheet is fully populated, the following command moves data into the pipeline in an organized manner for further analysis.

```
python mlip/dataflow.py flow
```

### Run the pipeline

With data situated, the pipeline can be ran as:
```
snakemake -j $NUMBER_OF_JOBS all
```

After this, the `data` directory should be filled with lots of files of various formats, many which contain relevant virological information.

### Visualize outputs

```
python mlip/visualization.py
```

All paths below are assumed to be relative to the `data` directory. Element enclosed in brackets are Snakemake wildcards with further explanation in documentation. Relevant outputs include:
| File description     | File path     |
|--------------|--------------|
| Consensus sequences for a given sample | `{sample}/consensus.fasta` |
| Project wide overview of coverage | `coverage-report.tsv` |
| Replicate specific coverage | `{sample}/replicate-{replicate}/{mapping_stage}/coverage.tsv` | 

## Installation

Requires [Bioconda](https://bioconda.github.io/) and [Git](https://git-scm.com/). We recommend [Miniconda](https://docs.anaconda.com/miniconda/) be used as your conda distribution.

```
git clone https://github.com/moncla-lab/illumina-pipeline
cd illumina-pipeline
conda create -n mlip python=3.12 pandas=2 altair biopython bedtools bcftools bowtie2 multiqc samtools trimmomatic snakemake=8.27 snpeff varscan entrez-direct seqkit sed csvkit perbase
```

## Conventions

This pipeline makes some assumptions about how data was generated in BaseSpace that can be overridden. In particular, we assume BaseSpace/sequencing experiment IDs are of the following format:

```
{SAMPLE}_Seq{#}
```
For example:

```
bv_w1_Seq1
bv_w1_Seq2
rf_Seq1
rf_Seq2
rf_Seq3
```

This gives an overall gist of the pipeline. For further explanation, please consult [our documentation](./DOCUMENTATION.md)
