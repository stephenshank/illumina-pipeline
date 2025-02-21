# Bioinformatics pipeline for Illumina viral deep sequencing

WARNING: This repository is a work in progress.

This `README.md` is intended to be a quickstart overview. For a deeper understanding of this pipeline, please see [our documentation](https://github.com/moncla-lab/illumina-pipeline/blob/main/DOCUMENTATION.md).

## Usage

Usage instructions assume you've successfully followed the installation instructions, read about and adhere to the conventions used by this software, and have a basic understanding of [conda](https://docs.conda.io/en/latest/) and [Snakemake](https://snakemake.readthedocs.io/en/stable/).

### Setup
Suppose you have several FASTQs you'd like to analyze for a project that you'll identify as `MyAnalysis`.

Get a copy of the code and set up your environment for analysis:

```
git clone https://github.com/monclalab/illumina-pipeline MyAnalysis
cd MyAnalysis
conda activate mlip
```

Make appropriate edits to `config.yml`. This will most likely `subtype` and `data_root_directory`... be sure you understand conventions and assumptions made about the latter.

### Situate data

Organize output from Illumina into the expected directory structure of this pipeline:

Newline delimited BaseSpace IDs are placed into `/path/to/basespaceIDs.txt` (this can also be a CSV or TSV, though there won't be any commas or tabs... just the IDs from BaseSpace metadata spreadsheets).

```
python mlip/dataflow.py preprocess -f /path/to/basespaceIDs.txt
python mlip/dataflow.py flow
```

### Run the pipeline

```
snakemake -j $NUMBER_OF_JOBS all
```

After this, `data` should be filled with lots of files of various formats, many which contain relevant virological information.

### Visualize outputs

```
python mlip/visualization.py
```

All paths below are assumed to be relative to the `data` directory. Relevant outputs include:
| File description     | File path     |
|--------------|--------------|
| Consensus sequences for a given sample | `{sample}/sequences.fasta` |
| Intrahost variant plot for multiple replicates | `{sample}/ml.html` |
| Project wide overview of coverage | `coverage-summary.tsv` |
| Specific sample coverage | `{sample}/replicate-{replicate}/{mapping_stage}/coverage-summary.tsv` | 

The `mapping_stage` [wildcard](https://snakemake.readthedocs.io/en/stable/snakefiles/rules.html#wildcards) is either `initial` or `remapping`.
## Installation

Requires [Bioconda](https://bioconda.github.io/) and [Git](https://git-scm.com/). Recommends [Miniconda](https://docs.anaconda.com/miniconda/).

```
git clone https://github.com/monclalab/illumina-pipeline
cd illumina-pipeline
conda create -n mlip python=3.12 pandas=2 altair biopython bedtools bcftools bowtie2 multiqc samtools trimmomatic snakemake=8.27 snpeff varscan entrez-direct seqkit sed csvkit
```

## Conventions

This repo is written for data that is fresh off of Illumina sequencers. It is assumed that the directory pointed to by `config['data_root_directory']` contains FASTQs for a particular project, as well as a metadata spreadsheet that was generated for that project. Note that not every FASTQ in the metadata spreadsheet need be present, though they certainly can be. The typical structure provided by Illumina is fine, and what this pipeline builds off of.

Moreover, it assumes that Illumina identifers that end in `_Rep1` and `_Rep2` denote replicates for the same sample. For example, `ABC_Rep1` and `ABC_Rep2` denote two technical replicates for a biological sample with ID `ABC`, whereas `XYZ` denotes one technical replicate for a biological sample with ID `XYZ`.

The first step of the pipeline (`dataflow`) will try to reconcile your metadata spreadsheet, what you've downloaded to a particular directory, and the sample ID conventions described above to situate your data for analysis.

That is the overall gist, though there is more detail in the (forthcoming) Wiki.
