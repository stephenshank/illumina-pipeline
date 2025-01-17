# Moncla Lab Illumina Pipeline Documentation

This document covers a more comprehensive overview of this pipeline.

## Table of Contents

## Table of Contents

1. [**Starting from reads**](#1-starting-from-reads) - Going from data fresh off of Illumina sequences to biological information
2. [**Tools utilized**](#2-tools-utilized) - Bioinformatics tools used in this pipeline with a brief description and rationale.
3. [**User stories**](#3-user-stories) - Various use cases for this pipeline and how to execute them.
4. [**Output produced**](#4-output-produced) - The various files that are produced by this pipeline.

## 1. Starting from reads

At this stage, we assume the lab technicians have done their part of the process and sequenced a few samples. There's raw read data available in Illumina Basespace waiting to be analyzed. It's now on the informatician to continue the work, and extract some information from this data. Begin by logging into BaseSpace and sure that you're associated to the correct group:

![PNG Image](documentation/accessing-different-groups.png)

Navigate to projects and select the appropriate one...

![PNG Image](documentation/navigate-to-projects.png)

Navigate to FASTQs, and select the IDs that you intend to analyze...

![PNG Image](documentation/selecting-fastqs-for-download.png)

Go back up to the file icon and select dataset from the following dropdown...

![PNG Image](documentation/download-menu-dropdown.png)

You'll be prompted with the files that you're going to download. You will need Illumina Basespace downloader to be installed on your local machine if you haven't already. Also, one quirk I've found is that this program needs to be open for me. The web application can start a download in the program if it's already open, but if it's not, it can't open it for you.

Click download to proceed...

![PNG Image](documentation/downloading-from-basespace.png)

Select the folder where you'd like your sequences locally, then click start download...

![PNG Image](documentation/selecting-local-folder.png)

If you've made it this far, and see the download in action as below, congratulations! You have completed the first step of informatics: get the data!

![PNG Image](documentation/download-in-action.png)

You can now refer back to the [README](https://github.com/moncla-lab/illumina-pipeline/blob/main/README.md) for how to set up and do a first run. Just remember the directory that you used above, as that will need to go into the configuration file for the repository.

## 2. Tools utilized

TODO: give a better description of the overall flow, but at this point it's [Snakemake](https://snakemake.readthedocs.io/en/stable/), [Bioconda](https://bioconda.github.io/), [Genbank records](https://www.ncbi.nlm.nih.gov/genbank/samplerecord/), [GTF format](https://genome.ucsc.edu/FAQ/FAQformat.html#format4), [Trimmomatic](https://academic.oup.com/bioinformatics/article/30/15/2114/2390096), [Bowtie2](https://www.nature.com/articles/nmeth.1923), [SAMtools](https://academic.oup.com/bioinformatics/article/25/16/2078/204688), [VarScan](https://genome.cshlp.org/content/22/3/568.short), [BEDTools](https://academic.oup.com/bioinformatics/article/26/6/841/244688),[SeqKit](https://journals.plos.org/plosone/article?id=10.1371/journal.pone.0163962), [MultiQC](https://academic.oup.com/bioinformatics/article/32/19/3047/2196507), [NCBI E-utilities](https://www.ncbi.nlm.nih.gov/books/NBK25501/), and [Bash](https://www.gnu.org/software/bash/).

## 3. User stories

TODO: write these up in more detail and confirm acceptability with Louise. Right now there's some functionality for individual samples:

![Individual pipeline](documentation/individual-pipeline.svg)

as well as 

![Duplicate pipeline](documentation/duplicate-pipeline.svg)

## 4. Output produced

TODO: write this up along with a description of each format