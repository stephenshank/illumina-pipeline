import os
import csv
import json
from itertools import product

import pandas as pd

from mlip import *


wildcard_constraints:
  segment="[^/]+",
  sample="[^/]+"

configfile: "config.yml"

REFERENCE = config['reference']
reference_dictionary = load_reference_dictionary(REFERENCE)
metadata_dictionary = load_metadata_dictionary()
SEGMENTS = reference_dictionary.keys()

rule fetch_reference_data:
    message:
        'Fetching reference data for segment {wildcards.segment}...'
    output:
        fasta='data/reference/{segment}/sequence.fasta',
        genbank='data/reference/{segment}/metadata.gb'
    resources:
        ncbi_fetches=1
    params:
        genbank_accession=(
            lambda wildcards:
            reference_dictionary[wildcards.segment]['genbank_accession']
        ),
    shell:
        '''
            efetch -db nuccore \
                -id {params.genbank_accession} \
                -format genbank \
                > {output.genbank}

            efetch -db nuccore \
                -id {params.genbank_accession} \
                -format fasta \
            | seqkit replace -p "^(.+)" -r "{wildcards.segment} genbank"\
                > {output.fasta}
        '''

rule build_genbank_reference:
    message:
        'Concatenating reference data into single FASTA...'
    input:
        expand(
            'data/reference/{segment}/sequence.fasta',
            segment=SEGMENTS
        )
    output:
        'data/reference/sequences.fasta',
    shell:
        'cat {input} > {output}'

rule genbank_to_gtf:
    message:
        'Converting Genbank data to GTF...'
    input:
        rules.fetch_reference_data.output.genbank
    output:
        'data/reference/{segment}/metadata.gtf'
    run:
        genbank_to_gtf(input[0], output[0])

rule full_gtf:
    input:
        expand('data/reference/{segment}/metadata.gtf', segment=SEGMENTS)
    output:
        'data/reference/metadata.gtf'
    shell:
        'cat {input} > {output}'

rule coding_regions:
    input:
        rules.full_gtf.output[0]
    output:
        'data/reference/coding_regions.json'
    run:
        coding_regions = define_coding_regions(input[0])
        gb_to_segkey = {
            v['genbank_accession']: k
            for k, v in reference_dictionary.items()
        }
        with open(output[0], 'w') as json_file:
            json.dump({
                gb_to_segkey[k]: v
                for k, v  in coding_regions.items()
            }, json_file, indent=2)



def forward_fastq_merge_inputs(wildcards):
    experiments = metadata_dictionary[wildcards.sample][wildcards.replicate]
    forward_path = 'data/%s/sequencing-{sequencing}/forward.fastq.gz' % wildcards.sample
    return expand(
        forward_path,
        sequencing=metadata_dictionary[wildcards.sample]
    )


def reverse_fastq_merge_inputs(wildcards):
    experiments = metadata_dictionary[wildcards.sample][wildcards.replicate]
    reverse_path = 'data/%s/sequencing-{sequencing}/reverse.fastq.gz' % wildcards.sample
    return expand(
        reverse_path,
        sequencing=metadata_dictionary[wildcards.sample]
    )


rule concatenate_experiments_to_replicates:
    input:
        forward=forward_fastq_merge_inputs,
        # reverse is reserved for internal Snakemake use! So we append an _
        reverse_=reverse_fastq_merge_inputs
    output:
        forward='data/{sample}/replicate-{replicate}/forward.fastq.gz',
        reverse_='data/{sample}/replicate-{replicate}/reverse.fastq.gz'
    shell:
        '''
            gzip -dc {input.forward} | gzip > {output.forward}
            gzip -dc {input.reverse_} | gzip > {output.reverse_}
        '''

rule trimmomatic:
    message:
        '''
            Trimming replicate {wildcards.replicate} of sample {wildcards.sample}...
            Parameters:
                Window size: {params.rule_trimmomatic_window_size}
                Q-score: {params.rule_trimmomatic_trim_qscore}
                Minimum length: {params.rule_trimmomatic_min_length}
        '''
    input:
        forward='data/{sample}/replicate-{replicate}/forward.fastq.gz',
        reverse_='data/{sample}/replicate-{replicate}/reverse.fastq.gz'
    output:
        forward_paired='data/{sample}/replicate-{replicate}/forward_paired.fastq',
        reverse_paired='data/{sample}/replicate-{replicate}/reverse_paired.fastq',
        forward_unpaired='data/{sample}/replicate-{replicate}/forward_unpaired.fastq',
        reverse_unpaired='data/{sample}/replicate-{replicate}/reverse_unpaired.fastq',
        stdout='data/{sample}/replicate-{replicate}/trimmomatic-stdout.txt',
        log='data/{sample}/replicate-{replicate}/trimmomatic.log',
    params: **config
    shell:
        '''
            trimmomatic PE \
                {input.forward} {input.reverse_} \
                {output.forward_paired} {output.forward_unpaired} \
                {output.reverse_paired} {output.reverse_unpaired} \
                SLIDINGWINDOW:{params.rule_trimmomatic_window_size}:{params.rule_trimmomatic_trim_qscore} \
                MINLEN:{params.rule_trimmomatic_min_length} \
                > {output.stdout} 2> {output.log}
        '''

rule index:
    message:
        'Indexing reference sequence...'
    input:
        'data/{sample}/sequences.fasta'
    params:
        'data/{sample}/index'
    output:
        index1='data/{sample}/index.1.bt2',
        index2='data/{sample}/index.2.bt2',
        index3='data/{sample}/index.3.bt2',
        index4='data/{sample}/index.4.bt2',
        indexrev1='data/{sample}/index.rev.1.bt2',
        indexrev2='data/{sample}/index.rev.2.bt2',
        stdout='data/{sample}/bowtie2-stdout.txt',
        stderr='data/{sample}/bowtie2-stderr.txt'
    shell:
        'bowtie2-build {input} {params} > {output.stdout} 2> {output.stderr}'


get_sample = lambda wildcards: 'reference' if wildcards.mapping_stage == 'initial' else wildcards.sample

rule mapping:
    message:
        'Mapping replicate {wildcards.replicate} of sample {wildcards.sample} to reference...'
    input:
        forward_paired=rules.trimmomatic.output.forward_paired,
        reverse_paired=rules.trimmomatic.output.reverse_paired,
        forward_unpaired=rules.trimmomatic.output.forward_unpaired,
        reverse_unpaired=rules.trimmomatic.output.reverse_unpaired,
        index=lambda wildcards: f'data/{get_sample(wildcards)}/index.1.bt2'
    params:
        index=lambda wildcards: f'data/{get_sample(wildcards)}/index',
        sensitivity=lambda wildcards: '--very-sensitive' if wildcards.mapping_stage == 'initial' else ''
    output:
        sam='data/{sample}/replicate-{replicate}/{mapping_stage}/mapped.sam',
        stdout='data/{sample}/replicate-{replicate}/{mapping_stage}/bowtie2-stdout.txt',
        stderr='data/{sample}/replicate-{replicate}/{mapping_stage}/bowtie2-stderr.txt'
    shell:
        '''
            bowtie2 {params.sensitivity} -x {params.index} \
                -1 {input.forward_paired} -2 {input.reverse_paired} \
                -U {input.forward_unpaired},{input.reverse_unpaired} \
                --local \
                -S {output.sam} \
                > {output.stdout} 2> {output.stderr}
        '''

rule samtools:
    message:
        'Running various samtools modules on {wildcards.replicate} of sample {wildcards.sample}...'
    input:
        rules.mapping.output.sam
    output:
        mapped='data/{sample}/replicate-{replicate}/{mapping_stage}/mapped.bam',
        sorted_='data/{sample}/replicate-{replicate}/{mapping_stage}/sorted.bam',
        index='data/{sample}/replicate-{replicate}/{mapping_stage}/sorted.bam.bai',
        stats='data/{sample}/replicate-{replicate}/{mapping_stage}/stats.txt',
        flagstat='data/{sample}/replicate-{replicate}/{mapping_stage}/flagstat.txt',
        depth='data/{sample}/replicate-{replicate}/{mapping_stage}/depth.txt',
        stdout='data/{sample}/replicate-{replicate}/{mapping_stage}/samtools-stdout.txt',
        stderr='data/{sample}/replicate-{replicate}/{mapping_stage}/samtools-stderr.txt'
    shell:
        '''
            samtools view -S -b {input} > {output.mapped} 2> {output.stderr}
            samtools sort {output.mapped} -o {output.sorted_} > {output.stdout} 2>> {output.stderr}
            samtools index {output.sorted_} >> {output.stdout} 2>> {output.stderr}
            samtools stats {output.sorted_} > {output.stats} 2>> {output.stderr}
            samtools flagstat {output.sorted_} > {output.flagstat} 2>> {output.stderr}
            samtools depth {output.sorted_} > {output.depth} 2>> {output.stderr}
        '''

get_reference = lambda wildcards: f'data/{get_sample(wildcards)}/sequences.fasta'

rule call_variants:
    message:
        '''
            Calling variants on replicate {wildcards.replicate} of sample {wildcards.sample}...
            Parameters:
                Mapping stage: {wildcards.mapping_stage}
                SNP frequency: {params.rule_call_variants_min_var_freq}
                Minimum coverage for variant calling: {params.rule_call_consensus_min_coverage}
                Strand filter: {params.tool_varscan_strand_filter}
                SNP quality threshold: {params.tool_varscan_min_avg_qual}
        '''
    input:
        bam=rules.samtools.output.sorted_,
        stderr=rules.samtools.output.stderr,
        reference=get_reference
    output:
        pileup= 'data/{sample}/replicate-{replicate}/{mapping_stage}/samtools.pileup',
        vcf=    'data/{sample}/replicate-{replicate}/{mapping_stage}/varscan.vcf',
        tsv=    'data/{sample}/replicate-{replicate}/{mapping_stage}/varscan.tsv',
        vcf_zip='data/{sample}/replicate-{replicate}/{mapping_stage}/varscan.vcf.gz',
        index=  'data/{sample}/replicate-{replicate}/{mapping_stage}/varscan.vcf.gz.tbi',
        stderr= 'data/{sample}/replicate-{replicate}/{mapping_stage}/varscan-stderr.txt'
    params:
        **config
    shell:
        '''
            samtools mpileup -A \
                -Q 0 \
                -d 1000000 \
                -f {input.reference} \
                {input.bam} > {output.pileup} 2>> {input.stderr}
            varscan mpileup2snp {output.pileup} \
                --min-coverage {params.rule_call_consensus_min_coverage} \
                --min-avg-qual {params.tool_varscan_min_avg_qual} \
                --min-var-freq {params.rule_call_variants_min_var_freq} \
                --strand-filter {params.tool_varscan_strand_filter} \
                --output-vcf 1 > {output.vcf} 2> {output.stderr}
            grep -v '^##' {output.vcf} > {output.tsv}
            bgzip -c {output.vcf} > {output.vcf_zip}
            tabix -p vcf {output.vcf_zip}
        '''

rule coverage:
    message:
        'Computing coverage of replicate {wildcards.replicate} for sample {wildcards.sample}...'
    input:
        rules.samtools.output.sorted_
    output:
        bg= 'data/{sample}/replicate-{replicate}/{mapping_stage}/coverage.bedGraph',
        tsv='data/{sample}/replicate-{replicate}/{mapping_stage}/coverage.tsv'
    shell:
        '''
            echo "segment\tstart\tend\tcoverage" > {output.tsv}
            bedtools genomecov -ibam {input} -bga > {output.bg}
            cat {output.bg} >> {output.tsv}
        '''

rule coverage_summary:
    message:
        'Computing coverage summary of replicate {wildcards.replicate} for sample {wildcards.sample}...'
    input:
        rules.coverage.output.tsv
    output:
        'data/{sample}/replicate-{replicate}/{mapping_stage}/coverage-summary.tsv'
    run:    
        compute_coverage_categories_io(input[0], output[0])

rule call_consensus:
    message:
        '''
            Calling consensus on replicate {wildcards.replicate} of sample {wildcards.sample}...
            Parameters:
                Mapping stage: {wildcards.mapping_stage}
                SNP frequency: {params.rule_call_consensus_min_var_freq}
                Minimum coverage for consensus calling: {params.rule_call_consensus_min_coverage}
                Strand filter: {params.tool_varscan_strand_filter}
                SNP quality threshold: {params.tool_varscan_min_avg_qual}
        '''
    input:
        reference=get_reference,
        pileup=rules.call_variants.output.pileup,
        vcf_zip=rules.call_variants.output.vcf_zip
    output:
        vcf='data/{sample}/replicate-{replicate}/{mapping_stage}/consensus.vcf',
        vcf_zip='data/{sample}/replicate-{replicate}/{mapping_stage}/consensus.vcf.gz',
        index='data/{sample}/replicate-{replicate}/{mapping_stage}/consensus.vcf.gz.tbi',
        fasta='data/{sample}/replicate-{replicate}/{mapping_stage}/consensus.fasta'
    params: **{ \
        **config, \
        'initial_mapping_stage': lambda wildcards: 'genbank' if wildcards.mapping_stage == 'initial' else 'initial' \
    }
    shell:
        '''
            varscan mpileup2cns {input.pileup} \
                --min-coverage {params.rule_call_consensus_min_coverage} \
                --min-avg-qual {params.tool_varscan_min_avg_qual} \
                --min-var-freq {params.rule_call_consensus_min_var_freq} \
                --strand-filter {params.tool_varscan_strand_filter} \
                --output-vcf 1 > {output.vcf}
            bgzip -c {output.vcf} > {output.vcf_zip}
            tabix -p vcf {output.vcf_zip}
            cat {input.reference} | bcftools consensus {output.vcf_zip} > {output.fasta}
            sed -i 's/{params.initial_mapping_stage}/{wildcards.mapping_stage}/g' {output.fasta}
        '''

rule mask_consensus:
    input:
        fasta=rules.call_consensus.output.fasta,
        coverage=rules.coverage.output.tsv
    output:
        'data/{sample}/replicate-{replicate}/{mapping_stage}/masked_consensus.fasta'
    run:
        mask_low_coverage(
            input.fasta,
            input.coverage,
            config['rule_call_consensus_min_coverage'],
            output[0]
        )

rule pluck_segment:
    input:
        rules.mask_consensus.output[0]
    output:
        'data/{sample}/replicate-{replicate}/{mapping_stage}/segments/{segment}.fasta',
    shell:
        '''
            seqkit grep -p "{wildcards.segment}" {input} | \
            sed 's/{wildcards.segment}/{wildcards.sample}/g' | \
            sed 's/remapping/{wildcards.segment}/g' > {output}
        '''

rule build_sample_reference:
    message:
        'Situating intrasample reference of sample {wildcards.sample} for remapping...'
    input:
        'data/{sample}/replicate-1/initial/consensus.fasta'
    output:
        'data/{sample}/sequences.fasta'
    shell:
        'cp {input} {output}'

rule multiqc:
    message:
        'Running Multi QC on {wildcards.replicate} of sample {wildcards.sample}...'
    input:
        rules.trimmomatic.output.log,
        rules.samtools.output.stats,
        rules.samtools.output.flagstat,
        rules.samtools.output.depth
    output:
        'data/{sample}/replicate-{replicate}/{mapping_stage}/multiqc_report.html'
    params:
        'data/{sample}/replicate-{replicate}/{mapping_stage}'
    shell:
        'multiqc -f {params} --outdir {params}'

rule annotate_varscan:
    input:
        coding_regions=rules.coding_regions.output[0],
        reference=rules.build_genbank_reference.output[0],
        varscan=rules.call_variants.output.vcf
    output:
        'data/{sample}/replicate-{replicate}/{mapping_stage}/varscan-annotated.tsv'
    run:
        with open(input.coding_regions) as json_file:
            coding_regions = json.load(json_file)
        transcripts = slice_fastas(coding_regions, input.reference)
        annotate_amino_acid_changes(
            coding_regions, transcripts, input.varscan, output[0]
        )

rule clean_varscan:
    message:
        'Cleaning varscan VCF from replicate {wildcards.replicate} of sample {wildcards.sample}...'
    input:
        rules.call_variants.output.tsv
    output:
        'data/{sample}/replicate-{replicate}/{mapping_stage}/ml.tsv'
    run:
        df = pd.read_csv(input[0], sep='\t')
        clean_varscan(df).to_csv(output[0], sep='\t', index=False)

def merge_varscan_inputs(wildcards):
    replicates = metadata_dictionary[wildcards.sample][wildcards.replicates]
    expand(
        'data/{{sample}}/replicate-{replicate}/remapping/ml.tsv',
        replicate=range(1, len(replicates) + 1)
    )

rule merge_varscan_across_replicates:
    message:
        'Merging variant calls of sample {wildcards.sample}...'
    input: merge_varscan_inputs
    output:
        'data/{sample}/ml.tsv'
    run:
        merge_varscan_io(input, output[0])

rule visualize_replicate_calls:
    message:
        'Visualizing replicate variant calls of sample {wildcards.sample}...'
    input:
        rules.merge_varscan_across_replicates.output[0]
    output:
        'data/{sample}/ml.html'
    run:
        replicate_variant_plot(input[0], output[0])


def full_coverage_summary_input(wildcards):
    coverage_filepaths = []
    for sample, replicates in metadata_dictionary.items():
        for replicate in replicates.keys():
            coverage_filepaths.append(
                f'data/{sample}/replicate-{replicate}/remapping/coverage-summary.tsv'
            )
    return coverage_filepaths


rule full_coverage_summary:
    input: full_coverage_summary_input
    output:
        'data/coverage-summary.tsv',
    run:
        coverage_summary(input, output[0])


def full_segment_input(wildcards):
    segment_filepaths = []
    for sample, replicates in metadata_dictionary.items():
        for replicate in replicates.keys():
            segment_filepaths.append(
                f'data/{sample}/replicate-{replicate}/remapping/segments/{wildcards.segment}.fasta'
            )
    return segment_filepaths

rule full_segment:
    input: full_segment_input
    output:
        'data/{segment}.fasta'
    shell:
        'cat {input} > {output}'

rule check_consensus:
    input:
        fasta=rules.call_consensus.output.fasta,
        pileup=rules.call_variants.output.pileup
    output:
        'data/{sample}/replicate-{replicate}/{mapping_stage}/consensus-report.tsv'
    run:
        check_consensus_io(input.fasta, input.pileup, output[0])

def full_consensus_summary_input(wildcards):
    consensus_filepaths = []
    for sample, replicates in metadata_dictionary.items():
        for replicate in replicates.keys():
            consensus_filepaths.append(
                f'data/{sample}/replicate-{replicate}/remapping/consensus-report.tsv'
            )
    return consensus_filepaths


rule full_consensus_summary:
    input: full_consensus_summary_input
    output:
        'data/consensus-summary.tsv',
    shell:
        'csvstack {input} > {output}'

rule all_full_segments:
    input:
        expand(
            'data/{segment}.fasta',
            segment=SEGMENTS
        )
    output:
        'data/all.fasta'
    shell:
        'cat {input} > {output}'

rule all:
    input:
        rules.full_coverage_summary.output[0]

rule clean:
    shell:
        '''
            rm $(find data -type f | grep -v forward.fastq.gz$ | grep -v reverse.fastq.gz)
            find data -type d -empty -delete
        '''