import os
import csv
import json
from itertools import product

import pandas as pd

from mlip import *


wildcard_constraints:
  segment="[^/]+",
  sample="[^/]+",
  replicate="[^/]+",
  mapping_stage="[^/]+"

configfile: "config.yml"

REFERENCE = config['reference']
reference_dictionary = load_reference_dictionary(REFERENCE)
metadata_dictionary = load_metadata_dictionary()
SAMPLES = samples_to_analyze()
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
    result = expand(
        forward_path,
        sequencing=metadata_dictionary[wildcards.sample][wildcards.replicate]
    )
    return result


def reverse_fastq_merge_inputs(wildcards):
    experiments = metadata_dictionary[wildcards.sample][wildcards.replicate]
    reverse_path = 'data/%s/sequencing-{sequencing}/reverse.fastq.gz' % wildcards.sample
    return expand(
        reverse_path,
        sequencing=metadata_dictionary[wildcards.sample][wildcards.replicate]
    )


rule concatenate_experiments_to_replicates:
    input:
        forward=forward_fastq_merge_inputs,
        # reverse is reserved for internal Snakemake use! So we append an _
        reverse_=reverse_fastq_merge_inputs
    output:
        forward='data/{sample}/{segment}/replicate-{replicate}/forward.fastq.gz',
        reverse_='data/{sample}/{segment}/replicate-{replicate}/reverse.fastq.gz'
    message:
        '''
            Concatenating
                {input.forward}
            to {output.forward}
            as well as 
                {input.reverse_}
            to {output.reverse_}...
        '''
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
                Window size: {params.trimming_window_size}
                Q-score: {params.minimum_quality_score}
                Minimum length: {params.trimming_minimum_length}
        '''
    input:
        forward='data/{sample}/{segment}/replicate-{replicate}/forward.fastq.gz',
        reverse_='data/{sample}/{segment}/replicate-{replicate}/reverse.fastq.gz'
    output:
        forward_paired='data/{sample}/{segment}/replicate-{replicate}/forward_paired.fastq',
        reverse_paired='data/{sample}/{segment}/replicate-{replicate}/reverse_paired.fastq',
        forward_unpaired='data/{sample}/{segment}/replicate-{replicate}/forward_unpaired.fastq',
        reverse_unpaired='data/{sample}/{segment}/replicate-{replicate}/reverse_unpaired.fastq',
        stdout='data/{sample}/{segment}/replicate-{replicate}/trimmomatic-stdout.txt',
        log='data/{sample}/{segment}/replicate-{replicate}/trimmomatic.log',
    params: **config
    shell:
        '''
            trimmomatic PE \
                {input.forward} {input.reverse_} \
                {output.forward_paired} {output.forward_unpaired} \
                {output.reverse_paired} {output.reverse_unpaired} \
                SLIDINGWINDOW:{params.trimming_window_size}:{params.minimum_quality_score} \
                MINLEN:{params.trimming_minimum_length} \
                ILLUMINACLIP:adapters.fasta:2:30:10 \
                > {output.stdout} 2> {output.log}
        '''


def situate_reference_input(wildcards):
    if wildcards.mapping_stage == 'initial':
        return f'data/reference/{wildcards.segment}/sequence.fasta'
    elif wildcards.mapping_stage == 'remapping':
        return f'data/{wildcards.sample}/{wildcards.segment}/replicate-{wildcards.replicate}/initial/consensus.fasta'
    else:
        return f'data/{wildcards.sample}/{wildcards.segment}/replicate-{wildcards.replicate}/remapping/consensus.fasta'


rule situate_reference:
    input:
        situate_reference_input
    output:
        'data/{sample}/{segment}/replicate-{replicate}/{mapping_stage}/reference/sequence.fasta'
    shell:
        'cp {input} {output}'

rule index:
    message:
        'Indexing reference sequence...'
    input:
        rules.situate_reference.output[0]
    params:
        'data/{sample}/{segment}/replicate-{replicate}/{mapping_stage}/reference/index'
    output:
        index1='data/{sample}/{segment}/replicate-{replicate}/{mapping_stage}/reference/index.1.bt2',
        index2='data/{sample}/{segment}/replicate-{replicate}/{mapping_stage}/reference/index.2.bt2',
        index3='data/{sample}/{segment}/replicate-{replicate}/{mapping_stage}/reference/index.3.bt2',
        index4='data/{sample}/{segment}/replicate-{replicate}/{mapping_stage}/reference/index.4.bt2',
        indexrev1='data/{sample}/{segment}/replicate-{replicate}/{mapping_stage}/reference/index.rev.1.bt2',
        indexrev2='data/{sample}/{segment}/replicate-{replicate}/{mapping_stage}/reference/index.rev.2.bt2',
        stdout='data/{sample}/{segment}/replicate-{replicate}/{mapping_stage}/reference/bowtie2-stdout.txt',
        stderr='data/{sample}/{segment}/replicate-{replicate}/{mapping_stage}/reference/bowtie2-stderr.txt'
    shell:
        'bowtie2-build {input} {params} > {output.stdout} 2> {output.stderr}'


rule mapping:
    message:
        'Mapping replicate {wildcards.replicate} of sample {wildcards.sample} to reference...'
    input:
        forward_paired=rules.trimmomatic.output.forward_paired,
        reverse_paired=rules.trimmomatic.output.reverse_paired,
        forward_unpaired=rules.trimmomatic.output.forward_unpaired,
        reverse_unpaired=rules.trimmomatic.output.reverse_unpaired,
        index=rules.index.output.index1
    params:
        'data/{sample}/{segment}/replicate-{replicate}/{mapping_stage}/reference/index'
    output:
        sam='data/{sample}/{segment}/replicate-{replicate}/{mapping_stage}/mapped.sam',
        stdout='data/{sample}/{segment}/replicate-{replicate}/{mapping_stage}/bowtie2-stdout.txt',
        stderr='data/{sample}/{segment}/replicate-{replicate}/{mapping_stage}/bowtie2-stderr.txt'
    shell:
        '''
            bowtie2 --very-sensitive-local -x {params} --local \
                -1 {input.forward_paired} -2 {input.reverse_paired} \
                -U {input.forward_unpaired},{input.reverse_unpaired} \
                -S {output.sam} > {output.stdout} 2> {output.stderr}
        '''

rule samtools:
    message:
        'Running various samtools modules on {wildcards.replicate} of sample {wildcards.sample}...'
    input:
        sam=rules.mapping.output.sam,
        reference=rules.situate_reference.output[0]
    output:
        mapped='data/{sample}/{segment}/replicate-{replicate}/{mapping_stage}/mapped.bam',
        sorted_='data/{sample}/{segment}/replicate-{replicate}/{mapping_stage}/sorted.bam',
        index='data/{sample}/{segment}/replicate-{replicate}/{mapping_stage}/sorted.bam.bai',
        pileup='data/{sample}/{segment}/replicate-{replicate}/{mapping_stage}/pileup.txt',
        stats='data/{sample}/{segment}/replicate-{replicate}/{mapping_stage}/stats.txt',
        flagstat='data/{sample}/{segment}/replicate-{replicate}/{mapping_stage}/flagstat.txt',
        depth='data/{sample}/{segment}/replicate-{replicate}/{mapping_stage}/depth.txt',
        stdout='data/{sample}/{segment}/replicate-{replicate}/{mapping_stage}/samtools-stdout.txt',
        stderr='data/{sample}/{segment}/replicate-{replicate}/{mapping_stage}/samtools-stderr.txt'
    shell:
        '''
            samtools view -S -b {input.sam} > {output.mapped} 2> {output.stderr}
            samtools sort {output.mapped} -o {output.sorted_} > {output.stdout} 2>> {output.stderr}
            samtools index {output.sorted_} >> {output.stdout} 2>> {output.stderr}
            samtools stats {output.sorted_} > {output.stats} 2>> {output.stderr}
            samtools flagstat {output.sorted_} > {output.flagstat} 2>> {output.stderr}
            samtools depth {output.sorted_} > {output.depth} 2>> {output.stderr}
            samtools mpileup -A -B -aa -Q 0 -d 0 -f {input.reference} {output.sorted_} > {output.pileup} 2>> {output.stderr}
        '''

rule call_variants:
    message:
        '''
            Calling variants on replicate {wildcards.replicate} of sample {wildcards.sample}...
            Parameters:
                Mapping stage: {wildcards.mapping_stage}
                SNP frequency: {params.variants_minimum_frequency}
                Minimum coverage for variant calling: {params.variants_minimum_coverage}
                Strand filter: {params.strand_filter}
                SNP quality threshold: {params.minimum_quality_score}
        '''
    input:
        bam=rules.samtools.output.sorted_,
        pileup=rules.samtools.output.pileup,
        stderr=rules.samtools.output.stderr,
        reference=situate_reference_input
    output:
        vcf=    'data/{sample}/{segment}/replicate-{replicate}/{mapping_stage}/varscan.vcf',
        tsv=    'data/{sample}/{segment}/replicate-{replicate}/{mapping_stage}/varscan.tsv',
        vcf_zip='data/{sample}/{segment}/replicate-{replicate}/{mapping_stage}/varscan.vcf.gz',
        index=  'data/{sample}/{segment}/replicate-{replicate}/{mapping_stage}/varscan.vcf.gz.tbi',
        stderr= 'data/{sample}/{segment}/replicate-{replicate}/{mapping_stage}/varscan-stderr.txt'
    params:
        **config
    shell:
        '''
            varscan mpileup2snp {input.pileup} \
                --min-coverage {params.variants_minimum_coverage} \
                --min-avg-qual {params.minimum_quality_score} \
                --min-var-freq {params.variants_minimum_frequency} \
                --strand-filter {params.strand_filter} \
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
        bg= 'data/{sample}/{segment}/replicate-{replicate}/{mapping_stage}/coverage.bedGraph',
        tsv='data/{sample}/{segment}/replicate-{replicate}/{mapping_stage}/coverage.tsv'
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
        'data/{sample}/{segment}/replicate-{replicate}/{mapping_stage}/coverage-summary.tsv'
    run:    
        compute_coverage_categories_io(input[0], output[0])

rule call_all_consensus:
    input:
        bam=rules.samtools.output.sorted_,
        pileup=rules.samtools.output.pileup,
        reference=situate_reference_input,
        original_reference=rules.fetch_reference_data.output.fasta
    output:
        vc_fasta='data/{sample}/{segment}/replicate-{replicate}/{mapping_stage}/vc_consensus.fasta',
        vc_counts='data/{sample}/{segment}/replicate-{replicate}/{mapping_stage}/vc_counts.tsv',
        vc_json='data/{sample}/{segment}/replicate-{replicate}/{mapping_stage}/vc_insertions.json',
        vcf='data/{sample}/{segment}/replicate-{replicate}/{mapping_stage}/varscan-consensus.vcf',
        tsv='data/{sample}/{segment}/replicate-{replicate}/{mapping_stage}/varscan-consensus.tsv',
        vcf_zip='data/{sample}/{segment}/replicate-{replicate}/{mapping_stage}/varscan-consensus.vcf.gz',
        index='data/{sample}/{segment}/replicate-{replicate}/{mapping_stage}/varscan-consensus.vcf.gz.tbi',
        varscan_fasta='data/{sample}/{segment}/replicate-{replicate}/{mapping_stage}/varscan-consensus.fasta',
        ivar_fasta='data/{sample}/{segment}/replicate-{replicate}/{mapping_stage}/ivar.fa',
        vc_unaligned='data/{sample}/{segment}/replicate-{replicate}/{mapping_stage}/vc_unaligned.fasta',
        vc_aligned='data/{sample}/{segment}/replicate-{replicate}/{mapping_stage}/vc_aligned.fasta',
        ivar_unaligned='data/{sample}/{segment}/replicate-{replicate}/{mapping_stage}/ivar_unaligned.fasta',
        ivar_aligned='data/{sample}/{segment}/replicate-{replicate}/{mapping_stage}/ivar_aligned.fasta'
    params: ** { \
        **config, \
        'ivar': 'data/{sample}/{segment}/replicate-{replicate}/{mapping_stage}/ivar' \
    }
    shell:
        '''
            set +e
            viral_consensus -i {input.bam} \
                --ref_genome {input.reference} \
                --out_consensus {output.vc_fasta} \
                --min_depth {params.consensus_minimum_coverage} \
                --min_qual {params.minimum_quality_score} \
                --out_pos_counts {output.vc_counts} \
                --out_ins_counts {output.vc_json}

            varscan mpileup2cns {input.pileup} \
                --min-coverage {params.consensus_minimum_coverage} \
                --min-avg-qual {params.minimum_quality_score} \
                --min-var-freq {params.consensus_minimum_frequency} \
                --strand-filter {params.strand_filter} \
                --output-vcf 1 > {output.vcf}
            bgzip -c {output.vcf} > {output.vcf_zip}
            tabix -p vcf {output.vcf_zip}
            cat {input.reference} | bcftools consensus {output.vcf_zip} > {output.varscan_fasta}
            grep -v '^##' {output.vcf} > {output.tsv}

            cat {input.pileup} | ivar consensus -p {params.ivar} \
                -m {params.consensus_minimum_coverage} \
                -q {params.minimum_quality_score} \
                -t {params.consensus_minimum_frequency}

            cat {input.original_reference} {output.vc_fasta} > {output.vc_unaligned}
            mafft {output.vc_unaligned} > {output.vc_aligned}
            cat {input.original_reference} {output.ivar_fasta} > {output.ivar_unaligned}
            mafft {output.ivar_unaligned} > {output.ivar_aligned}
        '''

#rule clip_consensus:
    #input:
        #vc=rules.call_all_consensus.output.vc_aligned,
        #ivar=rules.call_all_consensus.output.ivar_aligned
    #output:
        #vc='data/{sample}/{segment}/replicate-{replicate}/{mapping_stage}/vc_clipped.fasta',
        #ivar='data/{sample}/{segment}/replicate-{replicate}/{mapping_stage}/ivar_clipped.fasta'
    #run:
        #clip_consensus(input.vc, output.vc)
        #clip_consensus(input.ivar, output.ivar)

rule call_consensus:
    input:
        rules.call_all_consensus.output.ivar_fasta
    output:
        'data/{sample}/{segment}/replicate-{replicate}/{mapping_stage}/consensus.fasta'
    shell:
        '''
            echo ">{wildcards.sample} {wildcards.mapping_stage}" > {output}
            tail -n +2 {input} >> {output}
        '''

rule multiqc:
    message:
        'Running Multi QC on {wildcards.replicate} of sample {wildcards.sample}...'
    input:
        rules.trimmomatic.output.log,
        rules.samtools.output.stats,
        rules.samtools.output.flagstat,
        rules.samtools.output.depth
    output:
        'data/{sample}/{segment}/replicate-{replicate}/{mapping_stage}/multiqc_report.html'
    params:
        'data/{sample}/{segment}/replicate-{replicate}/{mapping_stage}'
    shell:
        'multiqc -f {params} --outdir {params}'

rule annotate_varscan:
    input:
        coding_regions=rules.coding_regions.output[0],
        reference=rules.build_genbank_reference.output[0],
        varscan=rules.call_variants.output.vcf
    output:
        'data/{sample}/{segment}/replicate-{replicate}/{mapping_stage}/varscan-annotated.tsv'
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
        'data/{sample}/{segment}/replicate-{replicate}/{mapping_stage}/ml.tsv'
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
                f'data/{sample}/{segment}/replicate-{replicate}/remapping/coverage-summary.tsv'
            )
    return coverage_filepaths


rule full_coverage_summary:
    input: full_coverage_summary_input
    output:
        'data/coverage-summary.tsv',
    run:
        coverage_summary(input, output[0])

rule compare_all_consensus:
    input:
        reference=rules.fetch_reference_data.output.fasta,
        vc_initial='data/{sample}/{segment}/replicate-{replicate}/initial/vc_consensus.fasta',
        ivar_initial='data/{sample}/{segment}/replicate-{replicate}/initial/ivar.fa',
        vc_remapping='data/{sample}/{segment}/replicate-{replicate}/remapping/vc_consensus.fasta',
        ivar_remapping='data/{sample}/{segment}/replicate-{replicate}/remapping/ivar.fa',
        vc_reremapping='data/{sample}/{segment}/replicate-{replicate}/reremapping/vc_consensus.fasta',
        ivar_reremapping='data/{sample}/{segment}/replicate-{replicate}/reremapping/ivar.fa',
    output:
        unaligned='data/{sample}/{segment}/replicate-{replicate}/all_consensus_unaligned.fasta',
        aligned='data/{sample}/{segment}/replicate-{replicate}/all_consensus_aligned.fasta'
    shell:
        'cat {input} > {output.unaligned}; mafft {output.unaligned} > {output.aligned}'


def sample_consensus_input(wildcards):
    replicates = metadata_dictionary[wildcards.sample]
    return expand(
        'data/{{sample}}/{{segment}}/replicate-{replicate}/reremapping/consensus.fasta',
        replicate=replicates
    )


rule sample_consensus:
    input:
        sample_consensus_input
    output:
        'data/{sample}/{segment}/consensus.fasta'
    run:
        call_sample_consensus(input, output[0])

def all_segments_input(wildcards):
    segment_filepaths = []
    for sample, replicates in metadata_dictionary.items():
        for replicate in replicates.keys():
            segment_filepaths.append(
                f'data/{sample}/{segment}/replicate-{replicate}/reremapping/consensus.fasta'
            )
    return segment_filepaths

rule all_segments:
    input: all_segments_input
    output:
        'data/{segment}.fasta'
    shell:
        'cat {input} > {output}'

rule check_consensus:
    input:
        fasta=rules.call_consensus.output[0],
        pileup=rules.samtools.output.pileup,
        alignment=rules.compare_all_consensus.output.aligned
    output:
        'data/{sample}/{segment}/replicate-{replicate}/{mapping_stage}/consensus-report.tsv'
    run:
        check_consensus_io(input.fasta, input.pileup, output[0], wildcards.sample, wildcards.replicate)


def full_consensus_summary_input(wildcards):
    consensus_filepaths = []
    replicates = metadata_dictionary[wildcards.sample]
    for replicate in replicates.keys():
        for segment in SEGMENTS:
            consensus_filepaths.append(
                f'data/{wildcards.sample}/{segment}/replicate-{replicate}/reremapping/consensus-report.tsv'
            )
    return consensus_filepaths


rule check_sample_consensus:
    input:
        full_consensus_summary_input
    output:
        'data/{sample}/consensus-report.tsv'
    shell:
        'csvstack {input} > {output}'

rule full_consensus_summary:
    input:
        expand('data/{sample}/consensus-report.tsv', sample=SAMPLES)
    output:
        'data/consensus-report.tsv',
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

rule all_protein:
    input:
        segments=expand('data/{segment}.fasta', segment=SEGMENTS),
        genbanks=expand('data/reference/{segment}/metadata.gb', segment=SEGMENTS)
    output:
        directory('data/protein')
    run:
        translate_all(input.segments, input.genbanks, output[0], SEGMENTS)

rule all_preliminary:
    input:
        rules.full_coverage_summary.output[0]

rule all:
    input:
        rules.full_consensus_summary.output[0]

rule zip:
    input:
        rules.full_coverage_summary.output[0],
        rules.full_consensus_summary.output[0]
    output:
        'data/project.zip'
    shell:
        'zip {output} $(find data -type f | grep -v -e fastq -e bam -e sam -e pileup)'