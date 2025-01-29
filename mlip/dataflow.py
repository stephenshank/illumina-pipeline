import os
import csv
import json
import shutil
import re
from pathlib import Path
import argparse

from Bio import SeqIO
from Bio.Seq import Seq
import pandas as pd
import yaml


def load_config():
    with open('config.yml') as config_file:
        config = yaml.safe_load(config_file)
    return config


def load_reference_dictionary(subtype):
    reference_dictionary = {}
    with open('references.tsv', 'r', newline='') as tsv_file:
        reader = csv.DictReader(tsv_file, delimiter='\t')
        for row in reader:
            if row['subtype'] == subtype:
                reference_dictionary[row['segment_key']] = row
    return reference_dictionary


def tokenize(identifier):
    return identifier.lower().replace('_', '').replace('-', '')


def preprocess(id_filepath, seq_key='Seq'):
    config = load_config()
    with open(id_filepath, 'r') as f:
        lines = f.read().splitlines()
    raw_keys = set(lines)
    key_hash = {}

    seq_key_pattern = re.compile(rf'_{seq_key}(\d+)')
    # process keys
    for raw_key in raw_keys:
        match = seq_key_pattern.search(raw_key)
        found_match = False
        if match:
            sample_key = raw_key.replace(match.group(0), '').lower()
            sequence_key = f"sequence-{match.group(1)}"
            found_match = True
        else:
            sample_key = raw_key.lower()
            sequence_key = 'sequence-1'
        new_value = {
            'original_key': raw_key,
            'found_match': found_match,
            'forward_count': 0,
            'forward_filepath': '',
            'reverse_count': 0,
            'reverse_filepath': ''
        }
        if sample_key in key_hash:
            key_hash[sample_key][sequence_key] = new_value
        else:
            key_hash[sample_key]= {sequence_key: new_value}

    os.makedirs('data', exist_ok=True)
    rows = []
    for sample_key, sample_value in key_hash.items():
        for sequence_key, sequence_value in sample_value.items():
            rows.append({
                'SequencingId': sequence_value['original_key'],
                'SampleID': sample_key if sequence_value['found_match'] else '',
                'Replicate': ''
            })

    sorted_rows = sorted(rows, key=lambda x: x['SequencingId'].lower())
    print(sorted_rows)
    with open("data/metadata.tsv", "w", newline="") as f:
        writer = csv.DictWriter(f, fieldnames=sorted_rows[0].keys(), delimiter="\t")
        writer.writeheader()
        for row in sorted_rows:
            writer.writerow(row)

    return
    # match to files
    for sample_key in key_hash.keys():
        run_key_token = tokenize(run_key)
        token_pattern = re.compile(rf'^{run_key_token}(?!\d)')
        fastq_paths = Path(config['data_root_directory']).expanduser().rglob('*.fastq.gz')
        for fastq_path in fastq_paths:
            fastq_basename = os.path.basename(fastq_path)
            fastq_token = tokenize(fastq_basename)
            if token_pattern.search(fastq_token):
                if 'rep1' in fastq_token:
                    replicate_key = 'replicate-1'
                elif 'rep2' in fastq_token:
                    replicate_key = 'replicate-2'
                else:
                    replicate_key = 'replicate-1'
                is_forward = fastq_basename.split('_')[-2] == 'R1'
                direction_key = 'forward' if is_forward else 'reverse'    
                count_key = f'{direction_key}_count'
                key_hash[run_key][replicate_key][count_key] += 1
                key_hash[run_key][replicate_key][f'{direction_key}_filepath'] = fastq_path
                if key_hash[run_key][replicate_key][count_key] > 1:
                    raise ValueError("Multiple FASTQs found for same token... aborting!!")
    
    # clean up
    keys_to_delete = []
    for key, value in key_hash.items():
        if value['replicate-1']['forward_count'] == 0 and value['replicate-1']['reverse_count'] == 0:
            keys_to_delete.append(key)
    for key in keys_to_delete:
        del key_hash[key]
    
    # move or copy
    for sample_key, sample_value in key_hash.items():
        for replicate_key, replicate_value in sample_value.items():
            action = 'mov'
            old_forward_path = replicate_value['forward_filepath']
            old_forward_name = os.path.basename(old_forward_path)
            new_forward_path = os.path.join('data', sample_key, replicate_key, 'forward.fastq.gz')
            print(f'{action}ing {old_forward_name} to {new_forward_path}')
            directory_path = os.path.dirname(new_forward_path)
            os.makedirs(directory_path, exist_ok=True)
            shutil.copy(old_forward_path, new_forward_path)

            old_reverse_path = replicate_value['reverse_filepath']
            old_reverse_name = os.path.basename(old_reverse_path)
            new_reverse_path = os.path.join('data', sample_key, replicate_key, 'reverse.fastq.gz')
            print(f'{action}ing {old_reverse_name} to {new_reverse_path}\n')
            directory_path = os.path.dirname(new_forward_path)
            os.makedirs(directory_path, exist_ok=True)
            shutil.copy(old_reverse_path, new_reverse_path)

    sample_keys_filepath = os.path.join('data', 'samples.txt')
    with open(sample_keys_filepath, 'w') as samples_file:
        samples_file.write('\n'.join(key_hash.keys()))


def preprocess_cli(args):
    if args.key:
        preprocess(args.file, args.key)
    else:
        preprocess(args.file)


def flow():
    print('flow')


def flow_cli():
    flow()


def command_line_interface():
    parser = argparse.ArgumentParser(description="MLIP dataflow")
    subparsers = parser.add_subparsers(dest="command", required=True)

    pre_parser = subparsers.add_parser("preprocess", help="Convert Sequencing ID list to metadata spreadsheet.")
    pre_parser.add_argument("-f", "--file", required=True, type=str, help="Path to the ID list")
    pre_parser.add_argument("-k", "--key", required=False, type=str, help="Key to denote sequencing experiment")
    pre_parser.set_defaults(func=preprocess_cli)

    flow_parser = subparsers.add_parser("flow", help="Move data based on metadata spreadsheet.")
    flow_parser.add_argument("-f", "--file", required=True, type=str, help="Path to the metadata spreadsheet.")
    flow_parser.set_defaults(func=flow_cli)

    args = parser.parse_args()
    args.func(args)


def genbank_to_gtf(gbk_file, gtf_file):
    with open(gtf_file, 'w') as gtf_out:
        for record in SeqIO.parse(gbk_file, "genbank"):
            for feature in record.features:
                if feature.type == "CDS":
                    # Extract required fields
                    chrom = record.id
                    source = "GenBank"
                    feature_type = "CDS"
                    start = int(feature.location.start) + 1  # GTF uses 1-based indexing
                    end = int(feature.location.end)  # Already 1-based in GenBank
                    score = "."
                    strand = "+" if feature.location.strand == 1 else "-"
                    frame = str(feature.qualifiers.get('codon_start', ['0'])[0])
                    gene_name = feature.qualifiers.get('gene', ['unknown'])[0]
                    transcript_id = feature.qualifiers.get('protein_id', ['unknown'])[0]

                    # Format the GTF line
                    gtf_line = f"{chrom}\t{source}\t{feature_type}\t{start}\t{end}\t{score}\t{strand}\t{frame}\tgene_id \"{gene_name}\"; transcript_id \"{transcript_id}\";\n"
                    gtf_out.write(gtf_line)


def clean_varscan(varscan_df):
    varscan_df['Frequency'] = [
        float(format_.split(':')[6][:-1]) / 100
        for format_ in varscan_df['Sample1']
    ]
    return varscan_df.loc[:, ['#CHROM', 'POS', 'REF', 'ALT', 'Frequency']]


def merge_varscan(clean_varscan_dfs):
    variant_dictionary = {}
    for idx, clean_varscan_df in enumerate(clean_varscan_dfs):
        sample_key = 'sample_%d' % (idx+1)
        for _, row in clean_varscan_df.iterrows():
            variant_key = (row['#CHROM'], row['POS'], row['ALT'])
            variant_data = {'Frequency': row['Frequency']}
            if variant_key in variant_dictionary:
                variant_dictionary[variant_key][sample_key] = variant_data
            else:
                variant_dictionary[variant_key] = {sample_key: variant_data}
    
    variant_list = []
    for variant_key, variant_value in variant_dictionary.items():
        for idx in range(len(clean_varscan_dfs)):
            sample_key = 'sample_%d' % (idx+1)
            if not sample_key in variant_value:
                variant_value[sample_key] = {'Frequency': 0}
        
        flattened = {}
        for idx, sample_attributes in enumerate(variant_value.values()):
            for attribute_key, attribute_value in sample_attributes.items():
                flattened[f"{attribute_key}_{idx+1}"] = attribute_value
        variant_list.append({
            'segment': variant_key[0],
            'position': variant_key[1],
            'alt': variant_key[2],
            **flattened
        })
    return pd.DataFrame(variant_list)


def merge_varscan_io(input_tsv_filepaths, output_tsv_filepath):
        dfs = [
            pd.read_csv(input_tsv_filepath, sep='\t')
            for input_tsv_filepath in input_tsv_filepaths
        ]
        merged_df = merge_varscan(dfs)
        merged_df.to_csv(output_tsv_filepath, sep='\t')


def assign_coverage_bucket(coverage):
    if coverage == 0:
        return '0x'
    elif 1 <= coverage <= 100:
        return '1-100x'
    elif 101 <= coverage <= 1000:
        return '100-1000x'
    else:
        return '1000x+'


def compute_coverage_categories(df):
    # Calculate the number of sites in each range
    df['num_sites'] = df['end'] - df['start']
    
    # Assign each coverage value to a bucket
    df['coverage_bucket'] = df['coverage'].apply(assign_coverage_bucket)
    
    # Group by segment and coverage bucket and sum the number of sites
    bucket_summary = df.groupby(
        ['segment', 'coverage_bucket']
    )['num_sites'].sum().reset_index(name='count')
    
    # Create all possible combinations of segments and coverage buckets
    all_buckets = pd.DataFrame({
        'coverage_bucket': ['0x', '1-100x', '100-1000x', '1000x+']
    })
    segments = df['segment'].unique()
    all_combinations = pd.MultiIndex.from_product(
        [segments, all_buckets['coverage_bucket']],
        names=['segment', 'coverage_bucket']
    ).to_frame(index=False)
    
    # Merge the combinations with the actual bucket summary and fill missing values with 0
    bucket_summary = pd.merge(
        all_combinations,
        bucket_summary,
        on=['segment', 'coverage_bucket'],
        how='left'
    ).fillna(0)
    
    # Pivot to the desired format
    desired_structure = bucket_summary.pivot_table(
        index='segment', 
        columns='coverage_bucket', 
        values='count', 
        fill_value=0
    ).reset_index()
   
    desired_structure.columns.name = None
    return desired_structure

def compute_coverage_categories_io(input_coverage, output_summary):
    coverage_df = pd.read_csv(input_coverage, sep='\t')
    coverage_summary_df = compute_coverage_categories(coverage_df)
    coverage_summary_df.to_csv(output_summary, sep='\t', index=False)


def define_coding_regions(gtf_file):
    with open(gtf_file, "r") as gtf:
        coding_regions = {}

        for line in gtf:
            if line.strip("\n") != "":      # ignore blank lines (otherwise throws an index error)
                sequence_name = line.split("\t")[0]
                annotation_type = line.split("\t")[2]
                start = int(line.split("\t")[3]) - 1  # adding the -1 here for 0 indexing
                stop = int(line.split("\t")[4]) - 1    # adding the -1 here for 0 indexing
                gene_name = line.split("\t")[8]
                gene_name = gene_name.split(";")[0]
                gene_name = gene_name.replace("gene_id ","")
                gene_name = gene_name.replace("\"","")

                if annotation_type.lower() == "cds":
                    if sequence_name not in coding_regions:
                        coding_regions[sequence_name] = {}
                        coding_regions[sequence_name][gene_name] = [start, stop]
                    elif sequence_name in coding_regions and gene_name not in coding_regions[sequence_name]:
                        coding_regions[sequence_name][gene_name] = [start, stop]
                    elif gene_name in coding_regions[sequence_name]:
                        coding_regions[sequence_name][gene_name].extend([start, stop])

            # sort coding region coordinates so that they are always in the correct order
            for sequence_name in coding_regions:
                for gene in coding_regions[sequence_name]:
                    coding_regions[sequence_name][gene] = sorted(coding_regions[sequence_name][gene])

    return(coding_regions)


amino_acid_abbreviations = {
    "A": "Ala", "R": "Arg", "N": "Asn", "D": "Asp", "C": "Cys",
    "Q": "Gln", "E": "Glu", "G": "Gly", "H": "His", "I": "Ile",
    "L": "Leu", "K": "Lys", "M": "Met", "F": "Phe", "P": "Pro",
    "O": "Pyl", "S": "Ser", "U": "Sec", "T": "Thr", "W": "Trp",
    "Y": "Tyr", "V": "Val", "B": "Asx", "Z": "Glx", "U": "Sec",
    "X": "Xaa", "J": "Xle", "*": "Stop"
}


def slice_fastas(coding_regions, reference_sequence_path):
    reference_sequence = {}

    for seq in SeqIO.parse(reference_sequence_path, "fasta"):
        sequence_name = str(seq.id)
        sequence = str(seq.seq).lower()
        reference_sequence[sequence_name] = sequence

    transcripts = {}
    stop_codons = ["taa","tag","tga"]
    refseqname = 'A sequence'

    for c in coding_regions:
        for gene in coding_regions[c]:
            transcripts[gene] = ""
            coordinates = coding_regions[c][gene]   # define the coding regions for each gene

            for i in range(0,int(len(coordinates)),2):   # go along the coordinates in chunks of 2 at a time
                sequence_chunk = reference_sequence[c][coordinates[i]:coordinates[i+1]+1]
                transcripts[gene] += sequence_chunk     # append each piece of the transcript together 

    # loop through each transcript to make sure that it begins with a start codon and ends with a stop codon
    for t in transcripts:
        if transcripts[t][0:3] != "atg":
            print("WARNING! " + refseqname + " " + t + "does not contain a start codon!")
        if transcripts[t][-3:] not in stop_codons:
            print("WARNING! " + refseqname + " " + t + " does not contain a stop codon! These are the last 3 nucleotides: " + transcripts[t][-3:])

    return transcripts


def annotate_amino_acid_changes(coding_regions, transcripts, vcf, outfilename):
    with open(vcf, "r") as csvfile:
        with open(outfilename, "w") as outfile:
            to_write = ["sample","gene","reference_position","reference_allele","variant_allele","coding_region_change","synonymous/nonsynonymous","frequency(%)","frequency","\n"]
            to_write2 = "\t".join(to_write)
            outfile.write(to_write2)

        reader = csv.reader(csvfile, delimiter="\t")
        for row in reader:
            # ignore comment lines
            if "##" not in row[0] and "#CHROM" not in row[0]:
                sequence_name = row[0]
                site = int(row[1]) - 1     # to make this 0 indexed
                reference_allele = row[3].lower()
                alternative_allele = row[4].lower()
 
                # pull out the frequency using a string search
                SearchStr = r'.+\:([0-9]{1,2}\.{0,1}[0-9]{0,2}\%)'
                result = re.search(SearchStr,row[9])
                if result:
                    frequency = '\t'.join(result.groups())
                    frequency_decimal = frequency.replace("%","")
                    frequency_decimal = (float(frequency_decimal))/100
                    frequency_decimal = "%.4f" % frequency_decimal    # only include 4 numbers after decimal

                else:
                    frequency = "none reported"

                # figure out whether the SNP lies within a coding region:
                for gene in coding_regions[sequence_name]:
                    coordinates = coding_regions[sequence_name][gene]

                    # go through gene coordinates 2 at a time; this is for genes with multiple regions
                    for i in range(0,int(len(coordinates)),2):
                        if site >= coordinates[i] and site <= coordinates[i+1]:  # if site is within the gene

                            # determine the coding region site, depending on if there are 2 frames or 1
                            if len(coordinates) == 2:
                                cds_site = site - coordinates[i]
                            elif len(coordinates) == 4:
                                cds_site = (coordinates[1] - coordinates[0]) + (site - coordinates[2] + 1)
                            # now determine whether the site is in the 1st, 2nd, or 3rd codon position
                            # if SNP is in 1st position in codon:
                            aa_site = int(cds_site/3)+1

                            if float(cds_site) % 3 == 0:
                                codon = transcripts[gene][cds_site:cds_site+3]
                                variant_codon = alternative_allele + codon[1:3]
                                variant_aa = Seq(variant_codon).translate()
                                ref_codon = reference_allele + codon[1:3]
                                ref_aa = Seq(ref_codon).translate()

                            # if variant is in the middle of the codon:
                            elif float(cds_site - 1) % 3 == 0:
                                codon = transcripts[gene][cds_site-1:cds_site+2]
                                variant_codon = codon[0] + alternative_allele + codon[2]
                                variant_aa = Seq(variant_codon).translate()
                                ref_codon = codon[0] + reference_allele + codon[2]
                                ref_aa = Seq(ref_codon).translate()

                            # if the variant is in the 3rd codon position
                            elif float(cds_site -2) % 3 == 0:
                                codon = transcripts[gene][cds_site-2:cds_site+1]
                                variant_codon = codon[0:2] + alternative_allele
                                variant_aa = Seq(variant_codon).translate()
                                ref_codon = codon[0:2] + reference_allele
                                ref_aa = Seq(ref_codon).translate()


                            # return the amino acid changes, with single letters converted to 3-letter aa codes
                            if ref_aa == variant_aa:
                                syn_nonsyn = "synonymous"
                            elif ref_aa != variant_aa and variant_aa == "*":
                                syn_nonsyn = "stop_gained"
                            elif ref_aa != variant_aa:
                                syn_nonsyn = "nonsynonymous"

                            amino_acid_change = amino_acid_abbreviations[ref_aa] + str(aa_site) + amino_acid_abbreviations[variant_aa]
                            with open(outfilename, "a") as outfile:
                                output = [sequence_name,gene,str(site+1),reference_allele.upper(),alternative_allele.upper(),amino_acid_change,syn_nonsyn,frequency,frequency_decimal,"\n"]
                                output2 = "\t".join(output)
                                outfile.write(output2)


def coverage_summary(input_tsvs, output_tsv):
    dfs = []
    for tsv_path in input_tsvs:
        split_path = tsv_path.split('/')
        sample_id = split_path[1]
        replicate = split_path[2]
        df = pd.read_csv(tsv_path, sep='\t')
        df['sample_id'] = f'{sample_id}-{replicate}'
        dfs.append(df)
    full_df = pd.concat(dfs, ignore_index=True)
    full_df['total_coverage'] = full_df['0x'] + \
        full_df['1-100x'] + \
        full_df['100-1000x'] + \
        full_df['1000x+']

    full_df.to_csv(output_tsv, sep='\t', index=False)


def mask_low_coverage(fasta_file, coverage_file, threshold, output_fasta):
    # Read the coverage file into a dictionary
    coverage = {}
    with open(coverage_file, 'r') as cf:
        reader = csv.reader(cf, delimiter='\t')
        next(reader)  # Skip header
        for segment, start, end, cov in reader:
            start, end, cov = int(start), int(end), int(cov)
            if segment not in coverage:
                coverage[segment] = []
            coverage[segment].append((start, end, cov))

    # Process FASTA and mask low coverage regions
    seq_records = []
    for seq_record in SeqIO.parse(fasta_file, 'fasta'):
        seq = list(str(seq_record.seq))
        for start, end, cov in coverage.get(seq_record.id, []):
            if cov < threshold:
                for i in range(start, end):
                    if i < len(seq):
                        seq[i] = 'N'
        seq_record.seq = ''.join(seq)
        seq_records.append(seq_record)
    SeqIO.write(seq_records, output_fasta, 'fasta')


if __name__ == '__main__':
    command_line_interface()